import glob
import os
import numpy as np
import torch
from torch.nn.functional import conv3d
import nibabel as nib
from SuperSynth.utils import myzoom_torch
from scipy.io.matlab import loadmat
import ext.interpol as interpol

def super_generator(datadir,
              resolution_sampler,
              label_list_segmentation,
              n_neutral_labels,
              produceMNIregistration=False,
              siz=[160, 160, 160],
              max_rotation=15,
              max_shear=0.2,
              max_scaling=0.2,
              nonlin_scale_min=0.03,
              nonlin_scale_max=0.06,
              nonlin_std_max=4,
              bf_scale_min=0.02,
              bf_scale_max=0.04,
              bf_std_min=0.1,
              bf_std_max=0.6,
              bag_scale_min=0.02,
              bag_scale_max=0.08,
              gamma_std=0.1,
              min_noise_std=5,
              max_noise_std=15,
              exvixo_prob=0.25,
              photo_prob=0.2,
              bag_prob=0.5,
              pv=True,
              deform_one_hots=True,
              integrate_deformation_fields=False,
              produce_surfaces=False,
              bspline_zooming=False,
              device='cpu'):

    # Paths to the different subdirectories
    gen_dir = os.path.join(datadir, 'label_maps_generation')
    seg_dir = os.path.join(datadir, 'label_maps_segmentation')
    dist_dir = os.path.join(datadir, 'Dmaps')
    im_dir = os.path.join(datadir, 'images')
    bag_dir = os.path.join(datadir, 'DmapsBag')
    surface_dir = os.path.join(datadir, 'surfaces')
    mni_dir = os.path.join(datadir, 'MNIreg')
    if os.path.exists(seg_dir) is False:
        print('Directory with target segmentations not found; target segmentations will no be generated')
        seg_dir = None
    if os.path.exists(dist_dir) is False:
        print('Directory with distance maps not found; distance maps will no be generated')
        dist_dir = None
    if os.path.exists(im_dir) is False:
        print('Directory with real images not found; real images will no be generated')
        im_dir = None
    if os.path.exists(bag_dir) is False:
        print('Directory with distance maps for bag simulation not found; fake bags will no be generated')
        bag_dir = None
    if produce_surfaces is False:
        print('Surface generation switched off by user')
        surface_dir = None
    else:
        if os.path.exists(surface_dir) is False:
            raise Exception('User is asking for surfaces but directory with surface files was not found!')
        if integrate_deformation_fields is False:
            raise Exception('Using surfaces requires integrating deformation fields; you need to switch on integrate_deformation_fields option')
    if produceMNIregistration is False:
        print('MNI registration switched off by user')
        mni_dir = None
    else:
        if os.path.exists(mni_dir) is False:
            raise Exception('User is asking for MNI registrations surfaces but directory with deformation files was not found!')

    names = glob.glob(os.path.join(gen_dir, '*.nii.gz')) + glob.glob(os.path.join(gen_dir, '*.nii'))
    n = len(names)

    # Get resolution of training data
    aff = nib.load(names[0]).affine
    res_training_data = np.sqrt(np.sum(aff[:-1, :-1], axis=0))
    n_steps_svf_integration = 8

    with torch.no_grad():
        # prepare grid
        print('Preparing grid...')
        xx, yy, zz = np.meshgrid(range(siz[0]), range(siz[1]), range(siz[2]), sparse=False, indexing='ij')
        xx = torch.tensor(xx, dtype=torch.float, device=device)
        yy = torch.tensor(yy, dtype=torch.float, device=device)
        zz = torch.tensor(zz, dtype=torch.float, device=device)
        c = torch.tensor((np.array(siz) - 1) / 2, dtype=torch.float, device=device)
        xc = xx - c[0]
        yc = yy - c[1]
        zc = zz - c[2]

        # Matrix for one-hot encoding (includes a lookup-table)
        n_labels = len(label_list_segmentation)
        lut = torch.zeros(10000, dtype=torch.long, device=device)
        for l in range(n_labels):
            lut[label_list_segmentation[l]] = l
        onehotmatrix = torch.eye(n_labels, dtype=torch.float, device=device)

        nlat = int((n_labels - n_neutral_labels) / 2.0)
        vflip = np.concatenate([np.array(range(n_neutral_labels)),
                                np.array(range(n_neutral_labels + nlat, n_labels)),
                                np.array(range(n_neutral_labels, n_neutral_labels + nlat))])

        print('Generator is ready!')


        while True:

            # Select random case
            idx = np.random.randint(n)
            photo_mode = np.random.rand()<photo_prob
            Gimg = nib.load(names[idx])

            # The first thing we do is sampling the resolution and deformation, as this will give us a bounding box
            # of the image region we need, so we don't have to read the whole thing from disk (only works for uncompressed niftis!

            # Sample resolution
            if photo_mode:
                spac = 2.0 + 10 * np.random.rand()
                resolution = np.array([res_training_data[0], spac , res_training_data[2]])
                thickness = np.array([res_training_data[0], 0.0001, res_training_data[2]])
            else:
                resolution, thickness = resolution_sampler()

            # sample affine deformation
            rotations = (2 * max_rotation * np.random.rand(3) - max_rotation) / 180.0 * np.pi
            shears = (2 * max_shear * np.random.rand(3) - max_shear)
            scalings = 1 + (2 * max_scaling * np.random.rand(3) - max_scaling)
            scaling_factor_distances = np.prod(scalings) ** .33333333333 # we divide distance maps by this, not perfect, but better than nothing
            A = torch.tensor(make_affine_matrix(rotations, shears, scalings), dtype=torch.float, device=device)

            # sample center
            max_shift = (torch.tensor(np.array(Gimg.shape[0:3]) - siz, dtype=torch.float, device=device)) / 2
            max_shift[max_shift < 0] = 0
            c2 = torch.tensor((np.array(Gimg.shape[0:3]) - 1)/2, dtype=torch.float, device=device) + (2 * (max_shift * torch.rand(3, dtype=float, device=device)) - max_shift)

            # sample nonlinear deformation
            nonlin_scale = nonlin_scale_min + np.random.rand(1) * (nonlin_scale_max - nonlin_scale_min)
            siz_F_small = np.round(nonlin_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_F_small[1] = np.round(siz[1]/spac).astype(int)
            nonlin_std = nonlin_std_max * np.random.rand()
            Fsmall = nonlin_std * torch.randn([*siz_F_small, 3], dtype=torch.float, device=device)
            F = myzoom_torch(Fsmall, np.array(siz) / siz_F_small, device)
            if photo_mode:
                F[:, :, :, 1] = 0

            if integrate_deformation_fields:
                steplength = 1.0 / (2.0 ** n_steps_svf_integration)
                Fsvf = F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf += fast_3D_interp_torch(Fsvf, xx + Fsvf[:, :, :, 0], yy + Fsvf[:, :, :, 1], zz + Fsvf[:, :, :, 2], 'linear', device)
                Fsvf_neg = -F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf_neg += fast_3D_interp_torch(Fsvf_neg, xx + Fsvf_neg[:, :, :, 0], yy + Fsvf_neg[:, :, :, 1], zz + Fsvf_neg[:, :, :, 2], 'linear', device)
                F = Fsvf
                Fneg = Fsvf_neg

            # Start by deforming surfaces if needed (we need the inverse transform!)
            if produce_surfaces:
                filename = os.path.basename(names[idx])
                if filename.endswith('.nii.gz'):
                    filename = filename[:-7] + '.mat'
                else:
                    filename = filename[:-4] + '.mat'
                mat = loadmat(os.path.join(surface_dir, filename ))
                Vlw = torch.tensor(mat['Vlw'], dtype=torch.float, device=device)
                Flw = torch.tensor(mat['Flw'], dtype=torch.int, device=device)
                Vrw = torch.tensor(mat['Vrw'], dtype=torch.float, device=device)
                Frw = torch.tensor(mat['Frw'], dtype=torch.int, device=device)
                Vlp = torch.tensor(mat['Vlp'], dtype=torch.float, device=device)
                Flp = torch.tensor(mat['Flp'], dtype=torch.int, device=device)
                Vrp = torch.tensor(mat['Vrp'], dtype=torch.float, device=device)
                Frp = torch.tensor(mat['Frp'], dtype=torch.int, device=device)

                Ainv = torch.inverse(A);
                Vlw -= c2[None, :]
                Vlw = Vlw @ torch.transpose(Ainv, 0, 1)
                Vlw += fast_3D_interp_torch(Fneg, Vlw[:, 0]+c[0], Vlw[:, 1]+c[1], Vlw[:, 2]+c[2], 'linear', device)
                Vlw += c[None, :]
                Vrw -= c2[None, :]
                Vrw = Vrw @ torch.transpose(Ainv, 0, 1)
                Vrw += fast_3D_interp_torch(Fneg, Vrw[:, 0]+c[0], Vrw[:, 1]+c[1], Vrw[:, 2]+c[2], 'linear', device)
                Vrw += c[None, :]
                Vlp -= c2[None, :]
                Vlp = Vlp @ torch.transpose(Ainv, 0, 1)
                Vlp += fast_3D_interp_torch(Fneg, Vlp[:, 0] + c[0], Vlp[:, 1] + c[1], Vlp[:, 2] + c[2], 'linear', device)
                Vlp += c[None, :]
                Vrp -= c2[None, :]
                Vrp = Vrp @ torch.transpose(Ainv, 0, 1)
                Vrp += fast_3D_interp_torch(Fneg, Vrp[:, 0] + c[0], Vrp[:, 1] + c[1], Vrp[:, 2] + c[2], 'linear', device)
                Vrp += c[None, :]


            # deform the images (we do nonlinear "first" ie after so we can do heavy coronal deformations in photo mode)
            xx1 = xc + F[:, :, :, 0]
            yy1 = yc + F[:, :, :, 1]
            zz1 = zc + F[:, :, :, 2]
            xx2 = A[0, 0] * xx1 + A[0, 1] * yy1 + A[0, 2] * zz1 + c2[0]
            yy2 = A[1, 0] * xx1 + A[1, 1] * yy1 + A[1, 2] * zz1 + c2[1]
            zz2 = A[2, 0] * xx1 + A[2, 1] * yy1 + A[2, 2] * zz1 + c2[2]
            xx2[xx2 < 0] = 0
            yy2[yy2 < 0] = 0
            zz2[zz2 < 0] = 0
            xx2[xx2 > (Gimg.shape[0] - 1)] = Gimg.shape[0] - 1
            yy2[yy2 > (Gimg.shape[1] - 1)] = Gimg.shape[1] - 1
            zz2[zz2 > (Gimg.shape[2] - 1)] = Gimg.shape[2] - 1

            # Get the margins for reading images
            x1 = torch.floor(torch.min(xx2))
            y1 = torch.floor(torch.min(yy2))
            z1 = torch.floor(torch.min(zz2))
            x2 = 1+torch.ceil(torch.max(xx2))
            y2 = 1 + torch.ceil(torch.max(yy2))
            z2 = 1 + torch.ceil(torch.max(zz2))
            xx2 -= x1
            yy2 -= y1
            zz2 -= z1

            x1 = x1.cpu().numpy().astype(int)
            y1 = y1.cpu().numpy().astype(int)
            z1 = z1.cpu().numpy().astype(int)
            x2 = x2.cpu().numpy().astype(int)
            y2 = y2.cpu().numpy().astype(int)
            z2 = z2.cpu().numpy().astype(int)


            # Read in data
            G = torch.squeeze(torch.tensor(Gimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
            MNI = S = D = I = B = None
            if seg_dir is not None:
                Simg = nib.load(os.path.join(seg_dir, os.path.basename(names[idx])))
                S = torch.squeeze(torch.tensor(Simg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(int), dtype=torch.int, device=device))
            if dist_dir is not None:
                Dimg = nib.load(os.path.join(dist_dir, os.path.basename(names[idx])))
                D = torch.squeeze(torch.tensor(Dimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
                D /= scaling_factor_distances
            if im_dir is not None:
                Iimg = nib.load(os.path.join(im_dir, os.path.basename(names[idx])))
                I = torch.squeeze(torch.tensor(Iimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
                I[I < 0] = 0
                I /= torch.median(I[G==2])
            if mni_dir is not None:
                MNIimg = nib.load(os.path.join(mni_dir, os.path.basename(names[idx])))
                MNI = torch.squeeze(torch.tensor(MNIimg.get_fdata()[x1:x2, y1:y2, z1:z2, :].astype(float), dtype=torch.float, device=device))
                MNI /= 100
            if bag_dir is not None:
                Bimg = nib.load(os.path.join(bag_dir, os.path.basename(names[idx])))
                B = torch.squeeze(torch.tensor(Bimg.get_fdata()[x1:x2, y1:y2, z1:z2], dtype=torch.float, device=device))
                B /= scaling_factor_distances

            # Decide if we're simulating ex vivo (and possibly a bag) or photos
            if photo_mode or (np.random.rand() < exvixo_prob):
                G[G>255] = 0 # kill extracerebral
                if photo_mode:
                    G[G == 7] = 0
                    G[G == 8] = 0
                    G[G == 16] = 0
                    S[S == 24] = 0
                    S[S == 7] = 0
                    S[S == 8] = 0
                    S[S == 46] = 0
                    S[S == 47] = 0
                    S[S == 15] = 0
                    S[S == 16] = 0
                    if D is None: # without distance maps, killing 4 is the best we can do
                        G[G == 4] = 0
                    else:
                        Dpial = torch.minimum(D[...,1], D[..., 3])
                        th = 1.5 * np.random.rand() # band of random width...
                        G[G==4] = 0
                        G[(G == 0) & (Dpial < th)] = 4

                elif ((B is not None) and (np.random.rand(1) < bag_prob)):
                    bag_scale = bag_scale_min + np.random.rand(1) * (bag_scale_max - bag_scale_min)
                    siz_TH_small = np.round(bag_scale * np.array(G.shape)).astype(int).tolist()
                    bag_tness = torch.tensor(np.sort(1.0 + 20 * np.random.rand(2)), dtype=torch.float, device=device)
                    THsmall = bag_tness[0] + (bag_tness[1] - bag_tness[0]) * torch.rand(siz_TH_small, dtype=torch.float, device=device)
                    TH = myzoom_torch(THsmall, np.array(G.shape) / siz_TH_small, device)
                    G[(B>0) & (B<TH)] = 4

            # Sample Gaussian image
            mus = 25 + 200 * torch.rand(10000, dtype=torch.float, device=device)
            sigmas = 5 + 20 * torch.rand(10000, dtype=torch.float, device=device)
            if photo_mode or np.random.rand(1)<0.5: # set the background to zero every once in a while (or always in photo mode)
                mus[0] = 0
            Gr = torch.round(G).long()
            SYN = mus[Gr] + sigmas[Gr] * torch.randn(Gr.shape, dtype=torch.float, device=device)
            if pv:
                mask = (G!=Gr)
                SYN[mask] = 0
                Gv = G[mask]
                isv = torch.zeros(Gv.shape, dtype=torch.float, device=device )
                pw = (Gv<=3) * (3-Gv)
                isv += pw * mus[2] + pw * sigmas[2] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pg = (Gv<=3) * (Gv-2) + (Gv>3) * (4-Gv)
                isv += pg * mus[3] + pg * sigmas[3] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pcsf = (Gv>=3) * (Gv-3)
                isv += pcsf * mus[4] + pcsf * sigmas[4] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                SYN[mask] = isv
            SYN[SYN < 0] = 0


            SYNdef = fast_3D_interp_torch(SYN, xx2, yy2, zz2, 'linear', device)
            SdefOneHot = Ddef = Idef = None
            if S is not None:
                if deform_one_hots:
                    Sonehot = onehotmatrix[lut[S.long()]]
                    SdefOneHot = fast_3D_interp_torch(Sonehot, xx2, yy2, zz2, 'linear', device)
                else:
                    Sdef = fast_3D_interp_torch(S, xx2, yy2, zz2, 'nearest', device)
                    SdefOneHot = onehotmatrix[lut[Sdef.long()]]

            if D is not None:
                Ddef = fast_3D_interp_torch(D, xx2, yy2, zz2, 'linear', device, default_value_linear=torch.max(D))
            if I is not None:
                Idef = fast_3D_interp_torch(I, xx2, yy2, zz2, 'linear', device)
            if MNI is not None:
                MNIdef = fast_3D_interp_torch(MNI, xx2, yy2, zz2, 'linear', device)
            else:
                MNIdef = None

            # Gamma transform
            gamma = torch.tensor(np.exp(gamma_std * np.random.randn(1)[0]), dtype=float, device=device)
            SYNgamma = 300.0 * (SYNdef / 300.0) ** gamma

            # Bias field
            bf_scale = bf_scale_min + np.random.rand(1) * (bf_scale_max - bf_scale_min)
            siz_BF_small = np.round(bf_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_BF_small[1] = np.round(siz[1]/spac).astype(int)
            BFsmall = torch.tensor(bf_std_min + (bf_std_max - bf_std_min) * np.random.rand(1), dtype=torch.float, device=device) * torch.randn(siz_BF_small, dtype=torch.float, device=device)
            BFlog = myzoom_torch(BFsmall, np.array(siz) / siz_BF_small, device)
            BF = torch.exp(BFlog)
            SYNbf = SYNgamma * BF

            # Model Resolution
            stds = (0.85 + 0.3 * np.random.rand()) * np.log(5) /np.pi * thickness / res_training_data
            stds[thickness<=res_training_data] = 0.0 # no blur if thickness is equal to the resolution of the training data
            SYNblur = gaussian_blur_3d(SYNbf, stds, device)
            newsize = (np.array(siz) * res_training_data / resolution).astype(int)

            factors = np.array(newsize) / np.array(siz)
            delta = (1.0 - factors) / (2.0 * factors)
            vx = np.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0])[:newsize[0]]
            vy = np.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1])[:newsize[1]]
            vz = np.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2])[:newsize[2]]
            II, JJ, KK = np.meshgrid(vx, vy, vz, sparse=False, indexing='ij')
            II = torch.tensor(II, dtype=torch.float, device=device)
            JJ = torch.tensor(JJ, dtype=torch.float, device=device)
            KK = torch.tensor(KK, dtype=torch.float, device=device)

            SYNsmall = fast_3D_interp_torch(SYNblur, II, JJ, KK, 'linear', device)
            noise_std = torch.tensor(min_noise_std + (max_noise_std - min_noise_std) * np.random.rand(1), dtype=torch.float, device=device)
            SYNnoisy = SYNsmall + noise_std * torch.randn(SYNsmall.shape, dtype=torch.float, device=device)
            SYNnoisy[SYNnoisy<0] = 0

            # Back to original resolution
            if bspline_zooming:
                SYNresized = interpol.resize(SYNnoisy, shape=siz, anchor='edge', interpolation=3, bound='dct2', prefilter=True)
            else:
                SYNresized = myzoom_torch(SYNnoisy, 1 / factors, device)
            maxi = torch.max(SYNresized)
            SYNfinal = SYNresized / maxi

            # Flip 50% of times
            if np.random.rand()<0.5:
                SYNfinal = torch.flip(SYNfinal, [0])
                SdefOneHot = torch.flip(SdefOneHot, [0])[:, :, :, vflip]
                Ddef = torch.flip(Ddef, [0])[:, :, :, [2,3,0,1]]
                Idef = torch.flip(Idef, [0])
                if MNIdef is not None:
                    MNIdef = torch.flip(MNIdef, [0])
                    MNIdef[:, :, :, 0] = -MNIdef[:, :, :, 0] # pretty easy thanks to symmetric template
                BFlog = torch.flip(BFlog, [0])
                if produce_surfaces:
                    Vlw[:, 0] = Idef.shape[0] - 1 - Vlw[:, 0]
                    Vrw[:, 0] = Idef.shape[0] - 1 - Vrw[:, 0]
                    Vlp[:, 0] = Idef.shape[0] - 1 - Vlp[:, 0]
                    Vrp[:, 0] = Idef.shape[0] - 1 - Vrp[:, 0]
                    Vlw, Vrw = Vrw, Vlw
                    Vlp, Vrp = Vrp, Vlp
                    Flw, Frw = Frw, Flw
                    Flp, Frp = Frp, Flp

            # mask real image and MNI coordinates if needed
            if Idef is not None:
                Idef *= (1.0 - SdefOneHot[:, :, :, 0])
            if MNIdef is not None:
                MNIdef *= (1.0 - SdefOneHot[:, :, :, 0, None])

            if produceMNIregistration:
                if produce_surfaces:
                    yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog, MNIdef, Vlw, Flw, Vrw, Frw, Vlp, Flp, Vrp, Frp]
                else:
                    yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog, MNIdef]
            else:
                if produce_surfaces:
                    yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog, Vlw, Flw, Vrw, Frw, Vlp, Flp, Vrp, Frp]
                else:
                    yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog]


def supervised_generator(datadir,
              label_list_segmentation,
              n_neutral_labels,
              siz=[160, 160, 160],
              max_rotation=15,
              max_shear=0.2,
              max_scaling=0.2,
              nonlin_scale_min=0.03,
              nonlin_scale_max=0.06,
              nonlin_std_max=4,
              min_noise_std=0, # some noise is ok, but not too much, images already have noise
              max_noise_std=0.02,
              photo_prob=0.2,
              deform_one_hots=True,
              device='cpu'):

    # Paths to the different subdirectories
    seg_dir = os.path.join(datadir, 'label_maps_segmentation')
    im_dir = os.path.join(datadir, 'images')
    if os.path.exists(seg_dir) is False:
        raise Exception('Directory with target segmentations not found; target segmentations will no be generated')
        seg_dir = None
    if os.path.exists(im_dir) is False:
        print('Directory with real images not found; real images will no be generated')
        im_dir = None

    names = glob.glob(os.path.join(seg_dir, '*.nii.gz')) + glob.glob(os.path.join(seg_dir, '*.nii'))
    n = len(names)

    with torch.no_grad():
        # prepare grid
        print('Preparing grid...')
        xx, yy, zz = np.meshgrid(range(siz[0]), range(siz[1]), range(siz[2]), sparse=False, indexing='ij')
        c = torch.tensor((np.array(siz) - 1) / 2, dtype=torch.float, device=device)
        xc = torch.tensor(xx, dtype=torch.float, device=device) - c[0]
        yc = torch.tensor(yy, dtype=torch.float, device=device) - c[1]
        zc = torch.tensor(zz, dtype=torch.float, device=device) - c[2]

        # Matrix for one-hot encoding (includes a lookup-table)
        n_labels = len(label_list_segmentation)
        lut = torch.zeros(10000, dtype=torch.long, device=device)
        for l in range(n_labels):
            lut[label_list_segmentation[l]] = l
        onehotmatrix = torch.eye(n_labels, dtype=torch.float, device=device)

        nlat = int((n_labels - n_neutral_labels) / 2.0)
        vflip = np.concatenate([np.array(range(n_neutral_labels)),
                                np.array(range(n_neutral_labels + nlat, n_labels)),
                                np.array(range(n_neutral_labels, n_neutral_labels + nlat))])

        print('Generator is ready!')


        while True:

            # Select random case
            idx = np.random.randint(n)
            Simg = nib.load(names[idx])
            photo_mode = np.random.rand() < photo_prob

            # The first thing we do is sampling the  deformation, as this will give us a bounding box
            # of the image region we need, so we don't have to read the whole thing from disk
            # (only works for uncompressed niftis!)

            # sample affine deformation
            rotations = (2 * max_rotation * np.random.rand(3) - max_rotation) / 180.0 * np.pi
            shears = (2 * max_shear * np.random.rand(3) - max_shear)
            scalings = 1 + (2 * max_scaling * np.random.rand(3) - max_scaling)
            A = torch.tensor(make_affine_matrix(rotations, shears, scalings), dtype=torch.float, device=device)

            # sample nonlinear deformation
            nonlin_scale = nonlin_scale_min + np.random.rand(1) * (nonlin_scale_max - nonlin_scale_min)
            siz_F_small = np.round(nonlin_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                spac = 2.0 + 10 * np.random.rand()
                siz_F_small[1] = np.round(siz[1]/spac).astype(int)
            nonlin_std = nonlin_std_max * np.random.rand()
            Fsmall = nonlin_std * torch.randn([*siz_F_small, 3], dtype=torch.float, device=device)
            F = myzoom_torch(Fsmall, np.array(siz) / siz_F_small, device)
            if photo_mode:
                F[:, :, :, 1] = 0

            # sample center
            max_shift = (torch.tensor(np.array(Simg.shape[0:3]) - siz, dtype=torch.float, device=device)) / 2
            max_shift[max_shift < 0] = 0
            c2 = torch.tensor((np.array(Simg.shape[0:3]) - 1)/2, dtype=torch.float, device=device) + (2 * (max_shift * torch.rand(3, dtype=float, device=device)) - max_shift)

            # deform (we do nonlinear "first" ie after so we can do heavy coronal deformations in photo mode)
            xx1 = xc + F[:, :, :, 0]
            yy1 = yc + F[:, :, :, 1]
            zz1 = zc + F[:, :, :, 2]
            xx2 = A[0, 0] * xx1 + A[0, 1] * yy1 + A[0, 2] * zz1 + c2[0]
            yy2 = A[1, 0] * xx1 + A[1, 1] * yy1 + A[1, 2] * zz1 + c2[1]
            zz2 = A[2, 0] * xx1 + A[2, 1] * yy1 + A[2, 2] * zz1 + c2[2]
            xx2[xx2 < 0] = 0
            yy2[yy2 < 0] = 0
            zz2[zz2 < 0] = 0
            xx2[xx2 > (Simg.shape[0] - 1)] = Simg.shape[0] - 1
            yy2[yy2 > (Simg.shape[1] - 1)] = Simg.shape[1] - 1
            zz2[zz2 > (Simg.shape[2] - 1)] = Simg.shape[2] - 1

            # Get the margins for reading images
            x1 = torch.floor(torch.min(xx2))
            y1 = torch.floor(torch.min(yy2))
            z1 = torch.floor(torch.min(zz2))
            x2 = 1+torch.ceil(torch.max(xx2))
            y2 = 1 + torch.ceil(torch.max(yy2))
            z2 = 1 + torch.ceil(torch.max(zz2))
            xx2 -= x1
            yy2 -= y1
            zz2 -= z1

            x1 = x1.cpu().numpy().astype(int)
            y1 = y1.cpu().numpy().astype(int)
            z1 = z1.cpu().numpy().astype(int)
            x2 = x2.cpu().numpy().astype(int)
            y2 = y2.cpu().numpy().astype(int)
            z2 = z2.cpu().numpy().astype(int)

            # Read in data
            S = torch.squeeze(torch.tensor(Simg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
            S[S == 24] = 0
            Iimg = nib.load(os.path.join(im_dir, os.path.basename(names[idx])))
            I = torch.squeeze(torch.tensor(Iimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))

            # Deform
            Idef = fast_3D_interp_torch(I, xx2, yy2, zz2, 'linear', device)
            if deform_one_hots:
                Sonehot = onehotmatrix[lut[S.long()]]
                SdefOneHot = fast_3D_interp_torch(Sonehot, xx2, yy2, zz2, 'linear', device)
            else:
                Sdef = fast_3D_interp_torch(S, xx2, yy2, zz2, 'nearest', device)
                SdefOneHot = onehotmatrix[lut[Sdef.long()]]

            # Noise
            noise_std = torch.tensor(min_noise_std + (max_noise_std - min_noise_std) * np.random.rand(1), dtype=torch.float, device=device)
            Inoisy = torch.clamp(Idef + noise_std * torch.randn(Idef.shape, dtype=torch.float, device=device), min=0)

            # Flip 50% of times
            if np.random.rand()<0.5:
                Inoisy = torch.flip(Inoisy, [0])
                SdefOneHot = torch.flip(SdefOneHot, [0])[:, :, :, vflip]

            Inoisy *= (1.0 - SdefOneHot[:, :, :, 0])

            yield [Inoisy, SdefOneHot]


###############################
# SINGLE HEMISPHERE GENERATOR #
###############################

def super_generator_hemi(datadir,
              resolution_sampler,
              label_list_segmentation,
              siz=[96, 128, 128],
              max_rotation=15,
              max_shear=0.2,
              max_scaling=0.2,
              nonlin_scale_min=0.03,
              nonlin_scale_max=0.06,
              nonlin_std_max=4,
              bf_scale_min=0.02,
              bf_scale_max=0.04,
              bf_std_min=0.1,
              bf_std_max=0.6,
              bag_scale_min=0.02,
              bag_scale_max=0.08,
              bag_probability=0.5,
              gamma_std=0.1,
              min_noise_std=5,
              max_noise_std=15,
              exvixo_prob_vs_photo=0.66666666666666,
              pv=True,
              deform_one_hots=True,
              integrate_deformation_fields=False,
              produce_surfaces=False,
              bspline_zooming=False,
              device='cpu'):

    # Paths to the different subdirectories
    gen_dir = os.path.join(datadir, 'label_maps_generation')
    seg_dir = os.path.join(datadir, 'label_maps_segmentation')
    dist_dir = os.path.join(datadir, 'Dmaps')
    im_dir = os.path.join(datadir, 'images')
    bag_dir = os.path.join(datadir, 'DmapsBag')
    surface_dir = os.path.join(datadir, 'surfaces')
    if os.path.exists(seg_dir) is False:
        print('Directory with target segmentations not found; target segmentations will no be generated')
        seg_dir = None
    if os.path.exists(dist_dir) is False:
        print('Directory with distance maps not found; distance maps will no be generated')
        dist_dir = None
    if os.path.exists(im_dir) is False:
        print('Directory with real images not found; real images will no be generated')
        im_dir = None
    if os.path.exists(bag_dir) is False:
        print('Directory with distance maps for bag simulation not found; fake bags will no be generated')
        bag_dir = None
    if produce_surfaces is False:
        print('Surface generation switched off by user')
        surface_dir = None
    else:
        if os.path.exists(surface_dir) is False:
            raise Exception('User is asking for surfaces but directory with surface files was not found!')
        if integrate_deformation_fields is False:
            raise Exception('Using surfaces requires integrating deformation fields; you need to switch on integrate_deformation_fields option')

    # find an empty label we can use to generate a bag
    auxiliary_label = 0
    while auxiliary_label in label_list_segmentation:
        auxiliary_label += 1

    names = glob.glob(os.path.join(gen_dir, '*.nii.gz')) + glob.glob(os.path.join(gen_dir, '*.nii'))
    n = len(names)

    # Get resolution of training data
    aff = nib.load(names[0]).affine
    res_training_data = np.sqrt(np.sum(aff[:-1, :-1], axis=0))
    n_steps_svf_integration = 8

    with torch.no_grad():
        # prepare grid
        print('Preparing grid...')
        xx, yy, zz = np.meshgrid(range(siz[0]), range(siz[1]), range(siz[2]), sparse=False, indexing='ij')
        xx = torch.tensor(xx, dtype=torch.float, device=device)
        yy = torch.tensor(yy, dtype=torch.float, device=device)
        zz = torch.tensor(zz, dtype=torch.float, device=device)
        c = torch.tensor((np.array(siz) - 1) / 2, dtype=torch.float, device=device)
        xc = xx - c[0]
        yc = yy - c[1]
        zc = zz - c[2]

        # Matrix for one-hot encoding (includes a lookup-table)
        n_labels = len(label_list_segmentation)
        lut = torch.zeros(10000, dtype=torch.long, device=device)
        for l in range(n_labels):
            lut[label_list_segmentation[l]] = l
        onehotmatrix = torch.eye(n_labels, dtype=torch.float, device=device)

        print('Generator is ready!')


        while True:

            # Select random case
            idx = np.random.randint(n)
            exvivo_mode = np.random.rand() < exvixo_prob_vs_photo
            photo_mode = (not exvivo_mode)
            Gimg = nib.load(names[idx])

            # The first thing we do is sampling the resolution and deformation, as this will give us a bounding box
            # of the image region we need, so we don't have to read the whole thing from disk (only works for uncompressed niftis!

            # Sample resolution
            if photo_mode:
                spac = 2.0 + 10 * np.random.rand()
                resolution = np.array([res_training_data[0], spac , res_training_data[2]])
                thickness = np.array([res_training_data[0], 0.0001, res_training_data[2]])
            else:
                resolution, thickness = resolution_sampler()

            # sample affine deformation
            rotations = (2 * max_rotation * np.random.rand(3) - max_rotation) / 180.0 * np.pi
            shears = (2 * max_shear * np.random.rand(3) - max_shear)
            scalings = 1 + (2 * max_scaling * np.random.rand(3) - max_scaling)
            scaling_factor_distances = np.prod(scalings) ** .33333333333 # we divide distance maps by this, not perfect, but better than nothing
            A = torch.tensor(make_affine_matrix(rotations, shears, scalings), dtype=torch.float, device=device)

            # sample center
            max_shift = (torch.tensor(np.array(Gimg.shape[0:3]) - siz, dtype=torch.float, device=device)) / 2
            max_shift[max_shift<0] = 0
            c2 = torch.tensor((np.array(Gimg.shape[0:3]) - 1)/2, dtype=torch.float, device=device) + (2 * (max_shift * torch.rand(3, dtype=float, device=device)) - max_shift)

            # sample nonlinear deformation
            nonlin_scale = nonlin_scale_min + np.random.rand(1) * (nonlin_scale_max - nonlin_scale_min)
            siz_F_small = np.round(nonlin_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_F_small[1] = np.round(siz[1]/spac).astype(int)
            nonlin_std = nonlin_std_max * np.random.rand()
            Fsmall = nonlin_std * torch.randn([*siz_F_small, 3], dtype=torch.float, device=device)
            F = myzoom_torch(Fsmall, np.array(siz) / siz_F_small, device)
            if photo_mode:
                F[:, :, :, 1] = 0

            if integrate_deformation_fields:
                steplength = 1.0 / (2.0 ** n_steps_svf_integration)
                Fsvf = F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf += fast_3D_interp_torch(Fsvf, xx + Fsvf[:, :, :, 0], yy + Fsvf[:, :, :, 1], zz + Fsvf[:, :, :, 2], 'linear', device)
                Fsvf_neg = -F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf_neg += fast_3D_interp_torch(Fsvf_neg, xx + Fsvf_neg[:, :, :, 0], yy + Fsvf_neg[:, :, :, 1], zz + Fsvf_neg[:, :, :, 2], 'linear', device)
                F = Fsvf
                Fneg = Fsvf_neg

            # Start by deforming surfaces if needed (we need the inverse transform!)
            if produce_surfaces:
                filename = os.path.basename(names[idx])
                if filename.endswith('.nii.gz'):
                    filename = filename[:-7] + '.mat'
                else:
                    filename = filename[:-4] + '.mat'
                mat = loadmat(os.path.join(surface_dir, filename ))
                Vlw = torch.tensor(mat['Vlw'], dtype=torch.float, device=device)
                Flw = torch.tensor(mat['Flw'], dtype=torch.int, device=device)
                Vlp = torch.tensor(mat['Vlp'], dtype=torch.float, device=device)
                Flp = torch.tensor(mat['Flp'], dtype=torch.int, device=device)

                Ainv = torch.inverse(A);
                Vlw -= c2[None, :]
                Vlw = Vlw @ torch.transpose(Ainv, 0, 1)
                Vlw += fast_3D_interp_torch(Fneg, Vlw[:, 0]+c[0], Vlw[:, 1]+c[1], Vlw[:, 2]+c[2], 'linear', device)
                Vlw += c[None, :]
                Vlp -= c2[None, :]
                Vlp = Vlp @ torch.transpose(Ainv, 0, 1)
                Vlp += fast_3D_interp_torch(Fneg, Vlp[:, 0] + c[0], Vlp[:, 1] + c[1], Vlp[:, 2] + c[2], 'linear', device)
                Vlp += c[None, :]


            # deform the images (we do nonlinear "first" ie after so we can do heavy coronal deformations in photo mode)
            xx1 = xc + F[:, :, :, 0]
            yy1 = yc + F[:, :, :, 1]
            zz1 = zc + F[:, :, :, 2]
            xx2 = A[0, 0] * xx1 + A[0, 1] * yy1 + A[0, 2] * zz1 + c2[0]
            yy2 = A[1, 0] * xx1 + A[1, 1] * yy1 + A[1, 2] * zz1 + c2[1]
            zz2 = A[2, 0] * xx1 + A[2, 1] * yy1 + A[2, 2] * zz1 + c2[2]
            xx2[xx2 < 0] = 0
            yy2[yy2 < 0] = 0
            zz2[zz2 < 0] = 0
            xx2[xx2 > (Gimg.shape[0] - 1)] = Gimg.shape[0] - 1
            yy2[yy2 > (Gimg.shape[1] - 1)] = Gimg.shape[1] - 1
            zz2[zz2 > (Gimg.shape[2] - 1)] = Gimg.shape[2] - 1

            # Get the margins for reading images
            x1 = torch.floor(torch.min(xx2))
            y1 = torch.floor(torch.min(yy2))
            z1 = torch.floor(torch.min(zz2))
            x2 = 1+torch.ceil(torch.max(xx2))
            y2 = 1 + torch.ceil(torch.max(yy2))
            z2 = 1 + torch.ceil(torch.max(zz2))
            xx2 -= x1
            yy2 -= y1
            zz2 -= z1

            x1 = x1.cpu().numpy().astype(int)
            y1 = y1.cpu().numpy().astype(int)
            z1 = z1.cpu().numpy().astype(int)
            x2 = x2.cpu().numpy().astype(int)
            y2 = y2.cpu().numpy().astype(int)
            z2 = z2.cpu().numpy().astype(int)


            # Read in data
            G = torch.squeeze(torch.tensor(Gimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
            S = D = I = B = None
            if seg_dir is not None:
                Simg = nib.load(os.path.join(seg_dir, os.path.basename(names[idx])))
                S = torch.squeeze(torch.tensor(Simg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(int), dtype=torch.int, device=device))
            if dist_dir is not None:
                Dimg = nib.load(os.path.join(dist_dir, os.path.basename(names[idx])))
                D = torch.squeeze(torch.tensor(Dimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
                D /= scaling_factor_distances
            if im_dir is not None:
                Iimg = nib.load(os.path.join(im_dir, os.path.basename(names[idx])))
                I = torch.squeeze(torch.tensor(Iimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
                I[I < 0] = 0
                I /= torch.median(I[G==2])
            if bag_dir is not None:
                Bimg = nib.load(os.path.join(bag_dir, os.path.basename(names[idx])))
                B = torch.squeeze(torch.tensor(Bimg.get_fdata()[x1:x2, y1:y2, z1:z2], dtype=torch.float, device=device))
                B /= scaling_factor_distances

            # If we have a distance transform, we simulate a thin layer of noise (photo) or a thicker bag (ex vivo)
            bg_label = 0
            if (B is not None):
                if exvivo_mode:
                    if np.random.rand() < bag_probability:
                        bag_tness = torch.tensor(np.sort(1.0 + 20 * np.random.rand(2)), dtype=torch.float, device=device)
                    else:
                        G[G==0] = 4 # change background label to 4 for proper transition / partial voluming
                        bg_label = 4
                    if np.random.rand() < 0.5:  # every once in a while, create contrast between ventricle and bag
                        G[S == 4] = auxiliary_label
                else: # photo mode
                    bag_tness = torch.tensor(np.sort(2.0 * np.random.rand(2)), dtype=torch.float, device=device)
                    if np.random.rand() < 0.5: # every once in a while, set the ventricles to background
                        G[S==4] = 0
                if bg_label==0:
                    bag_scale = bag_scale_min + np.random.rand(1) * (bag_scale_max - bag_scale_min)
                    siz_TH_small = np.round(bag_scale * np.array(G.shape)).astype(int).tolist()
                    THsmall = bag_tness[0] + (bag_tness[1] - bag_tness[0]) * torch.rand(siz_TH_small, dtype=torch.float, device=device)
                    TH = myzoom_torch(THsmall, np.array(G.shape) / siz_TH_small, device)
                    G[(B > 0) & (G == 0) & (B < TH)] = 4

            # Sample Gaussian image
            mus = 25 + 200 * torch.rand(10000, dtype=torch.float, device=device)
            sigmas = 5 + 20 * torch.rand(10000, dtype=torch.float, device=device)
            if photo_mode or np.random.rand(1)<0.5: # set the background to zero every once in a while (or always in photo mode)
                mus[bg_label] = 0
            Gr = torch.round(G).long()
            SYN = mus[Gr] + sigmas[Gr] * torch.randn(Gr.shape, dtype=torch.float, device=device)
            if pv:
                mask = (G!=Gr)
                SYN[mask] = 0
                Gv = G[mask]
                isv = torch.zeros(Gv.shape, dtype=torch.float, device=device )
                pw = (Gv<=3) * (3-Gv)
                isv += pw * mus[2] + pw * sigmas[2] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pg = (Gv<=3) * (Gv-2) + (Gv>3) * (4-Gv)
                isv += pg * mus[3] + pg * sigmas[3] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pcsf = (Gv>=3) * (Gv-3)
                isv += pcsf * mus[4] + pcsf * sigmas[4] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                SYN[mask] = isv
            SYN[SYN < 0] = 0


            SYNdef = fast_3D_interp_torch(SYN, xx2, yy2, zz2, 'linear', device)
            SdefOneHot = Ddef = Idef = None
            if S is not None:
                if deform_one_hots:
                    Sonehot = onehotmatrix[lut[S.long()]]
                    SdefOneHot = fast_3D_interp_torch(Sonehot, xx2, yy2, zz2, 'linear', device)
                else:
                    Sdef = fast_3D_interp_torch(S, xx2, yy2, zz2, 'nearest', device)
                    SdefOneHot = onehotmatrix[lut[Sdef.long()]]

            if D is not None:
                Ddef = fast_3D_interp_torch(D, xx2, yy2, zz2, 'linear', device, default_value_linear=torch.max(D))
            if I is not None:
                Idef = fast_3D_interp_torch(I, xx2, yy2, zz2, 'linear', device)

            # Gamma transform
            gamma = torch.tensor(np.exp(gamma_std * np.random.randn(1)[0]), dtype=float, device=device)
            SYNgamma = 300.0 * (SYNdef / 300.0) ** gamma

            # Bias field
            bf_scale = bf_scale_min + np.random.rand(1) * (bf_scale_max - bf_scale_min)
            siz_BF_small = np.round(bf_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_BF_small[1] = np.round(siz[1]/spac).astype(int)
            BFsmall = torch.tensor(bf_std_min + (bf_std_max - bf_std_min) * np.random.rand(1), dtype=torch.float, device=device) * torch.randn(siz_BF_small, dtype=torch.float, device=device)
            BFlog = myzoom_torch(BFsmall, np.array(siz) / siz_BF_small, device)
            BF = torch.exp(BFlog)
            SYNbf = SYNgamma * BF

            # Model Resolution
            stds = (0.85 + 0.3 * np.random.rand()) * np.log(5) /np.pi * thickness / res_training_data
            stds[thickness<=res_training_data] = 0.0 # no blur if thickness is equal to the resolution of the training data
            SYNblur = gaussian_blur_3d(SYNbf, stds, device)
            newsize = (np.array(siz) * res_training_data / resolution).astype(int)

            factors = np.array(newsize) / np.array(siz)
            delta = (1.0 - factors) / (2.0 * factors)
            vx = np.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0])[:newsize[0]]
            vy = np.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1])[:newsize[1]]
            vz = np.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2])[:newsize[2]]
            II, JJ, KK = np.meshgrid(vx, vy, vz, sparse=False, indexing='ij')
            II = torch.tensor(II, dtype=torch.float, device=device)
            JJ = torch.tensor(JJ, dtype=torch.float, device=device)
            KK = torch.tensor(KK, dtype=torch.float, device=device)

            SYNsmall = fast_3D_interp_torch(SYNblur, II, JJ, KK, 'linear', device)
            noise_std = torch.tensor(min_noise_std + (max_noise_std - min_noise_std) * np.random.rand(1), dtype=torch.float, device=device)
            SYNnoisy = SYNsmall + noise_std * torch.randn(SYNsmall.shape, dtype=torch.float, device=device)
            SYNnoisy[SYNnoisy<0] = 0

            # Back to original resolution
            if bspline_zooming:
                SYNresized = interpol.resize(SYNnoisy, shape=siz, anchor='edge', interpolation=3, bound='dct2', prefilter=True)
            else:
                SYNresized = myzoom_torch(SYNnoisy, 1 / factors, device)

            maxi = torch.max(SYNresized)
            SYNfinal = SYNresized / maxi

            if produce_surfaces:
                yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog, Vlw, Flw, Vlp, Flp]
            else:
                yield [SYNfinal, SdefOneHot, Ddef, Idef, BFlog]



def supervised_generator_hemi(datadir,
              label_list_segmentation,
              siz=[96, 128, 128],
              max_rotation=15,
              max_shear=0.2,
              max_scaling=0.2,
              nonlin_scale_min=0.03,
              nonlin_scale_max=0.06,
              nonlin_std_max=4,
              min_noise_std=0, # some noise is ok, but not too much, images already have noise
              max_noise_std=0.02,
              photo_prob=0.2,
              deform_one_hots=True,
              device='cpu'):

    # Paths to the different subdirectories
    seg_dir = os.path.join(datadir, 'label_maps_segmentation')
    im_dir = os.path.join(datadir, 'images')
    if os.path.exists(seg_dir) is False:
        raise Exception('Directory with target segmentations not found; target segmentations will no be generated')
        seg_dir = None
    if os.path.exists(im_dir) is False:
        print('Directory with real images not found; real images will no be generated')
        im_dir = None

    names = glob.glob(os.path.join(seg_dir, '*.nii.gz')) + glob.glob(os.path.join(seg_dir, '*.nii'))
    n = len(names)

    with torch.no_grad():
        # prepare grid
        print('Preparing grid...')
        xx, yy, zz = np.meshgrid(range(siz[0]), range(siz[1]), range(siz[2]), sparse=False, indexing='ij')
        c = torch.tensor((np.array(siz) - 1) / 2, dtype=torch.float, device=device)
        xc = torch.tensor(xx, dtype=torch.float, device=device) - c[0]
        yc = torch.tensor(yy, dtype=torch.float, device=device) - c[1]
        zc = torch.tensor(zz, dtype=torch.float, device=device) - c[2]

        # Matrix for one-hot encoding (includes a lookup-table)
        n_labels = len(label_list_segmentation)
        lut = torch.zeros(10000, dtype=torch.long, device=device)
        for l in range(n_labels):
            lut[label_list_segmentation[l]] = l
        onehotmatrix = torch.eye(n_labels, dtype=torch.float, device=device)

        print('Generator is ready!')


        while True:

            # Select random case
            idx = np.random.randint(n)
            Simg = nib.load(names[idx])
            photo_mode = np.random.rand() < photo_prob

            # The first thing we do is sampling the  deformation, as this will give us a bounding box
            # of the image region we need, so we don't have to read the whole thing from disk
            # (only works for uncompressed niftis!)

            # sample affine deformation
            rotations = (2 * max_rotation * np.random.rand(3) - max_rotation) / 180.0 * np.pi
            shears = (2 * max_shear * np.random.rand(3) - max_shear)
            scalings = 1 + (2 * max_scaling * np.random.rand(3) - max_scaling)
            A = torch.tensor(make_affine_matrix(rotations, shears, scalings), dtype=torch.float, device=device)

            # sample nonlinear deformation
            nonlin_scale = nonlin_scale_min + np.random.rand(1) * (nonlin_scale_max - nonlin_scale_min)
            siz_F_small = np.round(nonlin_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                spac = 2.0 + 10 * np.random.rand()
                siz_F_small[1] = np.round(siz[1]/spac).astype(int)
            nonlin_std = nonlin_std_max * np.random.rand()
            Fsmall = nonlin_std * torch.randn([*siz_F_small, 3], dtype=torch.float, device=device)
            F = myzoom_torch(Fsmall, np.array(siz) / siz_F_small, device)
            if photo_mode:
                F[:, :, :, 1] = 0

            # sample center
            max_shift = (torch.tensor(np.array(Simg.shape[0:3]) - siz, dtype=torch.float, device=device)) / 2
            max_shift[max_shift < 0] = 0
            c2 = torch.tensor((np.array(Simg.shape[0:3]) - 1)/2, dtype=torch.float, device=device) + (2 * (max_shift * torch.rand(3, dtype=float, device=device)) - max_shift)

            # deform (we do nonlinear "first" ie after so we can do heavy coronal deformations in photo mode)
            xx1 = xc + F[:, :, :, 0]
            yy1 = yc + F[:, :, :, 1]
            zz1 = zc + F[:, :, :, 2]
            xx2 = A[0, 0] * xx1 + A[0, 1] * yy1 + A[0, 2] * zz1 + c2[0]
            yy2 = A[1, 0] * xx1 + A[1, 1] * yy1 + A[1, 2] * zz1 + c2[1]
            zz2 = A[2, 0] * xx1 + A[2, 1] * yy1 + A[2, 2] * zz1 + c2[2]
            xx2[xx2 < 0] = 0
            yy2[yy2 < 0] = 0
            zz2[zz2 < 0] = 0
            xx2[xx2 > (Simg.shape[0] - 1)] = Simg.shape[0] - 1
            yy2[yy2 > (Simg.shape[1] - 1)] = Simg.shape[1] - 1
            zz2[zz2 > (Simg.shape[2] - 1)] = Simg.shape[2] - 1

            # Get the margins for reading images
            x1 = torch.floor(torch.min(xx2))
            y1 = torch.floor(torch.min(yy2))
            z1 = torch.floor(torch.min(zz2))
            x2 = 1+torch.ceil(torch.max(xx2))
            y2 = 1 + torch.ceil(torch.max(yy2))
            z2 = 1 + torch.ceil(torch.max(zz2))
            xx2 -= x1
            yy2 -= y1
            zz2 -= z1

            x1 = x1.cpu().numpy().astype(int)
            y1 = y1.cpu().numpy().astype(int)
            z1 = z1.cpu().numpy().astype(int)
            x2 = x2.cpu().numpy().astype(int)
            y2 = y2.cpu().numpy().astype(int)
            z2 = z2.cpu().numpy().astype(int)

            # Read in data
            S  = torch.squeeze(torch.tensor(Simg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
            S[S == 24] = 0
            Iimg = nib.load(os.path.join(im_dir, os.path.basename(names[idx])))
            I = torch.squeeze(torch.tensor(Iimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))

            # Deform
            Idef = fast_3D_interp_torch(I, xx2, yy2, zz2, 'linear', device)
            if deform_one_hots:
                Sonehot = onehotmatrix[lut[S.long()]]
                SdefOneHot = fast_3D_interp_torch(Sonehot, xx2, yy2, zz2, 'linear', device)
            else:
                Sdef = fast_3D_interp_torch(S, xx2, yy2, zz2, 'nearest', device)
                SdefOneHot = onehotmatrix[lut[Sdef.long()]]

            # Noise
            noise_std = torch.tensor(min_noise_std + (max_noise_std - min_noise_std) * np.random.rand(1), dtype=torch.float, device=device)
            Inoisy = torch.clamp(Idef + noise_std * torch.randn(Idef.shape, dtype=torch.float, device=device), min=0)

            yield [Inoisy, SdefOneHot]








#######################
# Auxiliary functions #
#######################

def make_affine_matrix(rot, sh, s):
    Rx = np.array([[1, 0, 0], [0, np.cos(rot[0]), -np.sin(rot[0])], [0, np.sin(rot[0]), np.cos(rot[0])]])
    Ry = np.array([[np.cos(rot[1]), 0, np.sin(rot[1])], [0, 1, 0], [-np.sin(rot[1]), 0, np.cos(rot[1])]])
    Rz = np.array([[np.cos(rot[2]), -np.sin(rot[2]), 0], [np.sin(rot[2]), np.cos(rot[2]), 0], [0, 0, 1]])

    SHx = np.array([[1, 0, 0], [sh[1], 1, 0], [sh[2], 0, 1]])
    SHy = np.array([[1, sh[0], 0], [0, 1, 0], [0, sh[2], 1]])
    SHz = np.array([[1, 0, sh[0]], [0, 1, sh[1]], [0, 0, 1]])

    A = SHx @ SHy @ SHz @ Rx @ Ry @ Rz
    A[0, :] = A[0, :] * s[0]
    A[1, :] = A[1, :] * s[1]
    A[2, :] = A[2, :] * s[2]

    return A

def make_gaussian_kernel(sigma, device):

    sl = int(np.ceil(3 * sigma))
    ts = torch.linspace(-sl, sl, 2*sl+1, dtype=torch.float, device=device)
    gauss = torch.exp((-(ts / sigma)**2 / 2))
    kernel = gauss / gauss.sum()

    return kernel

def gaussian_blur_3d(input, stds, device):
    blurred = input[None, None, :, :, :]
    if stds[0]>0:
        kx = make_gaussian_kernel(stds[0], device=device)
        blurred = conv3d(blurred, kx[None, None, :, None, None], stride=1, padding=(len(kx) // 2, 0, 0))
    if stds[1]>0:
        ky = make_gaussian_kernel(stds[1], device=device)
        blurred = conv3d(blurred, ky[None, None, None, :, None], stride=1, padding=(0, len(ky) // 2, 0))
    if stds[2]>0:
        kz = make_gaussian_kernel(stds[2], device=device)
        blurred = conv3d(blurred, kz[None, None, None, None, :], stride=1, padding=(0, 0, len(kz) // 2))
    return torch.squeeze(blurred)

def fast_3D_interp_torch(X, II, JJ, KK, mode, device, default_value_linear=0.0):
    if mode=='nearest':
        IIr = torch.round(II).long()
        JJr = torch.round(JJ).long()
        KKr = torch.round(KK).long()
        IIr[IIr < 0] = 0
        JJr[JJr < 0] = 0
        KKr[KKr < 0] = 0
        IIr[IIr > (X.shape[0] - 1)] = (X.shape[0] - 1)
        JJr[JJr > (X.shape[1] - 1)] = (X.shape[1] - 1)
        KKr[KKr > (X.shape[2] - 1)] = (X.shape[2] - 1)
        if len(X.shape)==3:
            X = X[..., None]
        Y = torch.zeros([*II.shape, X.shape[3]], dtype=torch.float, device=device)
        for channel in range(X.shape[3]):
            aux = X[:, :, :, channel]
            Y[:,:,:,channel] = aux[IIr, JJr, KKr]
        if Y.shape[3] == 1:
            Y = Y[:, :, :, 0]

    elif mode=='linear':
        ok = (II>0) & (JJ>0) & (KK>0) & (II<=X.shape[0]-1) & (JJ<=X.shape[1]-1) & (KK<=X.shape[2]-1)
        IIv = II[ok]
        JJv = JJ[ok]
        KKv = KK[ok]

        fx = torch.floor(IIv).long()
        cx = fx + 1
        cx[cx > (X.shape[0] - 1)] = (X.shape[0] - 1)
        wcx = IIv - fx
        wfx = 1 - wcx

        fy = torch.floor(JJv).long()
        cy = fy + 1
        cy[cy > (X.shape[1] - 1)] = (X.shape[1] - 1)
        wcy = JJv - fy
        wfy = 1 - wcy

        fz = torch.floor(KKv).long()
        cz = fz + 1
        cz[cz > (X.shape[2] - 1)] = (X.shape[2] - 1)
        wcz = KKv - fz
        wfz = 1 - wcz

        if len(X.shape)==3:
            X = X[..., None]

        Y = torch.zeros([*II.shape, X.shape[3]], dtype=torch.float, device=device)
        for channel in range(X.shape[3]):
            Xc = X[:, :, :, channel]

            c000 = Xc[fx, fy, fz]
            c100 = Xc[cx, fy, fz]
            c010 = Xc[fx, cy, fz]
            c110 = Xc[cx, cy, fz]
            c001 = Xc[fx, fy, cz]
            c101 = Xc[cx, fy, cz]
            c011 = Xc[fx, cy, cz]
            c111 = Xc[cx, cy, cz]

            c00 = c000 * wfx + c100 * wcx
            c01 = c001 * wfx + c101 * wcx
            c10 = c010 * wfx + c110 * wcx
            c11 = c011 * wfx + c111 * wcx

            c0 = c00 * wfy + c10 * wcy
            c1 = c01 * wfy + c11 * wcy

            c = c0 * wfz + c1 * wcz

            Yc = torch.zeros(II.shape, dtype=torch.float, device=device)
            Yc[ok] = c.float()
            Yc[~ok] = default_value_linear
            Y[...,channel] = Yc

        if Y.shape[-1]==1:
            Y = Y[...,0]

    else:
        raise Exception('mode must be linear or nearest')

    return Y



##### Generator for just MNI coordinates
def mni_generator(datadir,
              resolution_sampler,
              siz=[160, 160, 160],
              max_rotation=15,
              max_shear=0.2,
              max_scaling=0.2,
              nonlin_scale_min=0.03,
              nonlin_scale_max=0.06,
              nonlin_std_max=4,
              bf_scale_min=0.02,
              bf_scale_max=0.04,
              bf_std_min=0.1,
              bf_std_max=0.6,
              bag_scale_min=0.02,
              bag_scale_max=0.08,
              gamma_std=0.1,
              min_noise_std=5,
              max_noise_std=15,
              exvixo_prob=0.25,
              photo_prob=0.2,
              bag_prob=0.5,
              pv=True,
              integrate_deformation_fields=False,
              simulate_bags=True,
              bspline_zooming=False,
              device='cpu'):

    if photo_prob>0:
        print('CAREFUL!!!! Photo mode is not properly supported right now')

    # Paths to the different subdirectories
    gen_dir = os.path.join(datadir, 'label_maps_generation')
    mni_dir = os.path.join(datadir, 'MNIreg')
    if simulate_bags is False:
        bag_dir = None
    else:
        bag_dir = os.path.join(datadir, 'DmapsBag')
        if os.path.exists(bag_dir) is False:
            print('Directory with distance maps for bag simulation not found; fake bags will no be generated')
            bag_dir = None

    names = glob.glob(os.path.join(gen_dir, '*.nii.gz')) + glob.glob(os.path.join(gen_dir, '*.nii'))
    n = len(names)

    # Get resolution of training data
    aff = nib.load(names[0]).affine
    res_training_data = np.sqrt(np.sum(aff[:-1, :-1], axis=0))
    n_steps_svf_integration = 8

    with torch.no_grad():
        # prepare grid
        print('Preparing grid...')
        xx, yy, zz = np.meshgrid(range(siz[0]), range(siz[1]), range(siz[2]), sparse=False, indexing='ij')
        xx = torch.tensor(xx, dtype=torch.float, device=device)
        yy = torch.tensor(yy, dtype=torch.float, device=device)
        zz = torch.tensor(zz, dtype=torch.float, device=device)
        c = torch.tensor((np.array(siz) - 1) / 2, dtype=torch.float, device=device)
        xc = xx - c[0]
        yc = yy - c[1]
        zc = zz - c[2]

        print('Generator is ready!')

        while True:

            # Select random case
            idx = np.random.randint(n)
            photo_mode = np.random.rand()<photo_prob
            Gimg = nib.load(names[idx])

            # The first thing we do is sampling the resolution and deformation, as this will give us a bounding box
            # of the image region we need, so we don't have to read the whole thing from disk (only works for uncompressed niftis!

            # Sample resolution
            if photo_mode:
                spac = 2.0 + 10 * np.random.rand()
                resolution = np.array([res_training_data[0], spac , res_training_data[2]])
                thickness = np.array([res_training_data[0], 0.0001, res_training_data[2]])
            else:
                resolution, thickness = resolution_sampler()

            # sample affine deformation
            rotations = (2 * max_rotation * np.random.rand(3) - max_rotation) / 180.0 * np.pi
            shears = (2 * max_shear * np.random.rand(3) - max_shear)
            scalings = 1 + (2 * max_scaling * np.random.rand(3) - max_scaling)
            scaling_factor_distances = np.prod(scalings) ** .33333333333 # we divide distance maps by this, not perfect, but better than nothing
            A = torch.tensor(make_affine_matrix(rotations, shears, scalings), dtype=torch.float, device=device)

            # sample center
            max_shift = (torch.tensor(np.array(Gimg.shape[0:3]) - siz, dtype=torch.float, device=device)) / 2
            max_shift[max_shift < 0] = 0
            c2 = torch.tensor((np.array(Gimg.shape[0:3]) - 1)/2, dtype=torch.float, device=device) + (2 * (max_shift * torch.rand(3, dtype=float, device=device)) - max_shift)

            # sample nonlinear deformation
            nonlin_scale = nonlin_scale_min + np.random.rand(1) * (nonlin_scale_max - nonlin_scale_min)
            siz_F_small = np.round(nonlin_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_F_small[1] = np.round(siz[1]/spac).astype(int)
            nonlin_std = nonlin_std_max * np.random.rand()
            Fsmall = nonlin_std * torch.randn([*siz_F_small, 3], dtype=torch.float, device=device)
            F = myzoom_torch(Fsmall, np.array(siz) / siz_F_small, device)
            if photo_mode:
                F[:, :, :, 1] = 0

            if integrate_deformation_fields:
                steplength = 1.0 / (2.0 ** n_steps_svf_integration)
                Fsvf = F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf += fast_3D_interp_torch(Fsvf, xx + Fsvf[:, :, :, 0], yy + Fsvf[:, :, :, 1], zz + Fsvf[:, :, :, 2], 'linear', device)
                Fsvf_neg = -F * steplength
                for _ in range(n_steps_svf_integration):
                    Fsvf_neg += fast_3D_interp_torch(Fsvf_neg, xx + Fsvf_neg[:, :, :, 0], yy + Fsvf_neg[:, :, :, 1], zz + Fsvf_neg[:, :, :, 2], 'linear', device)
                F = Fsvf
                Fneg = Fsvf_neg

            # deform the images (we do nonlinear "first" ie after so we can do heavy coronal deformations in photo mode)
            xx1 = xc + F[:, :, :, 0]
            yy1 = yc + F[:, :, :, 1]
            zz1 = zc + F[:, :, :, 2]
            xx2 = A[0, 0] * xx1 + A[0, 1] * yy1 + A[0, 2] * zz1 + c2[0]
            yy2 = A[1, 0] * xx1 + A[1, 1] * yy1 + A[1, 2] * zz1 + c2[1]
            zz2 = A[2, 0] * xx1 + A[2, 1] * yy1 + A[2, 2] * zz1 + c2[2]
            xx2[xx2 < 0] = 0
            yy2[yy2 < 0] = 0
            zz2[zz2 < 0] = 0
            xx2[xx2 > (Gimg.shape[0] - 1)] = Gimg.shape[0] - 1
            yy2[yy2 > (Gimg.shape[1] - 1)] = Gimg.shape[1] - 1
            zz2[zz2 > (Gimg.shape[2] - 1)] = Gimg.shape[2] - 1

            # Get the margins for reading images
            x1 = torch.floor(torch.min(xx2))
            y1 = torch.floor(torch.min(yy2))
            z1 = torch.floor(torch.min(zz2))
            x2 = 1+torch.ceil(torch.max(xx2))
            y2 = 1 + torch.ceil(torch.max(yy2))
            z2 = 1 + torch.ceil(torch.max(zz2))
            xx2 -= x1
            yy2 -= y1
            zz2 -= z1

            x1 = x1.cpu().numpy().astype(int)
            y1 = y1.cpu().numpy().astype(int)
            z1 = z1.cpu().numpy().astype(int)
            x2 = x2.cpu().numpy().astype(int)
            y2 = y2.cpu().numpy().astype(int)
            z2 = z2.cpu().numpy().astype(int)


            # Read in data
            G = torch.squeeze(torch.tensor(Gimg.get_fdata()[x1:x2, y1:y2, z1:z2].astype(float), dtype=torch.float, device=device))
            MNIimg = nib.load(os.path.join(mni_dir, os.path.basename(names[idx])))
            MNI = torch.squeeze(torch.tensor(MNIimg.get_fdata()[x1:x2, y1:y2, z1:z2, :].astype(float), dtype=torch.float, device=device))
            MNI /= 100
            if bag_dir is not None:
                Bimg = nib.load(os.path.join(bag_dir, os.path.basename(names[idx])))
                B = torch.squeeze(torch.tensor(Bimg.get_fdata()[x1:x2, y1:y2, z1:z2], dtype=torch.float, device=device))
                B /= scaling_factor_distances

            # Decide if we're simulating ex vivo (and possibly a bag) or photos
            BRAINMASK = torch.ones_like(G)
            BRAINMASK[G == 0] = 0
            BRAINMASK[G > 255] = 0

            if photo_mode or (np.random.rand() < exvixo_prob):

                G[G>255] = 0 # kill extracerebral
                if photo_mode:
                    tokill = (G==7); G[tokill] = 0; BRAINMASK[tokill] = 0;
                    tokill = (G==8); G[tokill] = 0; BRAINMASK[tokill] = 0;
                    tokill = (G==16); G[tokill] = 0; BRAINMASK[tokill] = 0;
                    # without distance maps, killing 4 is the best we can do
                    G[G == 4] = 0

                elif ((B is not None) and (np.random.rand(1) < bag_prob)):
                    bag_scale = bag_scale_min + np.random.rand(1) * (bag_scale_max - bag_scale_min)
                    siz_TH_small = np.round(bag_scale * np.array(G.shape)).astype(int).tolist()
                    bag_tness = torch.tensor(np.sort(1.0 + 20 * np.random.rand(2)), dtype=torch.float, device=device)
                    THsmall = bag_tness[0] + (bag_tness[1] - bag_tness[0]) * torch.rand(siz_TH_small, dtype=torch.float, device=device)
                    TH = myzoom_torch(THsmall, np.array(G.shape) / siz_TH_small, device)
                    G[(B>0) & (B<TH)] = 4

            # Sample Gaussian image
            mus = 25 + 200 * torch.rand(10000, dtype=torch.float, device=device)
            sigmas = 5 + 20 * torch.rand(10000, dtype=torch.float, device=device)
            if photo_mode or np.random.rand(1)<0.5: # set the background to zero every once in a while (or always in photo mode)
                mus[0] = 0
            Gr = torch.round(G).long()
            SYN = mus[Gr] + sigmas[Gr] * torch.randn(Gr.shape, dtype=torch.float, device=device)
            if pv:
                mask = (G!=Gr)
                SYN[mask] = 0
                Gv = G[mask]
                isv = torch.zeros(Gv.shape, dtype=torch.float, device=device )
                pw = (Gv<=3) * (3-Gv)
                isv += pw * mus[2] + pw * sigmas[2] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pg = (Gv<=3) * (Gv-2) + (Gv>3) * (4-Gv)
                isv += pg * mus[3] + pg * sigmas[3] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                pcsf = (Gv>=3) * (Gv-3)
                isv += pcsf * mus[4] + pcsf * sigmas[4] * torch.randn(Gv.shape, dtype=torch.float, device=device)
                SYN[mask] = isv
            SYN[SYN < 0] = 0


            SYNdef = fast_3D_interp_torch(SYN, xx2, yy2, zz2, 'linear', device)
            MNIdef = fast_3D_interp_torch(MNI, xx2, yy2, zz2, 'linear', device)
            BRAINMASKdef = fast_3D_interp_torch(BRAINMASK, xx2, yy2, zz2, 'nearest', device)

            # Gamma transform
            gamma = torch.tensor(np.exp(gamma_std * np.random.randn(1)[0]), dtype=float, device=device)
            SYNgamma = 300.0 * (SYNdef / 300.0) ** gamma

            # Bias field
            bf_scale = bf_scale_min + np.random.rand(1) * (bf_scale_max - bf_scale_min)
            siz_BF_small = np.round(bf_scale * np.array(siz)).astype(int).tolist()
            if photo_mode:
                siz_BF_small[1] = np.round(siz[1]/spac).astype(int)
            BFsmall = torch.tensor(bf_std_min + (bf_std_max - bf_std_min) * np.random.rand(1), dtype=torch.float, device=device) * torch.randn(siz_BF_small, dtype=torch.float, device=device)
            BFlog = myzoom_torch(BFsmall, np.array(siz) / siz_BF_small, device)
            BF = torch.exp(BFlog)
            SYNbf = SYNgamma * BF

            # Model Resolution
            stds = (0.85 + 0.3 * np.random.rand()) * np.log(5) /np.pi * thickness / res_training_data
            stds[thickness<=res_training_data] = 0.0 # no blur if thickness is equal to the resolution of the training data
            SYNblur = gaussian_blur_3d(SYNbf, stds, device)
            newsize = (np.array(siz) * res_training_data / resolution).astype(int)

            factors = np.array(newsize) / np.array(siz)
            delta = (1.0 - factors) / (2.0 * factors)
            vx = np.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0])[:newsize[0]]
            vy = np.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1])[:newsize[1]]
            vz = np.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2])[:newsize[2]]
            II, JJ, KK = np.meshgrid(vx, vy, vz, sparse=False, indexing='ij')
            II = torch.tensor(II, dtype=torch.float, device=device)
            JJ = torch.tensor(JJ, dtype=torch.float, device=device)
            KK = torch.tensor(KK, dtype=torch.float, device=device)

            SYNsmall = fast_3D_interp_torch(SYNblur, II, JJ, KK, 'linear', device)
            noise_std = torch.tensor(min_noise_std + (max_noise_std - min_noise_std) * np.random.rand(1), dtype=torch.float, device=device)
            SYNnoisy = SYNsmall + noise_std * torch.randn(SYNsmall.shape, dtype=torch.float, device=device)
            SYNnoisy[SYNnoisy<0] = 0

            # Back to original resolution
            if bspline_zooming:
                SYNresized = interpol.resize(SYNnoisy, shape=siz, anchor='edge', interpolation=3, bound='dct2', prefilter=True)
            else:
                SYNresized = myzoom_torch(SYNnoisy, 1 / factors, device)
            maxi = torch.max(SYNresized)
            SYNfinal = SYNresized / maxi

            # mask real image and MNI coordinates if needed
            MNIdef *= BRAINMASKdef[..., None]

            # Flip 50% of times
            if np.random.rand()<0.5:
                SYNfinal = torch.flip(SYNfinal, [0])
                MNIdef = torch.flip(MNIdef, [0])
                MNIdef[:, :, :, 0] = -MNIdef[:, :, :, 0] # pretty easy thanks to symmetric template
                # BFlog = torch.flip(BFlog, [0]) we don't need this

            yield [SYNfinal, MNIdef]



