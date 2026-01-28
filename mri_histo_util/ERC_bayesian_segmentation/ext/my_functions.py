import os
import numpy as np
import nibabel as nib
import cv2
import torch
import scipy.ndimage
import scipy.sparse as sp
from torch.nn import functional
from torch.utils.data import Dataset, DataLoader
from skimage.measure import label

###############################

def get_ras_axes(aff, n_dims=3):
    """This function finds the RAS axes corresponding to each dimension of a volume, based on its affine matrix.
    :param aff: affine matrix Can be a 2d numpy array of size n_dims*n_dims, n_dims+1*n_dims+1, or n_dims*n_dims+1.
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: two numpy 1d arrays of lengtn n_dims, one with the axes corresponding to RAS orientations,
    and one with their corresponding direction.
    """
    aff_inverted = np.linalg.inv(aff)
    img_ras_axes = np.argmax(np.absolute(aff_inverted[0:n_dims, 0:n_dims]), axis=0)
    return img_ras_axes

###############################

def align_volume_to_ref(volume, aff, aff_ref=None, return_aff=False, n_dims=3):
    """This function aligns a volume to a reference orientation (axis and direction) specified by an affine matrix.
    :param volume: a numpy array
    :param aff: affine matrix of the floating volume
    :param aff_ref: (optional) affine matrix of the target orientation. Default is identity matrix.
    :param return_aff: (optional) whether to return the affine matrix of the aligned volume
    :param n_dims: number of dimensions (excluding channels) of the volume corresponding to the provided affine matrix.
    :return: aligned volume, with corresponding affine matrix if return_aff is True.
    """

    # work on copy
    aff_flo = aff.copy()

    # default value for aff_ref
    if aff_ref is None:
        aff_ref = np.eye(4)

    # extract ras axes
    ras_axes_ref = get_ras_axes(aff_ref, n_dims=n_dims)
    ras_axes_flo = get_ras_axes(aff_flo, n_dims=n_dims)

    # align axes
    aff_flo[:, ras_axes_ref] = aff_flo[:, ras_axes_flo]
    for i in range(n_dims):
        if ras_axes_flo[i] != ras_axes_ref[i]:
            volume = torch.swapaxes(volume, ras_axes_flo[i], ras_axes_ref[i])
            swapped_axis_idx = np.where(ras_axes_flo == ras_axes_ref[i])
            ras_axes_flo[swapped_axis_idx], ras_axes_flo[i] = ras_axes_flo[i], ras_axes_flo[swapped_axis_idx]

    # align directions
    dot_products = np.sum(aff_flo[:3, :3] * aff_ref[:3, :3], axis=0)
    for i in range(n_dims):
        if dot_products[i] < 0:
            volume = torch.flip(volume, [i])
            aff_flo[:, i] = - aff_flo[:, i]
            aff_flo[:3, 3] = aff_flo[:3, 3] - aff_flo[:3, i] * (volume.shape[i] - 1)

    if return_aff:
        return volume, aff_flo
    else:
        return volume

###############################3

# Crop label volume
def cropLabelVol(V,
                 margin=10,
                 threshold=0):

    # Make sure it's 3D
    margin = np.array(margin)
    if len(margin.shape) < 2:
        margin = [margin, margin, margin]

    if len(V.shape) < 2:
        V = V[..., np.newaxis]
    if len(V.shape) < 3:
        V = V[..., np.newaxis]

    # Now
    idx = np.where(V > threshold)
    i1 = np.max([0, np.min(idx[0]) - margin[0]]).astype('int')
    j1 = np.max([0, np.min(idx[1]) - margin[1]]).astype('int')
    k1 = np.max([0, np.min(idx[2]) - margin[2]]).astype('int')
    i2 = np.min([V.shape[0], np.max(idx[0]) + margin[0] + 1]).astype('int')
    j2 = np.min([V.shape[1], np.max(idx[1]) + margin[1] + 1]).astype('int')
    k2 = np.min([V.shape[2], np.max(idx[2]) + margin[2] + 1]).astype('int')

    cropping = [i1, j1, k1, i2, j2, k2]
    cropped = V[i1:i2, j1:j2, k1:k2]

    return cropped, cropping

##################
def cropLabelVolTorch(V,
                 margin=10,
                 threshold=0):

    margin = torch.tensor(margin, device=V.device, dtype=torch.long)
    if len(V.shape) < 2:
        V = V[..., None]
    if len(V.shape) < 3:
        V = V[..., None]
    # Now crop
    idx = torch.where(V > threshold)
    i1 = torch.min(idx[0]) - margin
    j1 = torch.min(idx[1]) - margin
    k1 = torch.min(idx[2]) - margin
    i2 = torch.max(idx[0]) + margin
    j2 = torch.max(idx[1]) + margin
    k2 = torch.max(idx[2]) + margin
    # out of bounds check
    i1 = i1 if i1>=0 else torch.tensor(0, device=V.device, dtype=torch.long)
    j1 = j1 if j1 >= 0 else torch.tensor(0, device=V.device, dtype=torch.long)
    k1 = k1 if k1 >= 0 else torch.tensor(0, device=V.device, dtype=torch.long)
    i2 = i2 if i2 < V.shape[0] else torch.tensor(V.shape[0] - 1, device=V.device, dtype=torch.long)
    j2 = j2 if j2 < V.shape[1] else torch.tensor(V.shape[0] - 1, device=V.device, dtype=torch.long)
    k2 = k2 if k2 < V.shape[2] else torch.tensor(V.shape[0] - 1, device=V.device, dtype=torch.long)

    cropping = [i1, j1, k1, i2, j2, k2]
    cropped = V[i1:i2, j1:j2, k1:k2]

    return cropped, cropping

###############################3

def applyCropping(V, cropping):
    i1 = cropping[0]
    j1 = cropping[1]
    k1 = cropping[2]
    i2 = cropping[3]
    j2 = cropping[4]
    k2 = cropping[5]

    if len(V.shape)>2:
        Vcropped = V[i1:i2, j1: j2, k1: k2, ...]
    else:
        Vcropped = V[i1:i2, j1: j2]

    return Vcropped

###############################3

def viewVolume(x, aff=None):

    if aff is None:
        aff = np.eye(4)
    else:
        if type(aff) == torch.Tensor:
            aff = aff.detach().cpu().numpy()

    if type(x) is not list:
        x = [x]

    cmd = 'source /usr/local/freesurfer/nmr-dev-env-bash && freeview '

    for n in np.arange(len(x)):
        vol = x[n]
        if type(vol) == torch.Tensor:
            vol = vol.detach().cpu().numpy()
        vol = np.squeeze(np.array(vol))
        name = '/tmp/' + str(n) + '.nii.gz'
        MRIwrite(vol, aff, name)
        cmd = cmd + ' ' + name

    os.system(cmd + ' &')

###############################3

def getLargestCC(segmentation):
    labels = label(segmentation)
    largestCC = labels == np.argmax(np.bincount(labels.flat, weights=segmentation.flat))
    return largestCC

###############################3

def MRIwrite(volume, aff, filename, dtype=None):

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)

    nib.save(nifty, filename)

###############################3

def MRIread(filename, dtype=None, im_only=False, as_closest_canonical=False):

    assert filename.endswith(('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    if as_closest_canonical:
        x = nib.as_closest_canonical(x)
    volume = x.get_fdata()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

###############################

def get_largest_connected_component(binary_numpy):
    labeled_array, num_features = scipy.ndimage.label(binary_numpy)
    if num_features == 0:
      return np.zeros_like(binary_numpy)
    component_sizes = np.bincount(labeled_array.flatten())
    largest_component_label = np.argmax(component_sizes[1:]) + 1
    largest_component_mask = (labeled_array == largest_component_label)
    return largest_component_mask

###############################

def getM(ref, mov):
    device = ref.device
    dtype = torch.float
    zmat = torch.zeros(ref.shape[::-1], device=device, dtype=dtype)
    zcol = torch.zeros([ref.shape[1], 1], device=device, dtype=dtype)
    ocol = torch.ones([ref.shape[1], 1], device=device, dtype=dtype)
    zero = torch.zeros(zmat.shape, device=device, dtype=dtype)
    A = torch.concatenate([
        torch.concatenate([torch.t(ref), zero, zero, ocol, zcol, zcol], axis=1),
        torch.concatenate([zero, torch.t(ref), zero, zcol, ocol, zcol], axis=1),
        torch.concatenate([zero, zero, torch.t(ref), zcol, zcol, ocol], axis=1)], axis=0)
    b = torch.concatenate([torch.t(mov[0, :]), torch.t(mov[1, :]), torch.t(mov[2, :])], axis=0)
    x = (torch.linalg.inv(torch.t(A) @ A))  @ (torch.t(A) @ b)
    M = torch.tensor([
        [x[0], x[1], x[2], x[9]],
        [x[3], x[4], x[5], x[10]],
        [x[6], x[7], x[8], x[11]],
        [0, 0, 0, 1]], device=device, dtype=dtype)
    return M


###############################
def fast_3D_interp_torch(X, II, JJ, KK, mode, pad_value=0):
    if mode=='nearest':
        IIr = torch.round(II).long()
        JJr = torch.round(JJ).long()
        KKr = torch.round(KK).long()

        ok = torch.full(II.shape, True, device=X.device, dtype=torch.bool)

        mask = (IIr < 0)
        ok[mask] = False
        IIr[mask] = 0

        mask = (JJr < 0)
        ok[mask] = False
        JJr[mask] = 0

        mask = (KKr < 0)
        ok[mask] = False
        KKr[mask] = 0

        mask = (IIr > (X.shape[0] - 1))
        ok[mask] = False
        IIr[mask] = (X.shape[0] - 1)

        mask = (JJr > (X.shape[1] - 1))
        ok[mask] = False
        JJr[mask] = (X.shape[1] - 1)

        mask = (KKr > (X.shape[2] - 1))
        ok[mask] = False
        KKr[mask] = (X.shape[2] - 1)

        Y = X[IIr, JJr, KKr]

        Y[ok==False] = pad_value

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

        c000 = X[fx, fy, fz]
        c100 = X[cx, fy, fz]
        c010 = X[fx, cy, fz]
        c110 = X[cx, cy, fz]
        c001 = X[fx, fy, cz]
        c101 = X[cx, fy, cz]
        c011 = X[fx, cy, cz]
        c111 = X[cx, cy, cz]

        c00 = c000 * wfx + c100 * wcx
        c01 = c001 * wfx + c101 * wcx
        c10 = c010 * wfx + c110 * wcx
        c11 = c011 * wfx + c111 * wcx

        c0 = c00 * wfy + c10 * wcy
        c1 = c01 * wfy + c11 * wcy

        c = c0 * wfz + c1 * wcz

        Y = torch.full(II.shape, pad_value, device=X.device, dtype=X.dtype)
        Y[ok] = c.float()

    else:
        raise Exception('mode must be linear or nearest')

    return Y


################################

def downsampleMRI2d(X, aff, shape, factors, mode='image'):

    assert False, 'Function not debugged/tested yet...'

    assert mode=='image' or mode=='labels', 'Mode must be image or labels'
    assert (shape is None) or (factors is None), 'Either shape or factors must be None'
    assert (shape is not None) or (factors is not None), 'Either shape or factors must be not None'

    if shape is not None:
        factors = np.array(shape) / X.shape[0:2]
    else:
        factors = np.array(factors)
        shape = np.round(X.shape[0:2] * factors).astype('int')

    if mode == 'image':
        if np.mean(factors) < 1: # shrink
            Y = cv2.resize(X, shape, interpolation=cv2.INTER_AREA)
        else:  # expan
            Y = cv2.resize(X, shape, interpolation=cv2.INTER_LINEAR)
    else:
        Y = cv2.resize(X, shape, interpolation=cv2.INTER_NEAREST)

    aff2 = aff
    aff2[:, 0] = aff2[:, 0] * factors[0]
    aff2[:, 1] = aff2[:, 1] * factors[1]
    aff2[0:3, 3] = aff2[0:3, 3] + aff[0:3, 0:3] * (0.5*np.array([[factors[0]], [factors[1]], [1]])-0.5)

    return Y, aff2

###############################3

def vox2ras(vox, vox2ras):

    vox2 = np.concatenate([vox, np.ones(shape=[1, vox.shape[1]])], axis=0)

    ras = np.matmul(vox2ras, vox2)[:-1, :]

    return ras

###############################

def ras2vox(ras, vox2ras):

    ras2 = np.concatenate([ras, np.ones(shape=[1, ras.shape[1]])], axis=0)

    vox = np.matmul(np.linalg.inv(vox2ras), ras2)[:-1, :]

    return vox


###############################3

def prepBiasFieldBase3d(siz, max_order):
    x = np.linspace(-1, 1, siz[0])
    y = np.linspace(-1, 1, siz[1])
    z = np.linspace(-1, 1, siz[2])
    xx, yy, zz = np.meshgrid(x, y, z, indexing='ij')
    PSI = []
    for o in range(max_order + 1):
        for ox in range(o + 1):
            for oy in range(o + 1):
                for oz in range(o + 1):
                    if (ox + oy + oz) == o:
                        psi = np.ones(siz)
                        for i in range(1, ox + 1):
                            psi = psi * xx
                        for j in range(1, oy + 1):
                            psi = psi * yy
                        for k in range(1, oz + 1):
                            psi = psi * zz
                        PSI.append(psi)

    PSI = np.stack(PSI, axis=-1)

    return PSI

###############################3

def grad3d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3, 1]))
    Gz = scipy.ndimage.convolve(X, np.reshape(h, [1, 1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy + Gz * Gz)

    if provide_gradients:
        return Gmodule, Gx, Gy, Gz
    else:
        return Gmodule

###############################3

def grad2d(X, provide_gradients=False):
    h = np.array([-1, 0, 1])
    Gx = scipy.ndimage.convolve(X, np.reshape(h, [3, 1]))
    Gy = scipy.ndimage.convolve(X, np.reshape(h, [1, 3]))
    Gmodule = np.sqrt(Gx * Gx + Gy * Gy)

    if provide_gradients:
        return Gmodule, Gx, Gy
    else:
        return Gmodule


########################
def torch_resize(I, aff, resolution, device, power_factor_at_half_width=5, dtype=torch.float32, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_resize(I, aff, resolution, device, power_factor_at_half_width, dtype, slow)

    slow = slow or device.type == 'cpu'
    voxsize = np.sqrt(np.sum(aff[:-1, :-1] ** 2, axis=0))
    newsize = np.round(I.shape[0:3] * (voxsize / resolution)).astype(int)
    factors = np.array(I.shape[0:3]) / np.array(newsize)
    k = np.log(power_factor_at_half_width) / np.pi
    sigmas = k * factors
    sigmas[sigmas<=k] = 0  # TODO: we could maybe remove this line, to make sure we always smooth a bit?

    if len(I.shape) not in (3, 4):
        raise Exception('torch_resize works with 3D or 3D+label volumes')
    no_channels = len(I.shape) == 3
    if no_channels:
        I = I[:, :, :, None]
    if torch.is_tensor(I):
        I = I.permute([3, 0, 1, 2])
    else:
        I = I.transpose([3, 0, 1, 2])

    It_lowres = None
    for c in range(len(I)):
        It = torch.as_tensor(I[c], device=device, dtype=dtype)[None, None]
        # Smoothen if needed
        for d in range(3):
            if sigmas[d]>0:
                sl = np.ceil(sigmas[d] * 2.5).astype(int)
                v = np.arange(-sl, sl + 1)
                gauss = np.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / np.sum(gauss)
                kernel = torch.tensor(kernel,  device=device, dtype=dtype)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None,  None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])

            It = It.permute([0, 1, 4, 2, 3])
        It = torch.squeeze(It)
        It, aff2 = myzoom_torch(It, aff, newsize, device)
        It = It.detach()
        if torch.is_tensor(I):
            It = It.to(I.device)
        else:
            It = It.cpu().numpy()
        if len(I) == 1:
            It_lowres = It[None]
        else:
            if It_lowres is None:
                if torch.is_tensor(It):
                    It_lowres = It.new_empty([len(I), *It.shape])
                else:
                    It_lowres = np.empty_like(It, shape=[len(I), *It.shape])
            It_lowres[c] = It

        torch.cuda.empty_cache()

    if not no_channels:
        if torch.is_tensor(I):
            It_lowres = It_lowres.permute([1, 2, 3, 0])
        else:
            It_lowres = It_lowres.transpose([1, 2, 3, 0])
    else:
        It_lowres = It_lowres[0]

    return It_lowres, aff2

def torch_gaussian_filter3d(I, sigmas, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_gaussial_filter3d(I, sigmas, slow)
    device = I.device
    dtype = I.dtype
    slow = slow or device.type == 'cpu'
    if len(sigmas)==1:
        sigmas = sigmas * torch.ones(3, device=device, dtype=dtype)

    if len(I.shape) not in (3, 4):
        raise Exception('torch_resize works with 3D or 3D+label volumes')
    no_channels = len(I.shape) == 3
    if no_channels:
        I = I[:, :, :, None]
    I = I.permute([3, 0, 1, 2])

    It_lowres = None
    for c in range(len(I)):
        It = torch.as_tensor(I[c], device=device, dtype=dtype)[None, None]
        # Smoothen if needed
        for d in range(3):
            if sigmas[d]>0:
                sl = torch.ceil(sigmas[d] * 2.5).to(device).to(int)
                v = torch.arange(-sl, sl + 1).to(device).to(dtype)
                gauss = torch.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / torch.sum(gauss)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None,  None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])

            It = It.permute([0, 1, 4, 2, 3])
        It = torch.squeeze(It)
        It = It.detach()
        It = It.to(I.device)
        if len(I) == 1:
            It_smooth = It[None]
        else:
            if It_smooth is None:
                It_smooth = It.new_empty([len(I), *It.shape])
            It_smooth[c] = It

        torch.cuda.empty_cache()

    if not no_channels:
        It_smooth = It_smooth.permute([1, 2, 3, 0])
    else:
        It_smooth = It_smooth[0]

    return It_smooth

@torch.jit.script
def conv_slow_fallback(x, kernel):
    """1D Conv along the last dimension with padding"""
    y = torch.zeros_like(x)
    x = torch.nn.functional.pad(x, [(len(kernel) - 1) // 2]*2)
    x = x.unfold(-1, size=len(kernel), step=1)
    x = x.movedim(-1, 0)
    for i in range(len(kernel)):
        y = y.addcmul_(x[i], kernel[i])
    return y

########################

def myzoom_torch(X, aff, newsize, device):

    if len(X.shape)==3:
        X = X[..., None]

    factors = np.array(newsize) / np.array(X.shape[:-1])
    delta = (1.0 - factors) / (2.0 * factors)

    vx = torch.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0], dtype=torch.float, device=device)[:newsize[0]]
    vy = torch.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1], dtype=torch.float, device=device)[:newsize[1]]
    vz = torch.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2], dtype=torch.float, device=device)[:newsize[2]]

    vx[vx < 0] = 0
    vy[vy < 0] = 0
    vz[vz < 0] = 0
    vx[vx > (X.shape[0]-1)] = (X.shape[0]-1)
    vy[vy > (X.shape[1] - 1)] = (X.shape[1] - 1)
    vz[vz > (X.shape[2] - 1)] = (X.shape[2] - 1)

    fx = torch.floor(vx).int()
    cx = fx + 1
    cx[cx > (X.shape[0]-1)] = (X.shape[0]-1)
    wcx = vx - fx
    wfx = 1 - wcx

    fy = torch.floor(vy).int()
    cy = fy + 1
    cy[cy > (X.shape[1]-1)] = (X.shape[1]-1)
    wcy = vy - fy
    wfy = 1 - wcy

    fz = torch.floor(vz).int()
    cz = fz + 1
    cz[cz > (X.shape[2]-1)] = (X.shape[2]-1)
    wcz = vz - fz
    wfz = 1 - wcz

    Y = torch.zeros([newsize[0], newsize[1], newsize[2], X.shape[3]], dtype=torch.float, device=device)

    dtype = X.dtype
    for channel in range(X.shape[3]):
        Xc = X[:,:,:,channel]

        tmp1 = torch.zeros([newsize[0], Xc.shape[1], Xc.shape[2]], dtype=dtype, device=device)
        for i in range(newsize[0]):
            tmp1[i, :, :] = wfx[i] * Xc[fx[i], :, :] +  wcx[i] * Xc[cx[i], :, :]
        tmp2 = torch.zeros([newsize[0], newsize[1], Xc.shape[2]], dtype=dtype, device=device)
        for j in range(newsize[1]):
            tmp2[:, j, :] = wfy[j] * tmp1[:, fy[j], :] +  wcy[j] * tmp1[:, cy[j], :]
        for k in range(newsize[2]):
            Y[:, :, k, channel] = wfz[k] * tmp2[:, :, fz[k]] +  wcz[k] * tmp2[:, :, cz[k]]

    if Y.shape[3] == 1:
        Y = Y[:,:,:, 0]

    if aff is not None:
        aff_new = aff.copy()
        for c in range(3):
            aff_new[:-1, c] = aff_new[:-1, c] / factors[c]
        aff_new[:-1, -1] = aff_new[:-1, -1] - aff[:-1, :-1] @ (0.5 - 0.5 / factors)
        return Y, aff_new
    else:
        return Y


def read_LUT(LUT_file):

    names = [None] * 10000
    colors = [None] * 10000
    with open(LUT_file, 'r') as file:
        lines = file.readlines()
    for line in lines:
        split = line.lstrip().split()
        if split and not split[0].startswith('#'):
            index, name = split[:2]
            color = np.asarray(list(map(int, split[2:5])), dtype=np.float64)
            idx = int(index)
            names[idx] = name
            colors[idx] = color

    return names, colors

####

def flip_left_right_labels(X):

    lut = np.zeros(3000)

    lut[2] = 41; lut[4] = 43; lut[5] = 44; lut[7] = 46; lut[8] = 47; lut[10] = 49; lut[11] = 50; lut[12] = 51; lut[13] = 52; lut[14] = 14;
    lut[15] = 15; lut[16] = 16; lut[17] = 53; lut[18] = 54; lut[24] = 24; lut[26] = 58; lut[28] = 60; lut[41] = 2; lut[43] = 4; lut[44] = 5;
    lut[46] = 7; lut[47] = 8; lut[49] = 10; lut[50] = 11; lut[51] = 12; lut[52] = 13; lut[53] = 17; lut[54] = 18; lut[58] = 26; lut[60] = 28;
    lut[1000:2000] = range(2000, 3000)
    lut[2000:3000] = range(1000, 2000)

    Y = lut[X.astype(int)]

    return Y

def cumprod(sequence, reverse=False, exclusive=False):
    """Perform the cumulative product of a sequence of elements.

    Parameters
    ----------
    sequence : any object that implements `__iter__`
        Sequence of elements for which the `__mul__` operator is defined.
    reverse : bool, default=False
        Compute cumulative product from right-to-left:
        `cumprod([a, b, c], reverse=True) -> [a*b*c, b*c, c]`
    exclusive : bool, default=False
        Exclude self from the cumulative product:
        `cumprod([a, b, c], exclusive=True) -> [1, a, a*b]`

    Returns
    -------
    product : list
        Product of the elements in the sequence.

    """
    if reverse:
        sequence = reversed(sequence)
    accumulate = None
    seq = [1] if exclusive else []
    for elem in sequence:
        if accumulate is None:
            accumulate = elem
        else:
            accumulate = accumulate * elem
        seq.append(accumulate)
    if exclusive:
        seq = seq[:-1]
    if reverse:
        seq = list(reversed(seq))
    return seq


def sub2ind(subs, shape, out=None):
    """Convert sub indices (i, j, k) into linear indices.

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    subs : (D, ...) tensor
        List of sub-indices. The first dimension is the number of dimension.
        Each element should have the same number of elements and shape.
    shape : (D,) vector_like
        Size of each dimension. Its length should be the same as the
        first dimension of ``subs``.
    out : tensor, optional
        Output placeholder

    Returns
    -------
    ind : (...) tensor
        Linear indices
    """
    *subs, ind = subs
    if out is None:
        ind = ind.clone()
    else:
        out.reshape(ind.shape).copy_(ind)
        ind = out
    backend = dict(dtype=ind.dtype, device=ind.device)
    stride = cumprod(shape[1:], reverse=True)
    for i, s in zip(subs, stride):
        ind += torch.as_tensor(i, **backend) * torch.as_tensor(s, **backend)
    return ind


def ind2sub(ind, shape, out=None):
    """Convert linear indices into sub indices (i, j, k).

    The rightmost dimension is the most rapidly changing one
    -> if shape == [D, H, W], the strides are therefore [H*W, W, 1]

    Parameters
    ----------
    ind : tensor
        Linear indices
    shape : (D,) vector_like
        Size of each dimension.
    out : tensor, optional
        Output placeholder

    Returns
    -------
    subs : (D, ...) tensor
        Sub-indices.
    """
    backend = dict(dtype=ind.dtype, device=ind.device)
    stride = cumprod(shape, reverse=True, exclusive=True)
    stride = torch.as_tensor(stride, **backend)
    if out is None:
        sub = ind.new_empty([len(shape), *ind.shape])
    else:
        sub = out.reshape([len(shape), *ind.shape])
    sub[:, ...] = ind
    for d in range(len(shape)):
        if d > 0:
            torch.remainder(sub[d], torch.as_tensor(stride[d-1], **backend), out=sub[d])
        sub[d] = torch.div(sub[d], stride[d], out=sub[d], rounding_mode='trunc')
    return sub


def affine_sub(affine, shape, indices):
    """Update an affine matrix according to a sub-indexing of the lattice.

    Notes
    -----
    .. Only sub-indexing that *keep an homogeneous voxel size* are allowed.
       Therefore, indices must be `None` or of type `int`, `slice`, `ellipsis`.

    Parameters
    ----------
    affine : (..., ndim_out[+1], ndim_in+1) tensor
        Input affine matrix.
    shape : (ndim_in,) sequence[int]
        Input shape.
    indices : tuple[slice or ellipsis]
        Subscripting indices.

    Returns
    -------
    affine : (..., ndim_out[+1], ndim_new+1) tensor
        Updated affine matrix.
    shape : (ndim_new,) tuple[int]
        Updated shape.

    """
    def is_int(elem):
        if torch.is_tensor(elem):
            return elem.dtype in (torch.int32, torch.int64)
        elif isinstance(elem, int):
            return True
        else:
            return False

    def to_int(elem):
        if torch.is_tensor(elem):
            return elem.item()
        else:
            assert isinstance(elem, int)
            return elem

    # check types
    nb_dim = affine.shape[-1] - 1
    backend = dict(dtype=affine.dtype, device=affine.device)
    if torch.is_tensor(shape):
        shape = shape.tolist()
    if len(shape) != nb_dim:
        raise ValueError('Expected shape of length {}. Got {}'
                         .format(nb_dim, len(shape)))
    if not isinstance(indices, tuple):
        raise TypeError('Indices should be a tuple.')
    indices = list(indices)

    # compute the number of input dimension that correspond to each index
    #   > slice index one dimension but eliipses index multiple dimension
    #     and their number must be computed.
    nb_dims_in = []
    ind_ellipsis = None
    for n_ind, ind in enumerate(indices):
        if isinstance(ind, slice):
            nb_dims_in.append(1)
        elif ind is Ellipsis:
            if ind_ellipsis is not None:
                raise ValueError('Cannot have more than one ellipsis.')
            ind_ellipsis = n_ind
            nb_dims_in.append(-1)
        elif is_int(ind):
            nb_dims_in.append(1)
        elif ind is None:
            nb_dims_in.append(0)
        else:
            raise TypeError('Indices should be None, integers, slices or '
                            'ellipses. Got {}.'.format(type(ind)))
    nb_known_dims = sum(nb_dims for nb_dims in nb_dims_in if nb_dims > 0)
    if ind_ellipsis is not None:
        nb_dims_in[ind_ellipsis] = max(0, nb_dim - nb_known_dims)

    # transform each index into a slice
    # note that we don't need to know "stop" to update the affine matrix
    nb_ind = 0
    indices0 = indices
    indices = []
    for d, ind in enumerate(indices0):
        if isinstance(ind, slice):
            start = ind.start
            step = ind.step
            step = 1 if step is None else step
            start = 0 if (start is None and step > 0) else \
                    shape[nb_ind] - 1 if (start is None and step < 0) else \
                    shape[nb_ind] + start if start < 0 else \
                    start
            indices.append(slice(start, None, step))
            nb_ind += 1
        elif ind is Ellipsis:
            for dd in range(nb_ind, nb_ind + nb_dims_in[d]):
                start = 0
                step = 1
                indices.append(slice(start, None, step))
                nb_ind += 1
        elif is_int(ind):
            indices.append(to_int(ind))
        elif ind is None:
            assert (ind is None), "Strange index of type {}".format(type(ind))
            indices.append(None)

    # Extract shift and scale in each dimension
    shifts = []
    scales = []
    slicer = []
    shape_out = []
    for d, ind in enumerate(indices):
        # translation + scale
        if isinstance(ind, slice):
            shifts.append(ind.start)
            scales.append(ind.step)
            shape_out.append(shape[d] // abs(ind.step))
            slicer.append(slice(None))
        elif isinstance(ind, int):
            scales.append(0)
            shifts.append(ind)
            slicer.append(0)
        else:
            slicer.append(None)
            assert (ind is None), "Strange index of type {}".format(type(ind))

    # build voxel-to-voxel transformation matrix
    lin = torch.diag(torch.as_tensor(scales, **backend))
    if any(not isinstance(s, slice) for s in slicer):
        # drop/add columns
        lin = torch.unbind(lin, dim=-1)
        zero = torch.zeros(len(shifts), **backend)
        new_lin = []
        for s in slicer:
            if isinstance(s, slice):
                col, *lin = lin
                new_lin.append(col)
            elif isinstance(s, int):
                col, *lin = lin
            elif s is None:
                new_lin.append(zero)
        lin = torch.stack(new_lin, dim=-1) if new_lin else []
    trl = torch.as_tensor(shifts, **backend)[..., None]
    trf34 = torch.cat((lin, trl), dim=1) if len(lin) else trl
    trf = torch.eye(4, **backend)
    trf[:-1, :] = trf34

    # compose
    affine = affine.matmul(trf)
    return affine, tuple(shape_out)


####

def make_kernels1d():
    F0 = 151 / 315
    F1 = 397 / 1680
    F2 = 1 / 42
    F3 = 1 / 5040
    G0 = 2 / 3
    G1 = -1 / 8
    G2 = -1 / 5
    G3 = -1 / 120
    H0 = 8 / 3
    H1 = -3 / 2
    H2 = 0
    H3 = 1 / 6
    FG0 = 0
    FG1 = -49/144
    FG2 = -7/90
    FG3 = -1/720
    F = [F3, F2, F1, F0, F1, F2, F3]
    G = [G3, G2, G1, G0, G1, G2, G3]
    H = [H3, H2, H1, H0, H1, H2, H3]
    FG = [-FG3, -FG2, -FG1, FG0, FG1, FG2, FG3]
    F = torch.as_tensor(F, dtype=torch.double)
    G = torch.as_tensor(G, dtype=torch.double)
    H = torch.as_tensor(H, dtype=torch.double)
    FG = torch.as_tensor(FG, dtype=torch.double)
    return F, G, H, FG


def make_absolute3_kernel():
    F, *_ = make_kernels1d()
    K = F[None, None, :] * F[None, :, None] * F[:, None, None]
    return K


def make_membrane3_kernel():
    F, G, *_ = make_kernels1d()
    K = (F[None, None, :] * F[None, :, None] * G[:, None, None] +
         F[None, None, :] * G[None, :, None] * F[:, None, None] +
         G[None, None, :] * F[None, :, None] * F[:, None, None])
    return K


def make_bending3_kernel():
    F, G, H, *_ = make_kernels1d()
    K = (F[None, None, :] * F[None, :, None] * H[:, None, None] +
         F[None, None, :] * H[None, :, None] * F[:, None, None] +
         H[None, None, :] * F[None, :, None] * F[:, None, None] +
         F[None, None, :] * G[None, :, None] * G[:, None, None] * 2 +
         G[None, None, :] * F[None, :, None] * G[:, None, None] * 2 +
         G[None, None, :] * G[None, :, None] * F[:, None, None] * 2)
    return K


def make_linearelastic3_kernel():
    FF, GG, HH, FG = make_kernels1d()
    # diagonal of lam (divergence)
    Kxx = GG[None, None, :] * FF[None, :, None] * FF[:, None, None]
    Kyy = FF[None, None, :] * GG[None, :, None] * FF[:, None, None]
    Kzz = FF[None, None, :] * FF[None, :, None] * GG[:, None, None]
    # off diagonal (common to lam and mu)
    Kxy = - FG[None, None, :] * FG[None, :, None] * FF[:, None, None]
    Kxz = - FG[None, None, :] * FF[None, :, None] * FG[:, None, None]
    Kyz = - FF[None, None, :] * FG[None, :, None] * FG[:, None, None]
    # diagonal of mu == membrane of each component
    return Kxx, Kyy, Kzz, Kxy, Kxz, Kyz


kernels1d = make_kernels1d()
absolute3_kernel = make_absolute3_kernel()
membrane3_kernel = make_membrane3_kernel()
bending3_kernel = make_bending3_kernel()
linearelastic3_kernel = make_linearelastic3_kernel()


# The convolution kernels are separable so could be applied as a series
# of 1D convolutions. However, I've done a quick benchmark and it does
# not look beneficial.
#
# Benchmark on a [192, 192, 192, 3] field
# GPU
#   separable = True:  230 ms
#   separable = False:  70 ms
# CPU (2 x 20 cores)
#   separable = True:  600 ms
#   separable = False: 400 ms
# so better to do one 7x7x7 convolution


def absolute3(x, bound='circular'):
    """Absolute energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * membrane3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, absolute3_kernel, bound)


def membrane3(x, bound='circular'):
    """Membrane energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * membrane3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, membrane3_kernel, bound)


def bending3(x, bound='circular'):
    """Bending energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * bending3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    return _conv3(x, bending3_kernel, bound)


def linearelastic3(x, mu=0.05, lam=0.2, bound='circular'):
    """Linear-elastic energy of a field encoded by cubic splines.

    Apply the forward matrix-vector product of the regularization: L @ x
    The full loss is computed by: loss = 0.5 * (x * linearelastic3(x)).mean()

    Assumes isotropic voxel size.

    Parameters
    ----------
    x : (nx, ny, nz, 3) tensor
        Spline coefficients of a displacement field, in voxels
    mu : float
        Second lame constant (penalty on shears)
    lam : float
        First lame constant (penalty on divergence)
    bound : {'circular', 'reflect', 'zeros', 'replicate'}
        Boundary conditions

    """
    Kxx, Kyy, Kzz, Kxy, Kxz, Kyz = linearelastic3_kernel
    M = membrane3_kernel.to(x)
    K = M.new_empty([3, 3, 7, 7, 7])
    K[0, 0] = lam * Kxx + mu * M
    K[1, 1] = lam * Kyy + mu * M
    K[2, 2] = lam * Kzz + mu * M
    K[0, 1] = K[1, 0] = (lam + mu) * Kxy
    K[0, 2] = K[2, 0] = (lam + mu) * Kxz
    K[1, 2] = K[2, 1] = (lam + mu) * Kyz
    y = _conv3(x, K, bound)
    return y


def _conv3(x, kernel, bound='circular'):
    kernel = kernel.to(x)
    if kernel.ndim == 5:
        # linear elastic -> (3, 3, 7, 7, 7) kernel
        x = x.movedim(-1, 0)[None]
        y = functional.pad(x, [3]*6, mode=bound)
        y = functional.conv3d(y, kernel)
        y= y[0].movedim(0, -1)
    else:
        # absolute/membrane/bending -> (7, 7, 7) kernel
        x = x[None].movedim(-1, 0)
        y = functional.pad(x, [3]*6, mode=bound)
        y = functional.conv3d(y, kernel[None, None])
        y = y.movedim(0, -1)[0]
    return y


# Jacobian determinant in Torch
def jacobian_det_torch(phi, spacing=(1.0, 1.0, 1.0)):
    sx, sy, sz = spacing
    ux, uy, uz = phi.unbind(-1)
    def diff(f, dim, h):
        d = torch.zeros_like(f)
        slc = [slice(None)] * 3
        slc[dim] = slice(0, 1)
        slc_f = slc.copy()
        slc_f[dim] = slice(1, 2)
        d[tuple(slc)] = (f[tuple(slc_f)] - f[tuple(slc)]) / h
        slc[dim] = slice(-1, None)
        slc_f[dim] = slice(-2, -1)
        d[tuple(slc)] = (f[tuple(slc)] - f[tuple(slc_f)]) / h
        slc[dim] = slice(1, -1)
        slc_f[dim] = slice(2, None)
        slc_b = slc.copy()
        slc_b[dim] = slice(None, -2)
        d[tuple(slc)] = 0.5 * (f[tuple(slc_f)] - f[tuple(slc_b)]) / h
        return d

    ux_x = diff(ux, 0, sx); ux_y = diff(ux, 1, sy); ux_z = diff(ux, 2, sz)
    uy_x = diff(uy, 0, sx); uy_y = diff(uy, 1, sy); uy_z = diff(uy, 2, sz)
    uz_x = diff(uz, 0, sx); uz_y = diff(uz, 1, sy); uz_z = diff(uz, 2, sz)
    return (ux_x * (uy_y * uz_z - uy_z * uz_y) - ux_y * (uy_x * uz_z - uy_z * uz_x) + ux_z * (uy_x * uz_y - uy_y * uz_x))

