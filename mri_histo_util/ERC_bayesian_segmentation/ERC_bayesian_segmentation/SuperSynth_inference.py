import torch
import numpy as np
from ext.unet3d.model import EugeniosResidualEncoderUNet3D
from ext.my_functions import get_largest_connected_component
from torch.nn import Softmax
from scipy.ndimage.morphology import binary_dilation
import csv

def run_inference(im, flipping, model_file, mode, volfile, device, force_tiling=False):

    # some constants
    f_maps = 96
    tile_size = 160
    label_list_segmentation_whole_freesurfer = [0, 14, 15, 16, 24, 77, 85, 99, 901, 902, 906, 907, 908, 909, 911, 912, 914, 915, 916,
                                                930, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869,
                                                41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844, 866, 870]
    label_list_segmentation_exvivo_freesurfer = [0, 14, 15, 16, 77, 85, 99, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26,
                                                 819, 821, 843, 865, 869, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58,
                                                 820, 822, 844, 866, 870]
    label_list_segmentation_cerebrum_freesurfer = [0,  77,  85,  99, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869,
                                                   41,  42,  43,  44,  49,  50,  51,  52,  53,  54,  58, 820, 822, 844, 866, 870]
    label_list_segmentation_hemi_freesurfer_left = [0, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 77, 99, 819, 821, 843, 865, 869]
    label_list_segmentation_hemi_freesurfer_right = [0, 41, 42, 43, 44, 49, 50, 51, 52, 54, 54, 58, 77, 99, 820, 822, 844, 866, 870]
    label_list_segmentation_whole = [0, 11, 12, 13, 16, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46,
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 47, 49, 51, 53, 55,
                                     18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 48, 50, 52, 54, 56]
    label_list_segmentation_hemis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    label_list_segmentation_exvivo = [0, 11, 12, 13, 31, 32, 33, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 34, 36, 38,
                                      40, 42, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 37, 39, 41, 43]
    n_neutral_labels_whole = 20
    n_neutral_labels_hemis = len(label_list_segmentation_hemis)
    n_neutral_labels_exvivo = 7
    n_neutral_labels_cerebrum = 4
    n_labels_whole = len(label_list_segmentation_whole)
    n_labels_hemis = len(label_list_segmentation_hemis)
    n_labels_exvivo = len(label_list_segmentation_exvivo)
    n_labels_cerebrum = len(label_list_segmentation_cerebrum_freesurfer)
    nlat = int((n_labels_whole - n_neutral_labels_whole) / 2.0)
    vflip_invivo = np.concatenate([np.array(range(n_neutral_labels_whole)),
                                   np.array(range(n_neutral_labels_whole + nlat, n_labels_whole)),
                                   np.array(range(n_neutral_labels_whole, n_neutral_labels_whole + nlat))])

    nlat = int((len(label_list_segmentation_exvivo) - n_neutral_labels_exvivo) / 2.0)
    vflip_exvivo = np.concatenate([np.array(range(n_neutral_labels_exvivo)),
                                   np.array(range(n_neutral_labels_exvivo + nlat, len(label_list_segmentation_exvivo))),
                                   np.array(range(n_neutral_labels_exvivo, n_neutral_labels_exvivo + nlat))])
    nlat = int((len(label_list_segmentation_cerebrum_freesurfer) - n_neutral_labels_cerebrum) / 2.0)
    vflip_cerebrum = np.concatenate([np.array(range(n_neutral_labels_cerebrum)),
                                   np.array(range(n_neutral_labels_cerebrum + nlat, len(label_list_segmentation_cerebrum_freesurfer))),
                                   np.array(range(n_neutral_labels_cerebrum, n_neutral_labels_cerebrum + nlat))])
    final_layers = ['reg', 'seg']
    final_layer_nf = [3, n_labels_whole]

    # down to business
    with torch.no_grad():

        # Some more variables that we put in the GPU
        list_to_kill_photo_whole = [5, 6, 11, 12, 13, 16, 22, 23, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46]
        mask_photo_or_cerebrum_whole = torch.ones(len(label_list_segmentation_whole), dtype=torch.bool, device=device)
        for l in range(len(label_list_segmentation_whole)):
            if np.sum(np.array(list_to_kill_photo_whole) == label_list_segmentation_whole[l]) > 0:
                mask_photo_or_cerebrum_whole[l] = False
        v_left = []
        for lab in label_list_segmentation_hemi_freesurfer_left:
            v_left.append(np.where(np.array(label_list_segmentation_whole_freesurfer) == lab)[0][0])
        v_left = torch.tensor(v_left, device=device, dtype=torch.int32)
        v_right = []
        for lab in label_list_segmentation_hemi_freesurfer_right:
            v_right.append(np.where(np.array(label_list_segmentation_whole_freesurfer) == lab)[0][0])
        v_right = torch.tensor(v_right, device=device, dtype=torch.int32)
        mask_exvivo_whole = torch.ones(len(label_list_segmentation_whole), dtype=torch.bool, device=device)
        for l in range(len(label_list_segmentation_whole)):
            if np.sum(np.array(label_list_segmentation_exvivo_freesurfer) == label_list_segmentation_whole_freesurfer[l]) == 0:
                mask_exvivo_whole[l] = False


        print('Preparing model and loading weights')
        if device.type == 'cpu':
            cp = torch.load(model_file, map_location=torch.device('cpu'))
        else:
            cp = torch.load(model_file)
        backbone = EugeniosResidualEncoderUNet3D(1, None, final_sigmoid=False, f_maps=f_maps, layer_order='cgl',
                        num_groups=8, num_levels=5, is_segmentation=False, is3d=True,
                        skip_final_convolution=True).to(device)
        backbone.load_state_dict(cp['backbone_state_dict'])
        final_convs = dict()
        for l in range(len(final_layers)):
            final_convs[final_layers[l]] = torch.nn.Conv3d(f_maps, final_layer_nf[l], 1, device=device)
            final_convs[final_layers[l]].load_state_dict(cp[final_layers[l] + '_state_dict'])


        print('SuperSynth is working on the input image')
        print('   Padding input image')
        im /= im.max() # in case
        W = (np.ceil(np.array(im.shape) / 32.0) * 32).astype('int')
        if tile_size is not None:
            W[W < tile_size] = tile_size
        idx = np.floor((W - im.shape) / 2).astype('int')
        S = torch.zeros(*W, dtype=torch.float32, device=device)
        S[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]] = im

        print('   Pushing data through the CNN')
        if (device.type=='cpu') and (force_tiling==False):
            print('   Working on CPU; inference without tiling')
            bb = backbone(S[None, None, ...])
        else:
            print('   Working on ' + str(device) + '; inference with tiling')
            bb = process_tile(backbone, S[None, None, ...], tile_size=tile_size)
        reg = torch.permute(final_convs['reg'](bb)[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]], [1, 2, 3, 0])
        activations = final_convs['seg'](bb)[0, :, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
        if mode=='cerebrum':
            activations = activations[mask_photo_or_cerebrum_whole, ...]
        elif mode=='left-hemi':
            activations = activations[v_left, ...]
        elif mode=='right-hemi':
            activations = activations[v_right, ...]
        elif mode=='exvivo':
            activations = activations[mask_exvivo_whole, ...]
        elif mode=='invivo':
            pass
        else:
            raise Exception('mode not supported: ' + mode)
        softmax = Softmax(dim=0)
        seg = softmax(activations)
        if flipping:
            S = torch.flip(S, [0])
            if (device.type=='cpu') and (force_tiling==False):
                bb = backbone(S[None, None, ...])
            else:
                bb = process_tile(backbone, S[None, None, ...], tile_size=tile_size)
            aux = torch.flip(torch.permute(final_convs['reg'](bb)[0, ...], [1, 2, 3, 0]), [0])
            aux[..., 0] = -aux[..., 0]
            reg = 0.5 * reg + 0.5 * aux[idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2], :]
            activations = torch.flip(final_convs['seg'](bb)[0, ...], [1])[:, idx[0]:idx[0] + im.shape[0], idx[1]:idx[1] + im.shape[1], idx[2]:idx[2] + im.shape[2]]
            if mode == 'cerebrum':
                activations = activations[mask_photo_or_cerebrum_whole, ...]
                activations = activations[vflip_cerebrum, ...]
            elif mode == 'left-hemi':
                activations = activations[v_right, ...]
            elif mode == 'right-hemi':
                activations = activations[v_left, ...]
            elif mode == 'exvivo':
                activations = activations[mask_exvivo_whole, ...]
                activations = activations[vflip_exvivo, ...]
            else: # 'invivo':
                activations = activations[vflip_invivo, ...]
            seg = 0.5 * seg + 0.5 * softmax(activations)

        if mode=='cerebrum':
            seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[mask_photo_or_cerebrum_whole][torch.argmax(seg, 0)]
        elif mode=='left-hemi':
            seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_left, device=device)[torch.argmax(seg, 0)]
        elif mode=='right-hemi':
            seg_discrete = torch.tensor(label_list_segmentation_hemi_freesurfer_right, device=device)[torch.argmax(seg, 0)]
        elif mode=='exvivo':
            seg_discrete = torch.tensor(label_list_segmentation_exvivo_freesurfer, device=device)[torch.argmax(seg, 0)]
        elif mode=='invivo':
            seg_discrete = torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[torch.argmax(seg, 0)]
        else:
            raise Exception('mode not supported: ' + mode)

        # get masks for postprocessing segmentations / volumes
        M = (seg_discrete > 0) & (seg_discrete != 24) & (seg_discrete != 99) & (seg_discrete < 900) # useful for later
        M = get_largest_connected_component(M.detach().cpu().numpy())
        Mdilated = binary_dilation(M, iterations=2)
        M = torch.tensor(M, device=device, dtype=torch.bool)
        Mdilated = torch.tensor(Mdilated, device=device, dtype=torch.bool)
        seg_discrete[~M] = 0

        # postprocess soft segmentations and compute volumes
        seg[0][~Mdilated] = 1
        for l in range(seg.shape[0]):
            seg[l][~Mdilated] = 0
        vols = seg.sum(dim=[1, 2, 3]).detach().cpu().numpy()
        with open(volfile, 'w') as csvfile:
            writer = csv.writer(csvfile)
            if mode == 'cerebrum':
                llist = (torch.tensor(label_list_segmentation_whole_freesurfer, device=device)[mask_photo_or_cerebrum_whole]).detach().cpu().numpy()
            elif mode == 'left-hemi':
                llist = label_list_segmentation_hemi_freesurfer_left
            elif mode == 'right-hemi':
                llist = label_list_segmentation_hemi_freesurfer_right
            elif mode == 'exvivo':
                llist = label_list_segmentation_exvivo_freesurfer
            elif mode == 'invivo':
                llist = label_list_segmentation_whole_freesurfer
            row1 = []
            row2 = []
            for l in range(len(llist)):
                lab = llist[l]
                if (lab > 0) and ((lab != 24) or (mode=='invivo')) and (lab < 900) and (lab != 99):
                    row1.append(str(lab))
                    row2.append(str(vols[l]))
            writer.writerow(row1)
            writer.writerow(row2)

    return seg_discrete, reg

def process_tile(model, input_tensor, tile_size=160):
    B, C, H, W, D = input_tensor.shape
    weight_map = torch.zeros_like(input_tensor)
    # Blending mask
    x = y = z = torch.linspace(-1, 1, tile_size, device=input_tensor.device)
    xx, yy, zz = torch.meshgrid(x, y, z, indexing="ij")
    blending_mask = torch.cos(xx * torch.pi / 2) * torch.cos(yy * torch.pi / 2) * torch.cos(zz * torch.pi / 2)
    blending_mask = blending_mask[None, None, :, :].expand(1, C, -1, -1, -1)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                i1 = 0 if i == 0 else (H - tile_size)
                i2 = tile_size if i == 0 else H
                j1 = 0 if j == 0 else (W - tile_size)
                j2 = tile_size if j == 0 else W
                k1 = 0 if k == 0 else (D - tile_size)
                k2 = tile_size if k == 0 else D

                tile = input_tensor[:, :, i1:i2, j1:j2, k1:k2]
                with torch.no_grad():
                    tile_output = model(tile)
                if (i == 0) and (j == 0) and (k == 0):
                    output = torch.zeros([B, tile_output.shape[1], H, W, D], device=input_tensor.device)
                output[:, :, i1:i2, j1:j2, k1:k2] += tile_output * blending_mask
                weight_map[:, :, i1:i2, j1:j2, k1:k2] += blending_mask

    output /= torch.clamp(weight_map, min=1e-8)
    return output
