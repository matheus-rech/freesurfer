import nibabel as nib
import torch
import numpy as np
import os
from scipy.ndimage import label


###############################
def get_largest_connected_component(binary_numpy):
    labeled_array, num_features = label(binary_numpy)
    if num_features == 0:
      return np.zeros_like(binary_numpy)
    component_sizes = np.bincount(labeled_array.flatten())
    largest_component_label = np.argmax(component_sizes[1:]) + 1
    largest_component_mask = (labeled_array == largest_component_label)
    return largest_component_mask
###############################

def MRIwrite(volume, aff, filename, dtype=None):

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if aff is None:
        aff = np.eye(4)
    header = nib.Nifti1Header()
    nifty = nib.Nifti1Image(volume, aff, header)

    nib.save(nifty, filename)

###############################

def MRIread(filename, dtype=None, im_only=False):

    assert filename.endswith(('.nii', '.nii.gz', '.mgz')), 'Unknown data file: %s' % filename

    x = nib.load(filename)
    volume = x.get_fdata()
    aff = x.affine

    if dtype is not None:
        volume = volume.astype(dtype=dtype)

    if im_only:
        return volume
    else:
        return volume, aff

############################

def myzoom_torch_anisotropic(X, aff, newsize, device):

    if len(X.shape)==3:
        X = X[..., None]

    dtype = X.dtype

    factors = np.array(newsize) / np.array(X.shape[:-1])
    delta = (1.0 - factors) / (2.0 * factors)

    vx = torch.arange(delta[0], delta[0] + newsize[0] / factors[0], 1 / factors[0], dtype=dtype, device=device)[:newsize[0]]
    vy = torch.arange(delta[1], delta[1] + newsize[1] / factors[1], 1 / factors[1], dtype=dtype, device=device)[:newsize[1]]
    vz = torch.arange(delta[2], delta[2] + newsize[2] / factors[2], 1 / factors[2], dtype=dtype, device=device)[:newsize[2]]

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

    Y = torch.zeros([newsize[0], newsize[1], newsize[2], X.shape[3]], dtype=dtype, device=device)

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

############################

def torch_resize(I, aff, resolution, device, power_factor_at_half_width=5, dtype=torch.float32, slow=False):

    if torch.is_grad_enabled():
        with torch.no_grad():
            return torch_resize(I, aff, resolution, device, power_factor_at_half_width, dtype, slow)

    slow = slow or (device == 'cpu')
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
            It = It.permute([0, 1, 3, 4, 2])
            if sigmas[d]>0:
                sl = np.ceil(sigmas[d] * 2.5).astype(int)
                v = np.arange(-sl, sl + 1)
                gauss = np.exp((-(v / sigmas[d]) ** 2 / 2))
                kernel = gauss / np.sum(gauss)
                kernel = torch.tensor(kernel,  device=device, dtype=dtype)
                if slow:
                    It = conv_slow_fallback(It, kernel)
                else:
                    kernel = kernel[None, None, None, None, :]
                    It = torch.conv3d(It, kernel, bias=None, stride=1, padding=[0, 0, int((kernel.shape[-1] - 1) / 2)])


        It = torch.squeeze(It)
        It, aff2 = myzoom_torch_anisotropic(It, aff, newsize, device)
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



#######


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

##############

def get_ras_axes(aff, n_dims=3):
    # Exhaustive search is fine, only 6 possible permutations...
    candidates = [[0,1,2], [0,2,1], [1,0,2], [2,0,1], [1,2,0], [2,1,0]]
    best_score = -100000
    best_candidate = None
    for c in candidates:
        score = np.abs(aff[0, c[0]] / np.linalg.norm(aff[:, c[0]])) \
              + np.abs(aff[1, c[1]] / np.linalg.norm(aff[:, c[1]])) \
              + np.abs(aff[2, c[2]] / np.linalg.norm(aff[:, c[2]]))
        if score>best_score:
            best_score = score
            best_candidate = c
    return np.array(best_candidate)

################

def make_gaussian_kernel(sigma, device):
    sl = int(np.ceil(3 * sigma))
    ts = torch.linspace(-sl, sl, 2*sl+1, dtype=torch.float, device=device)
    gauss = torch.exp((-(ts / sigma)**2 / 2))
    kernel = gauss / gauss.sum()
    return kernel

################

def gaussian_blur_3d(input, stds, device):
    from torch.nn.functional import conv3d
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

################

def get_label_lists_etc():
    # Just defines a bunch of constants
    label_list_segmentation_whole_freesurfer = [0, 14, 15, 16, 24, 77, 85, 99, 901, 902, 906, 907, 908, 909, 911,
                                                912, 914, 915, 916,
                                                930, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843,
                                                865, 869,
                                                41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844,
                                                866, 870]
    label_list_segmentation_exvivo_freesurfer = [0, 14, 15, 16, 77, 85, 99, 2, 3, 4, 5, 7, 8, 10, 11, 12, 13, 17,
                                                 18, 26,
                                                 819, 821, 843, 865, 869, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52,
                                                 53, 54, 58,
                                                 820, 822, 844, 866, 870]
    label_list_segmentation_cerebrum_freesurfer = [0, 77, 85, 99, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 819, 821,
                                                   843, 865, 869,
                                                   41, 42, 43, 44, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844, 866,
                                                   870]
    label_list_segmentation_hemi_freesurfer_left = [0, 2, 3, 4, 5, 10, 11, 12, 13, 17, 18, 26, 77, 99, 819, 821,
                                                    843, 865, 869]
    label_list_segmentation_hemi_freesurfer_right = [0, 41, 42, 43, 44, 49, 50, 51, 52, 53, 54, 58, 77, 99, 820,
                                                     822, 844, 866, 870]
    label_list_segmentation_whole = [0, 11, 12, 13, 16, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46,
                                     1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 47, 49, 51, 53, 55,
                                     18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 48, 50, 52, 54, 56]
    label_list_segmentation_hemis = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18]
    label_list_segmentation_exvivo = [0, 11, 12, 13, 31, 32, 33, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 14, 15, 17, 34, 36,
                                      38,
                                      40, 42, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 35, 37, 39, 41,
                                      43]
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
                                   np.array(
                                       range(n_neutral_labels_exvivo + nlat, len(label_list_segmentation_exvivo))),
                                   np.array(range(n_neutral_labels_exvivo, n_neutral_labels_exvivo + nlat))])
    nlat = int((len(label_list_segmentation_cerebrum_freesurfer) - n_neutral_labels_cerebrum) / 2.0)
    vflip_cerebrum = np.concatenate([np.array(range(n_neutral_labels_cerebrum)),
                                     np.array(range(n_neutral_labels_cerebrum + nlat,
                                                    len(label_list_segmentation_cerebrum_freesurfer))),
                                     np.array(range(n_neutral_labels_cerebrum, n_neutral_labels_cerebrum + nlat))])
    list_to_kill_photo_whole = [5, 6, 11, 12, 13, 16, 22, 23, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 46]

    return label_list_segmentation_whole_freesurfer, label_list_segmentation_exvivo_freesurfer, label_list_segmentation_cerebrum_freesurfer, \
           label_list_segmentation_hemi_freesurfer_left, label_list_segmentation_hemi_freesurfer_right, label_list_segmentation_whole, \
           label_list_segmentation_hemis, label_list_segmentation_exvivo, n_neutral_labels_whole, n_neutral_labels_hemis, n_neutral_labels_exvivo, \
           n_neutral_labels_cerebrum, n_labels_whole, n_labels_hemis, n_labels_exvivo, n_labels_cerebrum, vflip_invivo, vflip_exvivo, vflip_cerebrum, \
           list_to_kill_photo_whole
