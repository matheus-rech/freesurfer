#Imports
import os
import sys
BASE_PATH = os.path.join(os.environ.get('FREESURFER_HOME'),'python/packages/ERC_bayesian_segmentation/')
sys.path.insert(0, BASE_PATH)
import numpy as np
import torch
from torch.utils.data import Dataset, DataLoader
import ext.my_functions as my
from datetime import datetime
import ERC_bayesian_segmentation.relabeling as relab
import glob
import scipy.sparse as sp
import ext.bias_field_correction_torch as bf
import csv
import argparse
import math
from torch.nn.functional import grid_sample
from ext.fireants_trimmed import Image, BatchedImages
from ext.fireants_trimmed import GreedyRegistration
import SimpleITK as sitk
from time import time
from ext.fireants_trimmed import HybridDiceLabelDiffloss
from scipy.ndimage import binary_erosion, binary_fill_holes


########################################################

parser = argparse.ArgumentParser(description='NextBrain segmentation with SuperSynth and FireANTs.')
parser.add_argument("--i", help="Image to segment.", required=True)
parser.add_argument("--model_file", help="Multi-task deep learning model for preprocessing", required=True)
parser.add_argument("--atlas_dir", help="Atlas directory", required=True)
parser.add_argument("--o", help="Output directory.", required=True)
parser.add_argument("--mode", help="Type of input (invivo, exvivo, cerebrum, hemi).", required=True)
parser.add_argument("--side", help="Hemisphere to segment (left or right).", required=True)
parser.add_argument("--bf_mode", help="bias field basis function: dct, polynomial, or hybrid", default="dct")
parser.add_argument("--yaml_path", help="path of custom YAML files to define groups of ROIs", default=None)
parser.add_argument("--write_rgb", action="store_true", help="Write soft segmentation to dis as RGB file.")
parser.add_argument("--write_bias_corrected", action="store_true", help="Write bias field corrected image to disk")
parser.add_argument("--device", help="Device (cpu, cuda)")
parser.add_argument("--device_registration", help="Use this option if you want to use a different device for the registration (useful for high-res ex vivo)")
parser.add_argument("--threads", type=int, default=-1, help="(optional) Number of CPU cores to be used. Default is -1 (use all available cores")
parser.add_argument("--skip", type=int, default=1, help="(optional) Skipping factor to easy memory requirements of priors when estimating Gaussian parameters. Default is 1.")
parser.add_argument("--resolution", type=float, default=0.4, help="(optional) Resolution of output segmentation")
parser.add_argument("--smoothing_steps_HRmask", type=int, default=3, help="(optional) Number of smoothing iterations when upsampling mask from 1mm segmentation")
parser.add_argument("--skip_bf", action="store_true", help="Skip bias field correction altogether")
parser.add_argument("--smooth_grad_sigma", type=float, default=1.00, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--smooth_warp_sigma", type=float, default=0.25, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--optimizer_lr", type=float, default=0.5, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--cc_kernel_size", type=int, default=7, help="(optional) Parameter of Greedy FireANTs registration")
parser.add_argument("--rel_weight_labeldiff", type=float, default=2.5, help="(optional) Relative weight of labels in Greedy FireANTs registration")
parser.add_argument("--save_atlas_nonlinear_reg", action="store_true", help="Save nonlinearly registered atlas")
parser.add_argument("--save_field", action="store_true", help="Save nonlinear deformation field")
parser.add_argument("--save_jacobian", action="store_true", help="Save Jacobian determinant (log10)")
args = parser.parse_args()

########################################################

if args.i is None:
    raise Exception('Input image is required')
if args.model_file is None:
    raise Exception('Model file is required')
if args.atlas_dir is None:
    raise Exception('Atlas directory must be provided')
if args.o is None:
    raise Exception('Output directory must be provided')
if args.side is None:
    raise Exception('side must be provided (left or right)')
else:
    if (args.side!='left') and (args.side!='right'):
        raise Exception('Side must be left or right, but you specified: ' + args.side)
if args.mode is None:
    raise Exception('mode must be provided (invivo/cerebrum/hemi/exvivo)')
else:
    if (args.mode!='invivo') and (args.mode!='cerebrum') and (args.mode!='hemi') and (args.mode!='exvivo'):
        raise Exception('Mode must be invivo/cerebrum/hemi/exvivo, but you specified: ' + args.mode)
if (args.resolution<=0):
    raise Exception('Resolution must be non-negative; Exitting...')
if (args.skip<1):
    raise Exception('Skip cannot be less than 1; exitting...')


########################################################

if torch.cuda.is_available():
    print('GPU/CUDA seems to be available')
    if args.device is None:
        args.device = 'cpu'
    if args.device_registration is None:
        args.device_registration = args.device
    print('You selected use of the following devices:')
    print('  Registration: ' + args.device_registration)
    device_registration = torch.device(args.device_registration)
    print('  Rest of code: ' + args.device)
    device = torch.device(args.device)
    os.environ["PYTORCH_CUDA_ALLOC_CONF"] = "garbage_collection_threshold:0.5,max_split_size_mb:32"

else:
    print('GPU/CUDA not available; using the CPU')
    device = torch.device('cpu')
    device_registration = torch.device('cpu')

########################################################

# limit the number of CPU threads to be used
if args.threads<0:
    args.threads = os.cpu_count()
    print('Using all available CPU threads ( %s )' % args.threads)
else:
    print('Using %s CPU thread(s)' % args.threads)
torch.set_num_threads(args.threads)

########################################################

# Reproducibility: https://discuss.pytorch.org/t/reproducibility-with-all-the-bells-and-whistles/81097
seed = 0; torch.manual_seed(seed); torch.cuda.manual_seed_all(seed); torch.cuda.manual_seed(seed)
np.random.seed(seed); np.random.seed(seed)
torch.backends.cudnn.deterministic = True; torch.backends.cudnn.benchmark = False

########################################################

# Disable gradients
torch.set_grad_enabled(False)

################

# Input data
input_volume = args.i
model_file = args.model_file
atlas_dir = args.atlas_dir
LUT_file = os.path.join(BASE_PATH, 'data_simplified', 'AllenAtlasLUT')
output_dir = args.o
skip_bf = args.skip_bf
bf_mode = args.bf_mode
side = args.side
mode = args.mode
resolution = args.resolution
skip = args.skip

########################################################
# Detect problems with output directory right off the bat
if os.path.exists(output_dir):
    if len(os.listdir(output_dir))==0:
        print('Warning: output directory exists')
    else:
        pass
        #raise Exception('Ouput directory exists and is not empty; exitting...')
else:
    os.mkdir(output_dir)


########################################################

# Constants
dtype = torch.float32
SET_BG_TO_CSF = True # True = median of ventricles -> it seems much better than 0!
RESOLUTION_ATLAS = 0.2
TOL = 1e-9

########################################################

now = datetime.now()
current_time = now.strftime("%H:%M:%S")
print("Current Time =", current_time)

########################################################

supersynth_segmentation = output_dir + '/SuperSynth/segmentation.mgz'
if os.path.exists(supersynth_segmentation):
    print('No need to run SuperSynth; segmentation found:')
    print('   ' + supersynth_segmentation)
else:
    print('Analyzing image with SuperSynth...')
    mode_supersynth = (side + '-hemi') if mode=='hemi' else mode
    cmd = 'mri_super_synth --i ' + input_volume + ' --o ' + output_dir + '/SuperSynth/ --mode ' + \
           mode_supersynth + ' --threads ' + str(args.threads) + ' --device ' + args.device
    if os.system(cmd):
        raise Exception('Problem with SuperSynth; exitting...')

########################################################

print('Reading input image and SuperSynth outputs')
Iim, aff = my.MRIread(input_volume)
Iim = np.squeeze(Iim)
while len(Iim.shape) > 3:
    Iim = np.mean(Iim, axis=-1)
Iim = torch.tensor(Iim, dtype=dtype, device=device)
im, imaff = my.MRIread(output_dir + '/SuperSynth/input_resampled.mgz')
seg, _ = my.MRIread(output_dir + '/SuperSynth/segmentation.mgz')
reg, _ = my.MRIread(output_dir + '/SuperSynth/mni_deformation.mgz')
im = torch.tensor(im, dtype=dtype, device=device)
seg = torch.tensor(seg, dtype=torch.int, device=device)
reg = torch.tensor(reg, dtype=dtype, device=device)
im[torch.isnan(im)] = 0
seg[torch.isnan(seg)] = 0
reg[torch.isnan(reg)] = 0

########################################################
print('Reslice segmentation and coordinate map to the space of the input image')
II, JJ, KK = torch.meshgrid(torch.arange(Iim.shape[0], device=device),
                          torch.arange(Iim.shape[1], device=device),
                          torch.arange(Iim.shape[2], device=device), indexing='ij')
affine = np.linalg.inv(imaff) @ aff
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
Sim = my.fast_3D_interp_torch(seg, II2, JJ2, KK2, 'nearest')
LRmap = my.fast_3D_interp_torch(reg[...,0], II2, JJ2, KK2, 'linear')
APmap = my.fast_3D_interp_torch(reg[...,1], II2, JJ2, KK2, 'linear')
ISmap = my.fast_3D_interp_torch(reg[...,2], II2, JJ2, KK2, 'linear')
del II, JJ, KK, II2, JJ2, KK2, reg, seg

########################################################
print('Adding affine component to coordinate map')
with open(output_dir + '/SuperSynth/mni_affine.txt') as f:
    matrix = [list(map(float, line.split(','))) for line in f]
M_input_vox_to_mni_ras = matrix @ aff
M = M_input_vox_to_mni_ras.copy()
II, JJ, KK = torch.meshgrid(torch.arange(Iim.shape[0], device=device),
                          torch.arange(Iim.shape[1], device=device),
                          torch.arange(Iim.shape[2], device=device), indexing='ij')
LRmap += (M[0, 0] * II + M[0, 1] * JJ + M[0, 2] * KK + M[0, 3])
APmap += (M[1, 0] * II + M[1, 1] * JJ + M[1, 2] * KK + M[1, 3])
ISmap += (M[2, 0] * II + M[2, 1] * JJ + M[2, 2] * KK + M[2, 3])
del II, JJ, KK, M

########################################################
print('Killing labels as needed')
if mode=='hemi':
    labels_to_kill = [24, 99]
else:
    if side=='left':
        labels_to_kill =  [14, 15, 24, 99, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 820, 822, 844, 866, 870, 901, 902, 906, 907, 908, 909, 911, 912, 914, 915, 916, 930]
        Sim[LRmap >= 0] = 0
    else:
        labels_to_kill =  [14, 15, 24, 99,  2,  3,  4,  5,  7,  8, 10, 11, 12, 13, 17, 18, 26, 819, 821, 843, 865, 869, 901, 902, 906, 907, 908, 909, 911, 912, 914, 915, 916, 930]
        Sim[LRmap < 0] = 0
lut = np.arange(3000)
lut[labels_to_kill] = 0
lut = torch.tensor(lut, dtype=torch.int, device=device)
Sim = lut[Sim]

########################################################

if skip_bf==False:
    print('Correcting bias field')
    try:
        Iim, _, bfmask = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=6, device=device, dtype=dtype, basis=bf_mode, dontmask=True, cerebrum_only=((mode=='hemi') or (mode=='cerebrum')))
    except:
        if device.type=='cpu':
            raise Exception('Bias correction failed (out of memory?)')
        else:
            print('Bias correction on GPU failed; trying with CPU')
            Iim, _, bfmask = bf.correct_bias(Iim, Sim, maxit=100, penalty=0.1, order=4, device='cpu', dtype=dtype, basis=bf_mode, dontmask=True, cerebrum_only=((mode=='hemi') or (mode=='cerebrum')))
    if args.write_bias_corrected:
        bfmask = binary_fill_holes(bfmask.detach().cpu().numpy())
        aux = Iim.detach().cpu().numpy()
        aux[~bfmask] = 0
        aux = (aux / np.max(aux) * 255).astype(np.uint8)
        my.MRIwrite(aux, aff, output_dir + '/bias.corrected.' + side + '.mgz', dtype=np.uint8)
        del aux
    del bfmask

print('Normalizing intensities')
Iim = Iim * 110 / torch.median(Iim[(Sim==2) | (Sim==41)])


########################################################
# Kill bottom of medulla and subdivide brainstem if needed
if (mode=='invivo') or (mode=='exvivo'):
    print('Dealing with bilateral labels: subdividing brainstem, optic chiasm, lesions');
    print('  (we also crop the bottom of the brainstem a bit)')
    LEFT = (LRmap < 0)
    Sim[(Sim == 16) & (ISmap < (-60))] = 0
    if (side == 'left'):
        Sim[(Sim == 16) & LEFT] = 161
        Sim[(Sim == 16)] = 0
        Sim[(Sim == 77) & (~LEFT)] = 0
        Sim[(Sim == 85) & (~LEFT)] = 0
    else:
        Sim[(Sim == 16) & (~LEFT)] = 162
        Sim[(Sim == 16)] = 0
        Sim[(Sim == 77) & LEFT] = 0
        Sim[(Sim == 85) & LEFT] = 0
else:
    print('No need to deal with bilateral labels')

# Release memory
del LRmap, APmap, ISmap

#######################################

# Prepare data for hemisphere at hand
print('Creating mask for tissue to segment (leave out ventricles)')
M = (Sim>0)
M[Sim==4] = 0
M[Sim==5] = 0
M[Sim==43] = 0
M[Sim==44] = 0
M[Sim==14] = 0
M[Sim==15] = 0
M = torch.tensor(my.getLargestCC(M.detach().cpu().numpy()), device=M.device, dtype=torch.bool)

# I now do this with the resampled mask later on, to avoid blurring mask edges
Mim, cropping = my.cropLabelVolTorch(M, margin=5)
del M
Sim = my.applyCropping(Sim, cropping)
Iim = my.applyCropping(Iim, cropping)
Iim_before_masking = Iim.clone()
Iim[Mim==0] = 0
aff[:3, -1] = aff[:3, -1] + aff[:-1, :-1] @ np.array([cropping[0].detach().cpu().numpy(), cropping[1].detach().cpu().numpy(), cropping[2].detach().cpu().numpy()])


########################################################


# Read atlas
print('Reading in atlas')

# ASEG labels are hard coded for each hemi, no biggie
if (side=='left'):
    aseg_label_list = np.array(
        [2, 3, 7, 8, 10, 11, 12, 13, 17, 18, 26, 28, 77, 85, 161, 819, 821, 843, 865, 869]).astype(int)
else:
    aseg_label_list = np.array(
        [41, 42, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 77, 85, 162, 820, 822, 844, 866, 870]).astype(int)

# Get the label groupings and atlas labels from the config files
if args.yaml_path is None:
    yaml_path = os.path.join(BASE_PATH, 'data_simplified')
else:
    yaml_path = args.yaml_path
tissue_index, grouping_labels, label_list, number_of_gmm_components, cheating_recipe = relab.get_tissue_settings(
            os.path.join(yaml_path, 'atlas_names_and_labels.yaml'),
            os.path.join(yaml_path, 'combined_atlas_labels_fireants.yaml'),
            os.path.join(yaml_path, 'combined_aseg_labels_new_targets.yaml'),
            os.path.join(yaml_path, 'gmm_components_fireants.yaml'),
            os.path.join(yaml_path, 'recipe_intensities_cheating_image_fireants.yaml'),
            aseg_label_list
)

tidx = tissue_index[np.where(label_list == 0)[0][0]]
if tidx>0:
    raise Exception('First tissue class must be the background')
n_tissues = np.max(tissue_index) + 1
n_labels = len(label_list)
atlas_names = sorted(glob.glob(atlas_dir + '/label_*.npz'))
atlas_size = np.load(atlas_dir + '/size.npy')
atlas_bounds_all_labels = np.load(atlas_dir + '/bounds.npy')

class LabelDataset(Dataset):

    def __init__(self, fnames):
        self.fnames = fnames

    def __len__(self):
        return len(self.fnames)

    def __getitem__(self, item):
        # print(item, self.fnames[item])
        prior = sp.load_npz(self.fnames[item])
        prior_indices = torch.as_tensor(prior.row)
        prior_values = torch.as_tensor(prior.data)
        return prior_indices, prior_values

# TODO: without this line, I get weird runtime errors...
prefetch = 4
workers = 2
prefetch_factor = max(prefetch//workers, 1)
if sys.platform == "darwin":
    label_loader = DataLoader(LabelDataset(atlas_names), num_workers=0)
else:
    label_loader = DataLoader(LabelDataset(atlas_names), num_workers=workers, prefetch_factor=prefetch_factor)

# A and A_reg are now lists, where we only keep in memory bounding boxes with mass
A = []
atlas_bounds_tissues = np.zeros([n_tissues, 6], dtype=np.int32)
for t in range(n_tissues):
    bounds = atlas_bounds_all_labels[np.array(tissue_index) == t, :]
    atlas_bounds_tissues[t, 0] = np.min(bounds[:, 0])
    atlas_bounds_tissues[t, 1] = np.max(bounds[:, 1])
    atlas_bounds_tissues[t, 2] = np.min(bounds[:, 2])
    atlas_bounds_tissues[t, 3] = np.max(bounds[:, 3])
    atlas_bounds_tissues[t, 4] = np.min(bounds[:, 4])
    atlas_bounds_tissues[t, 5] = np.max(bounds[:, 5])
    siz = [atlas_bounds_tissues[t, 1] - atlas_bounds_tissues[t, 0],
           atlas_bounds_tissues[t, 3] - atlas_bounds_tissues[t, 2],
           atlas_bounds_tissues[t, 5] - atlas_bounds_tissues[t, 4]]
    A.append(np.zeros(siz, dtype=np.float32))

label_sets_reg = relab.get_label_sets_for_label_registration(mode) # this is the atlas with tissue types
A_reg = []
atlas_bounds_segs = np.zeros([len(label_sets_reg), 6], dtype=np.int32)
for t in range(len(label_sets_reg)):
    bounds = atlas_bounds_all_labels[np.where(np.in1d(label_list, label_sets_reg[t]))[0], :]
    atlas_bounds_segs[t, 0] = np.min(bounds[:, 0])
    atlas_bounds_segs[t, 1] = np.max(bounds[:, 1])
    atlas_bounds_segs[t, 2] = np.min(bounds[:, 2])
    atlas_bounds_segs[t, 3] = np.max(bounds[:, 3])
    atlas_bounds_segs[t, 4] = np.min(bounds[:, 4])
    atlas_bounds_segs[t, 5] = np.max(bounds[:, 5])
    siz = [atlas_bounds_segs[t, 1] - atlas_bounds_segs[t, 0],
           atlas_bounds_segs[t, 3] - atlas_bounds_segs[t, 2],
           atlas_bounds_segs[t, 5] - atlas_bounds_segs[t, 4]]
    A_reg.append(np.zeros(siz, dtype=np.float32))

for n, (prior_indices, prior_values) in enumerate(label_loader):
    print('  Reading in label ' + str(n+1) + ' of ' + str(n_labels), end='\r')
    if prior_indices.numel() == 0:
        continue
    prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
    prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()
    idx = tissue_index[n]

    idx_reg = -1
    for j in range(len(label_sets_reg)):
        if np.any(label_sets_reg[j] == label_list[n]):
            idx_reg = j

    if n == 0:
        prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                        [torch.Size(atlas_size).numel()]).to_dense()
        del prior_indices, prior_values
        prior = prior.reshape(torch.Size(atlas_size)).cpu().numpy()
        A[idx] = A[idx] + prior
        if idx_reg>-1:
            A_reg[idx_reg] = A_reg[idx_reg] + prior
    else:
        prior_indices = my.ind2sub(prior_indices, atlas_size)
        min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item() + 1
        min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item() + 1
        min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item() + 1
        crop_atlas_size = [max_x - min_x, max_y - min_y, max_z - min_z]
        prior_indices[0] -= min_x
        prior_indices[1] -= min_y
        prior_indices[2] -= min_z
        prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
        crop = (slice(min_x - atlas_bounds_tissues[idx,0], max_x - atlas_bounds_tissues[idx,0]),
                slice(min_y - atlas_bounds_tissues[idx,2], max_y - atlas_bounds_tissues[idx,2]),
                slice(min_z - atlas_bounds_tissues[idx,4], max_z - atlas_bounds_tissues[idx,4]))
        A[idx][crop] = A[idx][crop] + prior.cpu().numpy()
        if idx_reg>-1:
            crop = (slice(min_x - atlas_bounds_segs[idx_reg, 0], max_x - atlas_bounds_segs[idx_reg, 0]),
                    slice(min_y - atlas_bounds_segs[idx_reg, 2], max_y - atlas_bounds_segs[idx_reg, 2]),
                    slice(min_z - atlas_bounds_segs[idx_reg, 4], max_z - atlas_bounds_segs[idx_reg, 4]))
            A_reg[idx_reg][crop] = A_reg[idx_reg][crop] + prior.cpu().numpy()
print(' ')

# We keep A in the CPU for now, to same our GPU memory for registration
# A = torch.tensor(A, dtype=dtype, device=device)
if (side=='left') or (side=='left-c') or (side=='left-ccb'):
    aff_A = np.diag([.2, .2, .2, 1])
else:
    aff_A = np.diag([-.2, .2, .2, 1])

################
# fake image
MU_CSF = 0
if (side=='left'):
    MU_WM = torch.median(Iim[Sim==2])
    MU_GM = torch.median(Iim[(Sim==3) | (Sim==17) | (Sim==18)])
    MU_WM_CEREBELLUM = torch.median(Iim[Sim==7])
    MU_GM_CEREBELLUM = torch.median(Iim[Sim==8])
    MU_CAUDATE = torch.median(Iim[Sim==11])
    MU_PUTAMEN = torch.median(Iim[Sim==12])
    MU_PALLIDUM = torch.median(Iim[Sim == 13])
else:
    MU_WM = torch.median(Iim[Sim==41])
    MU_GM = torch.median(Iim[(Sim==42) | (Sim==53) | (Sim==54)])
    MU_WM_CEREBELLUM = torch.median(Iim[Sim==46])
    MU_GM_CEREBELLUM = torch.median(Iim[Sim==47])
    MU_CAUDATE = torch.median(Iim[Sim==50])
    MU_PUTAMEN = torch.median(Iim[Sim==51])
    MU_PALLIDUM = torch.median(Iim[Sim == 52])
# Kill cerebellum and brainstem if needed
if (mode=='cerebrum') or (mode=='hemi'): # We don't kill the brainstem (hard to know where to crop) and let the registration handle it
    MU_WM_CEREBELLUM = MU_GM_CEREBELLUM = MU_DG_CEREBELLUM = 0
# Make cheating means
cheating_recipe = torch.tensor(cheating_recipe, device=device, dtype=dtype)
cheating_recipe[:,0] *= MU_WM
cheating_recipe[:,1] *= MU_GM
cheating_recipe[:,2] *= MU_WM_CEREBELLUM
cheating_recipe[:,3] *= MU_GM_CEREBELLUM
cheating_recipe[:,4] *= MU_CAUDATE
cheating_recipe[:,5] *= MU_PUTAMEN
cheating_recipe[:,6] *= MU_PALLIDUM
cheating_means = cheating_recipe.sum(dim=1)
# And make the actual image
sigma =  10.0
muI = torch.zeros(*atlas_size, device=device, dtype=dtype)
for l in range(n_tissues):
    muI[atlas_bounds_tissues[l,0]:atlas_bounds_tissues[l,1],
        atlas_bounds_tissues[l,2]:atlas_bounds_tissues[l,3],
        atlas_bounds_tissues[l,4]:atlas_bounds_tissues[l,5]] += (torch.tensor(A[l], device=device, dtype=dtype) * cheating_means[l])
sigmaI = sigma * torch.ones(muI.shape, device=device, dtype=dtype)
Ifake = torch.normal(muI, sigmaI)
del muI, sigmaI
Ifake[Ifake<0] = 0
native_resolution = np.sqrt(np.sum(aff[:-1,:-1]**2, axis=0))
# we blur a bit less since the atlas is blurry already (default power factor at W/2 is 5.0)
Ifake, aff_fake = my.torch_resize(Ifake, aff_A, native_resolution, device, dtype=dtype, power_factor_at_half_width=(0.5*5.0))

# For the linear registration, we concatenate the subject-MNI transform with precomputed MNI-NextBrain transforms
if (side=='left') or (side=='left-c') or (side=='left-ccb'):
    Mmni = np.matrix('1.0959   -0.0131   -0.0233   79.9750;    0.0381    1.1359    0.0035  136.9361;  0.0735   -0.0307    1.0872   89.2956; 0 0 0 1')
else:
    Mmni = np.matrix('1.0959    0.0131    0.0233  -79.9750;   -0.0381    1.1359    0.0035  136.9361; -0.0735   -0.0307    1.0872   89.2956; 0 0 0 1')
shift_mat = np.eye(4); shift_mat[0,3]=cropping[0]; shift_mat[1,3]=cropping[1]; shift_mat[2,3]=cropping[2] # accounts for the cropping!
affine = np.linalg.inv(aff_fake) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat
# Resampling
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
RSlin = my.fast_3D_interp_torch(Ifake, II2, JJ2, KK2, 'linear')
fake_filename = output_dir + '/temp_fake_linear.nii.gz'
my.MRIwrite(RSlin.detach().cpu().numpy(), aff, fake_filename)
del II, JJ, KK, II2, JJ2, KK2
del RSlin


##############################################
# OK we're ready for nonlinear registration! #
##############################################
print('Nonlinear registration of fake image with FireANTs')
# It's sad we need to write to disk... but there's no easy way of biulding a SimpleITK object
# from a generic affine vox2ras matrix :-S
ref_filename = output_dir + '/temp_reference.nii.gz'
my.MRIwrite((Iim*Mim).detach().cpu().numpy(), aff, ref_filename)

# FireANTs!
# Prepare reference image, with intensities in 1st channel, and a bunch of segmentations concatenated
image1 = Image.load_file(ref_filename, device=device_registration)
SimP = (Sim*Mim).permute([2, 1, 0])[None, None, ...]
# GM, WM (with basal forebrain), [CGM unless cerebrum mode], CA/PU/AC, TH, PA, [BG unless cerebrum mode], AM;
# Note that I had to insert hypothalamus before the amygdala (the special cases go at the end)
# In Feb'25, hypothal [13, 14], <amygdala goes here> fornix, optic chiasm, mammillary body, septal nucleus;
# I also killed WM because of ambighuities in subthalamic area (plus, its contour is modeled by cortical+subcortical GM anyway
if (side=='left') and ((mode=='invivo') or (mode=='exvivo')):
    groups = [[3,  17],#[2, 77,865],
              [8], [11, 12, 26], [10], [13], [0,  4,  5, 14, 15, 24, 41, 42, 43, 44, 46, 47, 49, 50, 51, 52, 53, 54, 58, 60, 162], [819], [18], [821], [85], [843], [869]]
if (side=='left') and ((mode=='hemi') or (mode=='cerebrum')):
    groups = [[3, 17], #[2, 77,865],
              [11, 12, 26], [10], [13], [819], [18], [821], [85], [843], [869] ]
if (side=='right') and ((mode=='invivo') or (mode=='exvivo')):
    groups = [[42, 53], # [41, 77,866],
              [47], [50, 51, 58], [49], [52], [0, 43, 44, 14, 15, 24,  2,  3,  4,  5,  7,  8, 10, 11, 12, 13, 17, 18, 26, 28, 161], [820], [54], [822], [85], [844], [870]]
if (side=='right') and ((mode=='hemi') or (mode=='cerebrum')):
    groups = [[42, 53], #[41, 77,866],
              [50, 51, 58], [49], [52], [820], [54], [822], [85], [844], [870] ]
for group in groups:
    M = torch.zeros(SimP.shape, device=device_registration, dtype=dtype)
    for lab in group:
        M[SimP == lab] = 1.0
    image1.array = torch.cat([image1.array, M], dim=1)
del SimP, M

# Same for the fake image (requires resampling atlas)
image2 = Image.load_file(fake_filename, device=device_registration)

# Now let's apply the linear transform to the frames / labels of the atlas for the label loss
II, JJ, KK = np.meshgrid(np.arange(Iim.shape[0]), np.arange(Iim.shape[1]), np.arange(Iim.shape[2]), indexing='ij')
II = torch.tensor(II, device=device, dtype=dtype)
JJ = torch.tensor(JJ, device=device, dtype=dtype)
KK = torch.tensor(KK, device=device, dtype=dtype)
# fakeVox <- fakeRAS <- imageVox
affine = np.linalg.inv(aff_A) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat
II2 = affine[0, 0] * II + affine[0, 1] * JJ + affine[0, 2] * KK + affine[0, 3]
JJ2 = affine[1, 0] * II + affine[1, 1] * JJ + affine[1, 2] * KK + affine[1, 3]
KK2 = affine[2, 0] * II + affine[2, 1] * JJ + affine[2, 2] * KK + affine[2, 3]
del II, JJ, KK
for f in range(len(A_reg)):
    pad = 1.0 if (np.any(label_sets_reg[f]==0)) else 0.0
    M = my.fast_3D_interp_torch(torch.tensor(A_reg[f], device=device, dtype=dtype),
                                II2 - atlas_bounds_segs[f,0],
                                JJ2 - atlas_bounds_segs[f,2],
                                KK2 - atlas_bounds_segs[f,4], 'linear', pad_value=pad)
    image2.array = torch.cat([image2.array, M.permute([2, 1, 0])[None, None, ...].to(device_registration)], dim=1)
del II2, JJ2, KK2, M, A_reg

#FireANTs options
# maxi = max(image1.array[0,0].max(), image2.array[0,0].max())
maxi = 110 # we use the median of the white matter
image1.array[0,0] /= maxi
image2.array[0,0] /= maxi
torch.set_grad_enabled(True)
batch1 = BatchedImages([image1])
batch2 = BatchedImages([image2])
scales = [4, 2, 1]; iterations = [300, 100, 30]
if np.mean(native_resolution)<0.7:
    scales = [6, 4, 2, 1]; iterations = [300, 100, 30, 30]
if np.mean(native_resolution)<0.5:
    scales = [8, 4, 2, 1]; iterations = [300, 100, 30, 30]
if np.mean(native_resolution)<0.3:
    scales = [12, 8, 4, 2, 1]; iterations = [300, 100, 30, 30, 30]
if np.mean(native_resolution) < 0.2:
    scales = [16, 8, 4, 2, 1]; iterations = [300, 100, 30, 30, 30]

reg = GreedyRegistration(scales=scales, iterations=iterations,
            fixed_images=batch1, moving_images=batch2,
            loss_type='custom',
            custom_loss=HybridDiceLabelDiffloss(kernel_size=args.cc_kernel_size, rel_weight_labeldiff=args.rel_weight_labeldiff),
            deformation_type='compositive',
            smooth_grad_sigma=args.smooth_grad_sigma, smooth_warp_sigma=args.smooth_warp_sigma,
            optimizer='adam', optimizer_lr=args.optimizer_lr,
            reduction='mean')  # crucial option, for proper display of loss...

start = time()
reg.optimize(save_transformed=False)
end = time()
print("Runtime nonlinear registration: ", end - start, "seconds")
moved = reg.evaluate(batch1, batch2)
reference_img = image1.itk_image
moved_image_np = (maxi * moved[0, 0]).detach().cpu().numpy()
moved_sitk_image = sitk.GetImageFromArray(moved_image_np)
moved_sitk_image.SetOrigin(reference_img.GetOrigin())
moved_sitk_image.SetSpacing(reference_img.GetSpacing())
moved_sitk_image.SetDirection(reference_img.GetDirection())
if args.save_atlas_nonlinear_reg:
    fake_filename_deformed = output_dir + '/atlas_nonlinear_reg.' + side + '.nii.gz'
    sitk.WriteImage(moved_sitk_image, fake_filename_deformed)
warped_coords = reg.get_warped_coordinates(batch1, batch2).to(device).detach()

# Save deformation field and Jacobian if needed!
coords = None
field_sitk_image = None
jacdet_sitk_image = None
if args.save_field:
    coords = (1 + warped_coords[0]) * (0.5 * (torch.tensor(Iim.shape, device=device) - 1))[None, None, None, ...]
    field_sitk_image = sitk.GetImageFromArray(coords.cpu().numpy())
    field_sitk_image.SetOrigin(reference_img.GetOrigin())
    field_sitk_image.SetSpacing(reference_img.GetSpacing())
    field_sitk_image.SetDirection(reference_img.GetDirection())
    field_filename = output_dir + '/nonlinear_field.' + side + '.nii.gz'
    sitk.WriteImage(field_sitk_image, field_filename)
    del field_sitk_image
if args.save_jacobian:
    if coords is None: # compute coords only if needed
        coords = (1 + warped_coords[0]) * (0.5 * (torch.tensor(Iim.shape, device=device) - 1))[None, None, None, ...]
    coords = coords.permute([2, 1, 0, 3])
    jacdet = my.jacobian_det_torch(coords)
    jacdet_sitk_image = sitk.GetImageFromArray(torch.log10(jacdet.permute([2,1,0]).abs().clip(min=1e-6)).cpu().numpy())
    jacdet_sitk_image.SetOrigin(reference_img.GetOrigin())
    jacdet_sitk_image.SetSpacing(reference_img.GetSpacing())
    jacdet_sitk_image.SetDirection(reference_img.GetDirection())
    jacdet_filename = output_dir + '/nonlinear_jac_logdet.' + side + '.nii.gz'
    sitk.WriteImage(jacdet_sitk_image, jacdet_filename)
    del jacdet, jacdet_sitk_image

del image1, image2, batch1, batch2, reg, moved, moved_sitk_image, coords
torch.set_grad_enabled(False)

########################################################
print('Computing initial values for means and variances')
mus_ini = []
vars_ini = []
mixture_weights = []

# First, the background
if SET_BG_TO_CSF:
    x = []
    if (side=='left') or (side=='left-c') or (side=='left-ccb'):
        for l in [4, 5]: # , 24]:
            x.append(Iim_before_masking[Sim==l])
    else:
        for l in [43, 44]: # , 24]:
            x.append(Iim_before_masking[Sim==l])
    mu_bg  = torch.median(torch.concatenate(x))
    if torch.isnan(mu_bg):
        mu_bg = torch.tensor(0, dtype=dtype, device=device)
else:
    mu_bg = torch.tensor(0, dtype=dtype, device=device)

# Now, the rest of classes
for t in range(len(number_of_gmm_components)-1):
    x = []
    for l in grouping_labels[t+1]:
        x.append(Iim_before_masking[Sim==l])
    x = torch.concatenate(x) if len(x)>0 else []

    if len(x) == 0: # if there's nothing in there (eg, cerebrum mode), grab from WM or GM as needed
        nc = number_of_gmm_components[t + 1]
        if (grouping_labels[t+1][0]==7) or (grouping_labels[t+1][0]==161) or (grouping_labels[t+1][0]==46) or (grouping_labels[t+1][0]==162):
            mus_ini.append(torch.mean(mus_ini[0]) * torch.ones(nc, dtype=dtype, device=device))
            vars_ini.append(torch.mean(vars_ini[0]) * torch.ones(nc, dtype=dtype, device=device))
        elif (grouping_labels[t+1][0] == 8) or (grouping_labels[t+1][0] == 47):
            mus_ini.append(torch.mean(mus_ini[1]) * torch.ones(nc, dtype=dtype, device=device))
            vars_ini.append(torch.mean(vars_ini[1]) * torch.ones(nc, dtype=dtype, device=device))
        else:
            raise Exception('Could not initialize label group ' + str(t+1) + ' (this should not happen unless SynthSeg did really poorly)')
        mixture_weights.append((1 / float(nc)) * torch.ones(nc, dtype=dtype, device=device))
    else:

        mu = torch.median(x)
        std = 1.4826 * torch.median(torch.abs(x - mu))
        var = std ** 2
        if number_of_gmm_components[t+1]==1:
            mus_ini.append(mu[None])
            vars_ini.append(var[None])
            mixture_weights.append(torch.ones(1,dtype=dtype,device=device))
        else:
            # Estimate GMM with shared variance (avoids a component with tiny variance)
            nc = number_of_gmm_components[t+1]
            nx = len(x)
            gmm_mus = torch.linspace(mu - 0.5 * std, mu + 0.5 * std, nc, dtype=dtype, device=device)
            gmm_var= var * torch.ones(1, dtype=dtype, device=device)
            gmm_ws = (1 / float(nc)) * torch.ones(nc, dtype=dtype, device=device)
            W = torch.zeros([nx, nc], dtype=dtype, device=device)
            for its in range(200):
                # E step
                for c in range(nc):
                    W[:, c] = gmm_ws[c] / torch.sqrt(2.0 * torch.pi * torch.sqrt(gmm_var)) * torch.exp(-0.5 * (x - gmm_mus[c])**2 / gmm_var)
                normalizer = torch.sum(W + 1e-9, axis=1)
                # print(-torch.mean(torch.log(normalizer)))
                W /= normalizer[:, None]
                # M step
                denominators = torch.sum(W, axis=0)
                gmm_ws = denominators / torch.sum(denominators)
                gmm_var = 0
                for c in range(nc):
                    gmm_mus[c] = torch.sum(W[:, c] * x) / denominators[c]
                    aux = x - gmm_mus[c]
                    gmm_var += torch.sum(W[:, c] * aux * aux)
                gmm_var /= torch.sum(denominators)

            mus_ini.append(gmm_mus)
            vars_ini.append(gmm_var * torch.ones(nc, dtype=dtype, device=device))
            mixture_weights.append(gmm_ws)

mus_ini = torch.concatenate(mus_ini)
vars_ini = torch.concatenate(vars_ini)
mixture_weights = torch.concatenate(mixture_weights)
# Replace -1 variance with minimum otherwise
vars_ini[vars_ini < 0] = torch.min(vars_ini[vars_ini>0])

########################################################

# EM for GMM parameters
print('Estimating GMM parameters')

print('  Resizing image, mask, and coordinates')
I_r, aff_r = my.torch_resize(Iim_before_masking, aff, resolution, device, dtype=dtype)
M_r, _ = my.torch_resize(Mim, aff, resolution, device, dtype=dtype)
# smoothen M_r a bit
kernel = torch.zeros([3,3,3], device=device, dtype=dtype)
kernel[1,1,:] = 0.1; kernel[1,:,1] = 0.1; kernel[:,1,1] = 0.1; kernel[1,1,1] = 0.4
for _ in range(args.smoothing_steps_HRmask):
    M_r = torch.conv3d(M_r[None, None, ...], kernel[None, None, ...], bias=None, stride=1, padding=[1, 1, 1]).squeeze()
I_r[M_r<0.5] = mu_bg
Iim_shape = Iim.shape
del Iim, Iim_before_masking, Mim, Sim

print('  Deforming atlas')
# OK so we need to create a field of coordinates and then resize it.
# The tricky bit is concatenating the nonlinear transform (in [0,1]) with the linear
wc_r, _ = my.torch_resize(warped_coords[0].permute([2,1,0,3]), aff, resolution, device, dtype=dtype)
del warped_coords
# we'll worry about skip later
# first bit: from [-1,1] of linearly registered, to absolute coordinates
T1 = torch.eye(4, device=device, dtype=dtype)
for k in range(3):
    T1[k, k] = 0.5 * (Iim_shape[k] - 1)
    T1[k, -1] = T1[k, k]
# second bit: vox2vox transform to fake  space  fake_vox <- image_vox
T2 = torch.tensor(np.linalg.inv(aff_fake) @ Mmni @ M_input_vox_to_mni_ras @ shift_mat, device=device, dtype=dtype) # TODO: is this correct?
# third bit: from vox to [-1, 1] coordinates that can be used with atlas at any resolution
T3 = torch.eye(4, device=device, dtype=dtype)
for k in range(3):
    T3[k, k] = 2 / (Ifake.shape[k] - 1)
    T3[k, -1] = -1
T = T3 @ T2 @ T1
I = T[0, 0] * wc_r[..., 0] + T[0, 1] * wc_r[..., 1] + T[0, 2] * wc_r[..., 2] + T[0, 3]
J = T[1, 0] * wc_r[..., 0] + T[1, 1] * wc_r[..., 1] + T[1, 2] * wc_r[..., 2] + T[1, 3]
K = T[2, 0] * wc_r[..., 0] + T[2, 1] * wc_r[..., 1] + T[2, 2] * wc_r[..., 2] + T[2, 3]
del wc_r

# We can now resample

if False: # a bit faster, but much more GPU memory hungry. Moves A to the GPU
    A = torch.tensor(A, device=device, dtype=dtype)
    priors = grid_sample(A.permute([3,0,1,2])[None, ...],
                         torch.stack([K[::skip,::skip,::skip],
                                      J[::skip,::skip,::skip],
                                      I[::skip,::skip,::skip]], axis=-1)[None,...], align_corners=True)
else: # slower, but more memory frugal. A remains on the CPU!
    priors = []
    locs = torch.stack([K[::skip, ::skip, ::skip], J[::skip, ::skip, ::skip], I[::skip, ::skip, ::skip]], axis=-1)[None, ...]
    aux = torch.zeros(*atlas_size, device=device, dtype=dtype)
    for f in range(n_tissues):
        aux[:] = 0
        aux[atlas_bounds_tissues[f,0]:atlas_bounds_tissues[f,1],
            atlas_bounds_tissues[f,2]:atlas_bounds_tissues[f,3],
            atlas_bounds_tissues[f,4]:atlas_bounds_tissues[f,5]] = torch.tensor(A[f], device=device, dtype=dtype)
        priors.append(grid_sample(aux[None, None, ...], locs, align_corners=True))
    priors = torch.concat(priors, dim=1)
    del locs, aux

# Reshape and deal with voxels outside the FOV
priors = torch.permute(priors[0], [1, 2, 3, 0])
missing_mass = 1 - torch.sum(priors, axis=-1)
priors[..., 0] += missing_mass

####
print('  EM for parameter estimation')
# data
means = torch.tensor([mu_bg, *mus_ini], device=device, dtype=dtype)
var_bg = torch.min(vars_ini)
variances = torch.tensor([var_bg, *vars_ini], device=device, dtype=dtype)
weights = torch.tensor([1.0, *mixture_weights], dtype=dtype, device=device)
W = torch.zeros([*priors.shape[:-1], number_of_gmm_components.sum()], dtype=dtype, device=device)
x = I_r[::skip, ::skip, ::skip]

# We now put a Scaled inverse chi-squared prior on the variances to prevent them from going to zero
prior_count = 100 / (resolution ** 3) / (skip ** 3)
prior_variance = var_bg
loglhood_old = -10000
eps = 1e-7
for em_it in range(100):
    # E step
    for c in range(n_tissues):
        prior = priors[:, :, :, c]
        num_components = number_of_gmm_components[c]
        for g in range(num_components):
            gaussian_number = sum(number_of_gmm_components[:c]) + g
            d = x - means[gaussian_number]
            W[:, :, :, gaussian_number] = weights[gaussian_number] * prior * torch.exp(
                -d * d / (2 * variances[gaussian_number])) / torch.sqrt(
                2.0 * torch.pi * variances[gaussian_number])

    normalizer = eps + torch.sum(W, dim=-1, keepdim=True)
    loglhood = torch.mean(torch.log(normalizer)).detach().cpu().numpy()
    W = W / normalizer

    # M step
    prior_loglhood = torch.zeros(1,dtype=dtype,device=device)
    for c in range(number_of_gmm_components.sum()):
        # crucially, we skip the background when we update the parameters (but we still add it to the cost)
        if c > 0:
            norm = eps + torch.sum(W[:, :, :, c])
            means[c] = torch.sum(x * W[:, :, :, c]) / norm
            d = x - means[c]
            variances[c] = (torch.sum(d * d * W[:, :, :, c]) + prior_count * prior_variance) / ( norm + prior_count + 2)
        v = variances[c]
        prior_loglhood = prior_loglhood - ((1 + 0.5 * prior_count) * torch.log(v) + 0.5 * prior_count * prior_variance / v) / torch.numel(normalizer)
    loglhood = loglhood + prior_loglhood.detach().cpu().numpy()

    mixture_weights = torch.sum(W[:, :, :, 1:].reshape([np.prod(priors.shape[:-1]), number_of_gmm_components.sum() - 1]) + eps, axis=0)
    for c in range(n_tissues - 1):
        # mixture weights are normalized (those belonging to one mixture sum to one)
        num_components = number_of_gmm_components[c + 1]
        gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[1:c + 1]) + \
                                        np.array(range(num_components)), device=device, dtype=dtype).long()

        mixture_weights[gaussian_numbers] /= torch.sum(mixture_weights[gaussian_numbers])

    weights[1:] = mixture_weights

    if (torch.sum(torch.isnan(means)) > 0) or (torch.sum(torch.isnan(variances)) > 0):
        print('nan in Gaussian parameters...')
        import pdb;

        pdb.set_trace()

    print('         Step %d of EM, -loglhood = %.6f' % (em_it + 1, -loglhood), flush=True, end='\r')
    if (loglhood - loglhood_old) < TOL:
        print(' ')
        print('         Decrease in loss below tolerance limit')
        break
    else:
        loglhood_old = loglhood
print(' ')


###########

print('Computing Gaussians at full resolution (we will reuse over and over)')
GAUSSIAN_LHOODS = torch.zeros([*I_r.shape, sum(number_of_gmm_components)], dtype=dtype, device=device)
for c in range(sum(number_of_gmm_components)):
    # The 1e-9 ensures no zeros were prior is not zero
    GAUSSIAN_LHOODS[..., c] = 1.0 / torch.sqrt(2 * math.pi * variances[c]) * torch.exp(
    -0.5 * torch.pow(I_r - means[c], 2.0) / variances[c]) + eps

print('Computing normalizers (faster to do now with clustered priors)')
# We deform one class at the time; slower, but less memory
normalizers = torch.zeros(GAUSSIAN_LHOODS.shape[:-1], dtype=dtype, device=device)
# normalizers[h] = eps
gaussian_number = 0
locs = torch.stack([K, J, I], axis=-1)[None, ...]
aux = torch.zeros(*atlas_size, device=device, dtype=dtype)
for c in range(len(A)):
    aux[:] = 0
    aux[atlas_bounds_tissues[c, 0]:atlas_bounds_tissues[c, 1],
        atlas_bounds_tissues[c, 2]:atlas_bounds_tissues[c, 3],
        atlas_bounds_tissues[c, 4]:atlas_bounds_tissues[c, 5]] = torch.tensor(A[c], device=device, dtype=dtype)
    prior = grid_sample(aux[None, None, ...], locs, align_corners=True)[0,0,...]

    if c==0: # background
        prior[(I < (-1)) | (I > 1) | (J < (-1)) | (J > 1) | (K < (-1)) | (K > 1)] = 1.0
    lhood = torch.zeros_like(prior)
    for g in range(number_of_gmm_components[c]):
        lhood += (weights[gaussian_number] * GAUSSIAN_LHOODS[..., gaussian_number])
        gaussian_number += 1
    normalizers += (prior * lhood)
del A, locs, aux

########

print('Deforming one label at a time')
names, colors = my.read_LUT(LUT_file)
seg = torch.zeros(normalizers.shape, dtype=torch.int, device=device)
seg_rgb = torch.zeros([*normalizers.shape, 3], dtype=dtype, device=device)
max_p = torch.zeros(normalizers.shape, dtype=dtype, device=device)
vols = torch.zeros(n_labels, device=device, dtype=dtype)

# TODO: choose good number of workers/prefetch factor
for n, (prior_indices, prior_values) in enumerate(label_loader):
    print('  Deforming label ' + str(n + 1) + ' of ' + str(n_labels), end='\r')

    if prior_indices.numel() == 0:
        continue
    prior_indices = torch.as_tensor(prior_indices, device=device, dtype=torch.long).squeeze()
    prior_values = torch.as_tensor(prior_values, device=device, dtype=dtype).squeeze()

    if n == 0:
        # background
        prior = torch.sparse_coo_tensor(prior_indices[None], prior_values,
                                        [torch.Size(atlas_size).numel()]).to_dense()
        del prior_indices, prior_values
        prior = prior.reshape(torch.Size(atlas_size))

    else:

        # find bounding box of label in atlas space
        prior_indices = my.ind2sub(prior_indices, atlas_size)
        min_x, max_x = prior_indices[0].min().item(), prior_indices[0].max().item()
        min_y, max_y = prior_indices[1].min().item(), prior_indices[1].max().item()
        min_z, max_z = prior_indices[2].min().item(), prior_indices[2].max().item()
        crop_atlas_size = [max_x - min_x + 1, max_y - min_y + 1, max_z - min_z + 1]
        prior_indices[0] -= min_x
        prior_indices[1] -= min_y
        prior_indices[2] -= min_z
        prior = torch.sparse_coo_tensor(prior_indices, prior_values, crop_atlas_size).to_dense()
        del prior_indices, prior_values

    skip_this_label = False
    if n==0:
        Irescaled = I
        Jrescaled = J
        Krescaled = K
        lr_crop = (slice(None),) * 3
    else:
        # find bounding box of label in MRI space
        Irescaled = I * ((atlas_size[0] - 1) / (max_x - min_x)) + ( (atlas_size[0] - 1 - min_x - max_x) / (max_x - min_x) )
        Jrescaled = J * ((atlas_size[1] - 1) / (max_y - min_y)) + ( (atlas_size[1] - 1 - min_y - max_y) / (max_y - min_y) )
        Krescaled = K * ((atlas_size[2] - 1) / (max_z - min_z)) + ( (atlas_size[2] - 1 - min_z - max_z) / (max_z - min_z) )
        mask = (Irescaled >= (-1))
        mask &= (Irescaled <= 1)
        mask &= (Jrescaled >= (-1))
        mask &= (Jrescaled <= 1)
        mask &= (Krescaled >= (-1))
        mask &= (Krescaled <= 1)
        if mask.any()==False:
            skip_this_label = True
        else:
            nx, ny, nz = mask.shape
            tmp = mask.reshape([nx, -1]).any(-1).nonzero()
            lr_min_x, lr_max_x = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.movedim(0, -1).reshape([ny, -1]).any(-1).nonzero()
            lr_min_y, lr_max_y = tmp.min().item(), tmp.max().item() + 1
            tmp = mask.reshape([-1, nz]).any(0).nonzero()
            lr_min_z, lr_max_z = tmp.min().item(), tmp.max().item() + 1
            del tmp, mask
            lr_crop = (slice(lr_min_x, lr_max_x), slice(lr_min_y, lr_max_y), slice(lr_min_z, lr_max_z))

    if skip_this_label==False:
        prior_resampled = grid_sample(prior[None, None, ...], torch.stack([Krescaled[lr_crop], Jrescaled[lr_crop], Irescaled[lr_crop]], axis=-1)[None, ...], align_corners=True)[0, 0, ...]
        if n==0: # background
            prior_resampled[(Krescaled[lr_crop]<(-1)) | (Krescaled[lr_crop]>1) | (Jrescaled[lr_crop]<(-1)) | (Jrescaled[lr_crop]>1)  | (Irescaled[lr_crop]<(-1)) | (Irescaled[lr_crop]>1)  ] = 1.0
        del Irescaled; del Jrescaled; del Krescaled
        num_components = number_of_gmm_components[tissue_index[n]]
        gaussian_numbers = torch.tensor(np.sum(number_of_gmm_components[:tissue_index[n]]) + \
                                np.array(range(num_components)), device=device, dtype=dtype).long()
        lhood = torch.sum(GAUSSIAN_LHOODS[:, :, :, gaussian_numbers] * weights[None, None, None, gaussian_numbers], 3)
        post = torch.clone(prior_resampled)
        post *= lhood[lr_crop]
        post /= normalizers[lr_crop]
        if n==0:
            post[torch.isnan(post)] = 1.0
        else:
            post[torch.isnan(post)] = 0.0
            if (mode=='hemi') or (mode=='cerebrum'): # cerebrum only: kill labels outside mask
                post *= M_r[lr_crop]
        # Also, for cerebrum only, kill cerebellar labels (not brainstem as we don't know where to stop)
        if ((mode == 'hemi') or (mode == 'cerebrum ')) and (np.sum(np.array([595, 597, 715, 721, 751, 752, 846])==label_list[n])):
            post[:] = 0
        del prior_resampled
        vols[n] = torch.sum(post) * (resolution ** 3)
        mask = (post > max_p[lr_crop])
        max_p[lr_crop][mask] = post[mask]
        lab = int(label_list[n])
        seg[lr_crop].masked_fill_(mask, lab)
        del mask
        for c in range(3):
            seg_rgb[(*lr_crop, c)].add_(post, alpha=colors[lab][c])
        del post
print('\n\n')

########################################################

print('Writing results to disk')
my.MRIwrite(seg.detach().cpu().numpy().astype(np.uint16), aff_r, output_dir + '/seg.' + side + '.nii.gz')
if args.write_rgb:
    my.MRIwrite(seg_rgb.detach().cpu().numpy().astype(np.uint8), aff_r, output_dir + '/seg.' + side + '.rgb.nii.gz')
vols = vols.detach().cpu().numpy()
with open(output_dir + '/vols.' + side + '.csv', 'w') as csvfile:
    writer = csv.writer(csvfile)
    aux = label_list[1:]
    row = []
    for l in aux:
        row.append(names[int(l)])
    writer.writerow(row)
    row = []
    for j in range(1, len(vols)):
        row.append(str(vols[j]))
    writer.writerow(row)
# Copy LUT file, for convenience
a = os.system('cp ' + LUT_file + ' ' + output_dir +  '/lut.txt >/dev/null')
if a==0:
    LUT_file = output_dir +  '/lut.txt'

# Clean up
os.system('rm -rf '  + output_dir +  '/temp*.nii.gz >/dev/null')

# Print commands to visualize output
print('You can try these commands to visualize outputs:')
cmd = '  freeview -v ' + input_volume + ' -v ' + output_dir + '/SuperSynth/segmentation.mgz:colormap=lut '
if args.write_bias_corrected and (skip_bf==False):
    cmd = cmd + ' -v ' + output_dir + '/bias.corrected.' + side + '.mgz'
cmd = cmd + ' -v ' + output_dir + '/seg.' + side + '.nii.gz:colormap=lut:lut=' + LUT_file
if args.write_rgb:
    cmd = cmd + ' -v ' + output_dir + '/seg.' + side + '.rgb.nii.gz:rgb=true'
if args.save_atlas_nonlinear_reg:
    cmd = cmd + ' ' + fake_filename_deformed
if args.save_jacobian:
    cmd = cmd + ' ' + jacdet_filename
print(cmd)
print('  oocalc ' + output_dir + '/vols.' + side + '.csv')
print(' ')
print('All done!')


##################

now2 = datetime.now()

current_time = now2.strftime("%H:%M:%S")
print("Current Time =", current_time)

runtime = now2 - now

print("Running Time =", runtime)

