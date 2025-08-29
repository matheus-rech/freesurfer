This is an implementation of a Bayesian segmentation method that relies on the histological atlas presented in the article:
"Next-Generation histological atlas and segmentation tool for echo "high-resolution in vivo human neuroimaging", by 
Casamitjana et al. 
-- preprint available at https://www.biorxiv.org/content/10.1101/2024.02.05.579016v1 

And the method presented in the article:
"Fast segmentation with the NextBrain histological atlas", by
Puonti et al. (under revision)

The code also relies on:
"A Modality-agnostic Multi-task Foundation Model for Human Brain Imaging"
Liu et al. (under revision)

## Prerequisites:

The first time you run the method, it will prompt you to download the atlas files, which are not distributed with the code.
If you use the FireANTs version (highly recommended), it will also prompt you to download a machine learning model file.s


## Usage for 'full' Bayesian version (slow, not recommended):

To run the code, please use the script segment.sh as follows:

mri_histo_atlas_segment INPUT_SCAN OUTPUT_DIRECTORY ATLAS_MODE GPU THREADS [BF_MODE] [GMM_MODE]

- INPUT SCAN: scan to process, in nii(.gz) or mgz format
- OUTPUT_DIRECTORY: directory where segmentations, volume files, etc will be written (more on this below).
- ATLAS_MODE: must be full (all 333 labels) or simplified (simpler brainstem protocol; recommended)
- GPU: set to 1 to use the GPU (*highly* recommended but requires a 24GB GPU!)
- THREADS: number of CPU threads to use (use -1 for all available threads)
- BF_MODE (optional): bias field mode: dct (default), polynomial, or hybrid
- GMM_MODE (optional): gaussian mixture model (GMM) model must be 1mm unless you define your own (see documentation)

Note that the first time that you run the code, you may be prompted you to download the atlas separately.

Also, Using a GPU (minimum memory: 24GB) is highly recommended. On the GPU, the code runs in about an hour (30 mins/hemisphere).
On the CPU, the running time depends heavily on the number of threads, but it can easily take over 10 hours if you do not
use many (>10) threads! Even if you use the GPU, we recommend using a bunch of CPU threads (e.g., 8) if possible, so the CPU 
parts of the algorithm run faster.

The default bias field mode (dct) uses a set of discrete cosine transform basis functions to model the bias field. The
polynomial mode uses a set of low-order 3D polynomials. The hybrid mode uses a combination of dct and polynomials.

The GMM model is crucial as it determines how different brain regions are grouped into tissue types for the purpose of 
image intensity modeling. This is specified though a set of files that should be found under data:

- data_[full/simplified]/gmm_components_[GMM_MODE].yaml: defines tissue classes and specificies the number of components of the corresponding GMM
- data_[full/simplified]/combined_aseg_labels_[GMM_MODE].yaml: defines the labels that belong to each tissue class
- data_[full/simplified]/combined_atlas_labels_[GMM_MODE].yaml: defines FreeSurfer ("aseg") labels that are used to initialize the parameters of each class.

We distribute a GMM_MODE named "1mm" that we have used in our experiments, and which is the default mode of the code. If you 
want to use your own model, you will need to create another triplet of files of your own (use the 1mm version as template).


## Output:

The output directory will contain the following files:

- bf_corrected.mgz: bias field corrected version of the input scan
- SynthSeg.mgz: SynthSeg segmentation of the scan at the whole structure level
- MNI_registration.mgz: deformation file with registration to MNI atlas (which can be found under data/mni.nii.gz)
- seg_[left/right].mgz: segmentation files (one per hemisphere).
- vols_[left/right].csv: files with volumes of the brain regions segmented by the atlas, in CSV format.
- lookup_table.txt: the lookup table to visualize seg_[left/right].mgz, for convenience
- done: this is an empty file that gets written upon successful completion of the pipeline.

You can visualize the output by CDing into the results directory and running the command:

freeview -v bf_corrected.mgz -v seg_left.mgz:colormap=lut:lut=lookup_table.txt -v seg_right.mgz:colormap=lut:lut=lookup_table.txt
 
## Alternative 'fast' versions:

We also distribute two faster version, where the atlas deformation is pre-computed
and then kept constant during the optimization, such that we only need to run the EM algorithm once for
the Gaussian parameters and that is it.

There are two sub-versions. One uses FireANTs (Jena et al., https://arxiv.org/abs/2404.01249) to register to the target
scan a fake cartoon derived from the atlas. This version actively tries to fit smaller structures and
subregions. The command line is:

mri_histo_atlas_segment_fireants INPUT_SCAN OUTPUT_DIRECTORY GPU THREADS [BF_MODE]

The other version uses SynthMorph (Hoffmann et al., Imaging Neuroscience, 2024), a neural network for 
image registration. This version is faster and does OK for 1mm isotropic scans. However, SynthMorph 
relies heavily on fitting the boundaries of whole structures and does not the map smaller regions 
as well (e.g., thalamic nuclei). Therefore, we do not recommend it for images with resolution better 
than 1mm isotropic (e.g., ex vivo scans).

mri_histo_atlas_segment_fast INPUT_SCAN OUTPUT_DIRECTORY GPU THREADS [BF_MODE]

For both commands, the options are similar to mri_histo_atlas_segment, but the atlas and gmm modes are always 'simplified' and '1mm', respectively. 
The output files in the output directory follow the same convention.

These faster versions are particularly useful if you are running the code on the CPU rather than CPU. 
On a semi-modern desktop, the run time should be less than an hour; note that mri_histo_atlas_segment_fireants
runs once per hemisphere, whhile mri_histo_atlas_segment_fast segments both hemispheres in a single run.




