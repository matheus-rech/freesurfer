This is an implementation of a Bayesian segmentation method that relies on the histological atlas presented in the article:
"A probabilistic histological atlas of the human brain for MRI segmentation", by 
Casamitjana et al. 
-- preprint available at https://www.biorxiv.org/content/10.1101/2024.02.05.579016v1 

And the method presented in the article:
"Fast segmentation with the NextBrain histological atlas", by
Puonti et al. (under revision).

The code also relies on a foundation model based on our article:
"A Modality-agnostic Multi-task Foundation Model for Human Brain Imaging",
Liu et al. (under revision),

As well as on FireANTS, an efficient registration algorith from UPenn:
"FireANTs: Adaptive Riemannian Optimization for Multi-Scale Diffeomorphic Registration",
Jena el al. (under revision).

## Prerequisites

The first time you run the method, it will prompt you to download the atlas and foundation model files, which are not
distributed with the code.

## Usage

mri_histo_atlas_segment_fireants INPUT_SCAN OUTPUT_DIRECTORY GPU THREADS MODE SIDE

- INPUT SCAN: scan to process, in nii(.gz) or mgz format
- OUTPUT_DIRECTORY: directory where segmentations, volume files, etc will be written (more on this below).
- GPU: set to 1 to use the GPU (*highly* recommended but requires a ~40GB GPU!)
- THREADS: number of CPU threads to use (use -1 for all available threads)
- MODE: type of scan: invivo, exvivo, cerebrum (ex vivo without brainstem or cerebellum), 
hemi (ex vivo with single cerebral hemisphere)
- SIDE: left or right


## Outputs

The output directory will contain the following files:

- seg.[left/right].nii.gz: segmentation of left/right hemisphere
- lookup_table.txt: the lookup table to visualize seg.[left/right].nii.gz, for convenience
- vols.[left/right].csv: files with volumes of the brain regions segmented by the atlas, in CSV format.
- supersynth.nii.gz: segmentation of the scan at the whole structure level (from the foundation model)
- supersynth.vols.csv: volumes estimated by foundation model


## Legacy versions

#### Synthmorph version

There is a version of the code that estimates the registration with SynthMorph (Hoffmann et al., Imaging Neuroscience, 2024),
- a neural network for image registration - rather than FireANTs. This version is pretty good for 1mm isotropic scans. However, 
SynthMorph  relies heavily on fitting the boundaries of whole structures and does not the map smaller regions 
(e.g., thalamic nuclei) as well as the FireANTs version. Therefore, we do not recommend it for images with resolution better 
than 1mm isotropic (e.g., ex vivo scans). Also, it does not support ex vivo scans, cerebra, or single hemispheres. 
the command is:

mri_histo_atlas_segment_synthmorph INPUT_SCAN OUTPUT_DIRECTORY GPU THREADS

The input arguments and structure of the output directory are very similar to those of mri_histo_atlas_segment_fireants.

##$# 'Full' Bayesian version (slow, not recommended):

We also distribute an implementation of the alrogithm described in the original paper, for reproducibility purposes,
but we do not recommend its use since it is really slow. Also, it does not support ex vivo scans, cerebra, or single hemispheres.
The command is:

mri_histo_atlas_segment_fullbayesian INPUT_SCAN OUTPUT_DIRECTORY ATLAS_MODE GPU THREADS [BF_MODE] [GMM_MODE]

- INPUT SCAN: scan to process, in nii(.gz) or mgz format
- OUTPUT_DIRECTORY: directory where segmentations, volume files, etc will be written (more on this below).
- ATLAS_MODE: must be full (all 333 labels) or simplified (simpler brainstem protocol; recommended)
- GPU: set to 1 to use the GPU (*highly* recommended but requires a 24GB GPU!)
- THREADS: number of CPU threads to use (use -1 for all available threads)
- BF_MODE (optional): bias field mode: dct (default), polynomial, or hybrid
- GMM_MODE (optional): gaussian mixture model (GMM) model must be 1mm unless you define your own (see documentation)

The GMM model is crucial as it determines how different brain regions are grouped into tissue types for the purpose of 
image intensity modeling. This is specified though a set of files that should be found under data:

- data_[full/simplified]/gmm_components_[GMM_MODE].yaml: defines tissue classes and specificies the number of components of the corresponding GMM
- data_[full/simplified]/combined_aseg_labels_[GMM_MODE].yaml: defines the labels that belong to each tissue class
- data_[full/simplified]/combined_atlas_labels_[GMM_MODE].yaml: defines FreeSurfer ("aseg") labels that are used to initialize the parameters of each class.

We distribute a GMM_MODE named "1mm" that we have used in our experiments, and which is the default mode of the code. If you 
want to use your own model, you will need to create another triplet of files of your own (use the 1mm version as template).
