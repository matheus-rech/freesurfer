This is an implementation of an approximate Bayesian segmentation method that relies 
on the histological  atlas "NextBrain", presented in the article:
"A probabilistic histological atlas of the human brain for MRI segmentation", 
Casamitjana et al., Nature, 2025 (https://www.nature.com/articles/s41586-025-09708-2)

The actual approximate inference method is presented in the article:
"Fast segmentation with the NextBrain histological atlas", 
Puonti et al. (under revision).

The code also relies on a foundation model based on our article:
"A Modality-agnostic Multi-task Foundation Model for Human Brain Imaging",
Liu et al. (under revision),

As well as on FireANTS, an efficient registration algorith from UPenn:
"FireANTs: Adaptive Riemannian Optimization for Multi-Scale Diffeomorphic Registration",
Jena el al. (under revision).

## Prerequisites

The first time you run the method, it will prompt you to download the atlas.
If you have never used the command mri_super_synth before, it will also prompt
you to download a model file.

## Basic Usage

mri_histo_atlas_segment_fireants --i INPUT_SCAN --o OUTPUT_DIRECTORY --device [cpu/cuda]  --side [left/right] --mode [invivo/cerebrum/hemi/exvivo]

## Outputs

The output directory will contain the following files:

- seg.[left/right].nii.gz: segmentation of left/right hemisphere
- lut.txt: the lookup table to visualize seg.[left/right].nii.gz, for convenience
- vols.[left/right].csv: files with volumes of the brain regions segmented by the atlas, in CSV format.
- SuperSynth: directory segmentation of the scan at the whole structure level (from the foundation model)

and additional files, e.g., nonlinearly registered atlas, depending on the specific flags.


## Advanced options 
The code also accepts the following optional flags:

- `--bf_mode`: Decides the bias field basis function model. Options: dct (default), polynomial, hybrid.
- `--write_rgb`: Save an RGB image based on the posterior probabilites to disk.
- `--write_bias_corrected`: Save the bias corrected input image to disk.
- `--device_registration`: Define a different device for the registration. Can be used to save GPU memory when working with an GPU with limited memory. Options: cpu, cuda. Default is the same as --device.
- `--threads`: Control the number of cpu thread used to run the algorithm. Default value is -1, which uses all available threads.
- `--skip`: An integer skipping (downsampling) factor for estimating the model parameters. More skipping saves memory, but sacrifices accuracy. Default: 1 (no skipping).
- `--resolution`: The resolution of the output segmentation. By default 0.4mm, which is higher than the typical input scan, to reduce aliasing. 
- `--smoothing_steps_HRmask`: Number of smoothing steps used when upsampling the 1mm brain mask from BrainFM. More smoothing makes the outer border less jagged, but too much smoothing reduces accuracy. Default: 3.
- `--skip_bf`: Skip the bias field correction. Can be used to save memory if the input scan is already bias corrected or does not have a bias field (non MRI modality).
- `--smooth_grad_sigma`: Gradient field smoothing parameter for the nonlinear FireAnts registration. Default: 1.0.
- `--smooth_warp_sigma`: Warp field smoothing parameter for the nonlinear FireAnts registration. Default: 0.25.
- `--optimizer_lr`: Learning rate for the nonlinear FireAnts registration optimizer. Default: 0.5.
- `--cc_kernel_size`: Size of the window for calculating the cross-correlation registration metric. Default: 7.
- `--rel_weight_labeldiff`: Relative weight for the Dice loss metric in the nonlinear registration. Default: 2.5.
- `--save_atlas_nonlinear_reg`: Save the nonlinearly registered atlas. Default: false.
- `--save_field`: Save the nonlinear deformation field. Default: false.
- `--save_jacobian`: Save the Jacobian determinant (log10) of the deformation field. Default: false.

Some notes:

- If you are running out of memory, using skip=2 can help without sacrificing much accuracy.
- The defaults smooth_grad_sigma=1,smooth_warp_sigma=0.25 are pretty liberal and can cope with massive deformation, e.g., as in the Hip-CT images shown in the paper "Fast segmentation with the NextBrain". If you are working with a population without very strong atropy or deformation, you can multiply those values by 2 or 3 and get more regular atlas deformation fields (you can explore the deformation with the --save_jacobian option.

Also: you can flexibly change the groupings of the modeled structures using the .yaml files under the /data_simplified folder. 
The structure groupings for the Gaussian Mixture modeled are controlled by two files: `gmm_components_fireants.yaml` and `combined_atlas_labels_fireants.yaml`.
Let's say, as an example, that you wanted to add the internal segment of globus pallidus (label 206) as its own structure.
To model it separately, you would first create a new class, called e.g., Internal Segment Pallidum, in the `combined_atlas_labels_fireants.yaml` file, and list label 206 under that structure (while removing it from the pallidum class).
Next, you would add the class, with exactly the same name, to the `gmm_components_fireants.yaml` file and decide how many Gaussian distributions should be used to model its intensities.
To make the non-linear registration aware of the contrast, you would add the structure, again with exactly the same name, to the file called `recipe_intensities_cheating_image_fireants.yaml`, 
and decide how its intensity should be generated from the seven structures than can be always reliably segmentation using BrainFM (see the file for examples).





