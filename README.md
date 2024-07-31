# umvsasl
Written by David Frey for the fMRI lab at the University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our umvsasl project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for ASL (Arterial Spin Labeling) MRI experiments.

## Table of Contents
1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
2. [Usage](#usage)
    - [Pulse Sequence](#pulse-sequence)
        - [Basic Sequence Control](#basic-sequence-control)
            - [Looping Structure](#looping-structure)
            - [K-space Sampling Scheme](#kspace-sampling-scheme)
        - [Controlling the Readout](#controlling-the-readout)
            - [Fast Spin Echo (FSE) Mode](#fse-mode)
            - [Spoiled Gradient Echo (SPGR) Mode](#spgr-mode)
            - [Balanced Steady State Free Precession (bSSFP) Mode](#bssfp-mode)
            - [EPIC Readout Structure](#epic-readout-structure)
        - [Controlling the Preparation Pulses](#prep-pulses)
            - [Fat Suppression](#fat-suppression)
            - [Custom Prep Pulses](#custom-prep-pulses)
            - [PCASL Pulse](#pcasl-pulse)
            - [Pre-saturation Pulse](#pre-saturation-pulse)
3. [License](#license)
4. [Contact](#contact)

## Introduction

## Usage

### Pulse Sequence

#### Basic Sequence Control

##### Looping Structure
The main sequence follows a looping structure over three readout dimensions: frames, shots (through-plane interleaves), and arms (in-plane interleaves), ordered from fastest to slowest, respectively. Each frame contains all preparation pulses, delay times, and the readout echo train. [The entire kspace volume is encoded](#kspace-encoding-scheme) between the echo train, shots, and arms dimensions.

 <img width="500" src="https://github.com/user-attachments/assets/3e301dd8-2dc4-4241-bfb1-1efa7e5062d6">

| Parameter | Controlled by |
| - | - |
| Number of frames (M) | cv: `nframes` |
| Number of shots (N) | cv: `opnshots` |
| Number of arms (P) | cv: `narms` |
| Shot interval (Long TR) (us) | cv: `optr` |
| Sequence duration (us) | cv: `pitscan` (do not fix this variable - view only) |

##### Kspace sampling scheme
This sequence uses 3D stack of spirals and spiral projection imaging (SPI) to encode kspace.
An initial 2D [variable-density spiral](https://www.researchgate.net/publication/5925197_Fast_3D_imaging_using_variable-density_spiral_trajectories_with_applications_to_limb_perfusion) is generated and then transformed in 3D using rotation matrices.
Each readout dimension also transforms the kspace trajectory such that each frame contains all information for a single volume in kspace.
Each arm interleaves the initial 2D spiral in-plane by uniform rotations as shown on the left of the figure below.
Each shot interleaves the projected spiral volume acquired in a single echo echo train. For stack of spirals, this looks like through-plane interleaving along z as shown in the top right of the figure below. For spiral projection imaging, this will interleave the 3D rotations as shown in the bottom right of the figure below.

<img width="600" src="https://github.com/user-attachments/assets/7e6ae57e-08ae-4ccd-9616-3332802ce290">

| Parameter | Controlled by |
| - | - |
| Number of navigator points (N<sub>nav</sub>) | cv: `nnav` |
| VDS center acceleration (undersampling) factor (R<sub>center</sub>) | cv: `vds_acc0` |
| VDS edge acceleration (undersampling) factor (R<sub>edge</sub>) | cv: `vds_acc1` |
| kz acceleration (undersampling) factor (SOS only) (R<sub>z</sub>) | cv: `kz_acc1` |
| 3D projection mode | cv: `spi_mode` = (1) Stack of spirals, (2) TGA-SPI, (3), 3DTGA-SPI |

#### Controlling the Readout
The following variables control all types of readouts:

| Parameter | Controlled by |
| - | - |
| Echo train length (ETL) | cv: `opetl` |
| Readout type | cv: `ro_type` = (1) FSE, (2) SPGR, or (3) bSSFP |
| Area under each gradient spoiler (as % of kmax) | cv: `crushfac` |

##### FSE mode
<img width="600" src="https://github.com/user-attachments/assets/77e34722-0bd2-45ab-a37b-f20061b934b1">

| Parameter | Controlled by |
| - | - |
| Echo time (TE: time from 90 to center of spiral in/out) (us) | cv: `opte` |
| flip angle (α: refocuser angle) (deg) | cv: `opflip` |

##### SPGR mode

| Parameter | Controlled by |
| - | - |
| RF spoiling (117deg phase cycling | cv: `rfspoil_flag` |

##### bSSFP mode
<img width="600" src="https://github.com/user-attachments/assets/1707c63b-604d-4255-bd45-f4bdc5737dee">


##### EPIC Readout Structure
The following pulse sequence diagram shows how the readout is structured in the EPIC sequence code:

<img width="800" src="https://github.com/user-attachments/assets/45684af3-7832-4e64-a7a8-24c5d8a39dd5">
