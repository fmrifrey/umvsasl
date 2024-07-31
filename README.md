# umvsasl
Written by David Frey for the fMRI lab at the University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our umvsasl project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for ASL (Arterial Spin Labeling) MRI experiments.

## Table of Contents
1. [Introduction](#introduction)
    - [Directory structure](#directory-structure)
2. [Usage](#usage)
    - [Pulse sequence](#pulse-sequence)
        - [Basic sequence control](#basic-sequence-control)
            - [Looping structure](#looping-structure)
            - [K-space sampling scheme](#kspace-sampling-scheme)
        - [Controlling the readout](#controlling-the-readout)
            - [Fast spin echo (FSE) mode](#fse-mode)
            - [Spoiled gradient echo (SPGR) Mode](#spgr-mode)
            - [Balanced steady state free precession (bSSFP) Mode](#bssfp-mode)
            - [EPIC readout structure](#epic-readout-structure)
        - [Controlling the preparation pulses](#prep-pulses)
            - [Fat suppression](#fat-suppression)
            - [Custom prep pulses](#custom-prep-pulses)
            - [PCASL pulse](#pcasl-pulse)
            - [Pre-saturation pulse](#pre-saturation-pulse)
3. [License](#license)
4. [Contact](#contact)

## Introduction

## Usage

### Pulse sequence

#### Basic sequence control

##### Looping structure
The main sequence follows a looping structure over three readout dimensions: frames, shots (through-plane interleaves), and arms (in-plane interleaves), ordered from fastest to slowest, respectively. Each TR (shot interval) contains all preparation pulses & delay times as well as the readout echo train. [The entire kspace volume is encoded](#kspace-encoding-scheme) between the echo train, shots, and arms dimensions. Each frame contains information from one volume in kspace.

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

Each arm interleaves the initial 2D spiral in-plane by uniform rotations as shown on the left of the figure below. Each shot interleaves the projected spiral volume acquired in a single echo echo train. For stack of spirals, this looks like through-plane interleaving along z as shown in the top right of the figure below. For spiral projection imaging, this will interleave the 3D rotations as shown in the bottom right of the figure below.

<img width="600" src="https://github.com/user-attachments/assets/7e6ae57e-08ae-4ccd-9616-3332802ce290">

The following parameters control the kspace trajectory:

| Parameter | Controlled by |
| - | - |
| Number of navigator points (N<sub>nav</sub>) | cv: `nnav` |
| VDS center acceleration (undersampling) factor (R<sub>center</sub>) | cv: `vds_acc0` |
| VDS edge acceleration (undersampling) factor (R<sub>edge</sub>) | cv: `vds_acc1` |
| kz acceleration (undersampling) factor (SOS only) (R<sub>z</sub>) | cv: `kz_acc1` |
| 3D projection mode | cv: `spi_mode` = (1) Stack of spirals, (2) Tiny Golden Angle, (3), 3D Tiny Golden Angle |
| Maximum slew rate (G/cm/s) | cv: `SLEWMAX` |
| Maximum gradient amplitude (G/cm) | cv: `GMAX` |

#### Controlling the readout
This sequence has the ability to run fast spin echo (FSE), spoiled GRE (SPGR), and balanced steady-state free precession (bSSFP) readouts. The following parameters control all types of readouts:

| Parameter | Controlled by |
| - | - |
| Readout type | cv: `ro_type` = (1) FSE, (2) SPGR, or (3) bSSFP |
| Echo train length (ETL) | cv: `opetl` |
| Number of disdaq echoes | cv: `ndisdaqechoes` |
| Number of disdaq echo trains | cv: `ndisdaqtrains` |
| Area under each gradient spoiler (as % of kmax) | cv: `crushfac` |
| Flow compensation | cv: `flowcomp_flag` = (0) off, (1) on |
| FID mode (no gradients) | cv: `kill_grads` = (0) off, (1) on |

Readout type-specific diagrams and parameters are as follows:

##### FSE mode
<img width="600" src="https://github.com/user-attachments/assets/77e34722-0bd2-45ab-a37b-f20061b934b1">

| Parameter | Controlled by |
| - | - |
| Echo time (TE) (us) | cv: `opte` |
| Refocuser flip angle (α) (deg) | cv: `opflip` |

##### SPGR mode
<img width="600" src="https://github.com/user-attachments/assets/d4ccecf0-994f-450e-b598-fdba2a4aaf68">

| Parameter | Controlled by |
| - | - |
| Echo spacing (ESP, Short TR) (us) | cv: `esp` |
| Echo time (TE) (us) | cv: `opte` |
| Tipdown flip angle (α) (deg) | cv: `opflip` |
| RF spoiling (117deg phase cycling | cv: `rfspoil_flag` |

##### bSSFP mode
<img width="600" src="https://github.com/user-attachments/assets/1707c63b-604d-4255-bd45-f4bdc5737dee">

| Parameter | Controlled by |
| - | - |
| Echo spacing (ESP, Short TR) (us) | cv: `esp` |
| Echo time (TE) (us) | cv: `opte` |
| Tipdown flip angle (α) (deg) | cv: `opflip` |

##### EPIC readout structure
The following pulse sequence diagram shows how the readout is structured in the EPIC sequence code:

<img width="800" src="https://github.com/user-attachments/assets/45684af3-7832-4e64-a7a8-24c5d8a39dd5">

note: (fixing any of the values without care will likely cause the sequence to crash)

#### Controlling the preparation pulses
