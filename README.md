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
            - [K-space Encoding Scheme](#kspace-encoding-scheme)
        - [Controlling the Readout](#controlling-the-readout)
            - [Fast Spin Echo (FSE) Mode](#fse-mode)
            - [Spoiled Gradient Echo (SPGR) Mode](#spgr-mode)
            - [Balanced Steady State Free Precession (bSSFP) Mode](#bssfp-mode)
            - [High-Level Readout Structure](#high-level-readout-structure)
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

<img width="300" src="https://github.com/user-attachments/assets/3e301dd8-2dc4-4241-bfb1-1efa7e5062d6">

| Variable | Description | Controlled by |
| - | - | - |
| M | number of frames | cv: `nframes` |
| N | number of shots | cv: `opnshots` |
| P | number of arms | cv: `narms` |
| Shot interval | - | cv: `optr` |
| Sequence duration | - | n/a |

##### Kspace encoding scheme
This sequence uses 3D stack of spirals and spiral projection imaging (SPI) to encode kspace.
An initial 2D [variable-density spiral](https://www.researchgate.net/publication/5925197_Fast_3D_imaging_using_variable-density_spiral_trajectories_with_applications_to_limb_perfusion) is generated and then transformed in 3D using rotation matrices.
Each readout dimension also transforms the kspace trajectory such that each frame contains all information for a single volume in kspace.
Each arm interleaves the initial 2D spiral in-plane by uniform rotations as shown on the left of the figure below.
Each shot interleaves the projected spiral volume acquired in a single echo echo train. For stack of spirals, this looks like through-plane interleaving along z as shown in the top right of the figure below. For spiral projection imaging, this will interleave the 3D rotations as shown in the bottom right of the figure below.

<img width="500" src="https://github.com/user-attachments/assets/7e6ae57e-08ae-4ccd-9616-3332802ce290">

| Variable | Description | Controlled by |
| - | - | - |
| N<sub>nav</sub> | number of points sampled at center of kspace | cv: `nnav` |
| R<sub>center</sub> | VDS center acceleration (undersampling) factor | cv: `vds_acc0` |
| R<sub>edge</sub> | VDS edge acceleration (undersampling) factor | cv: `vds_acc1` |
| R<sub>z</sub> | kz acceleration (undersampling) factor (SOS only) | cv: `kz_acc1` |
| Projection mode | (1) Stack of spirals, (2) TGA-SPI, (3), 3DTGA-SPI | cv: `spi_mode` |

#### Controlling the Readout
The following variables control all types of readouts:
| Variable | Description | Controlled by |
| - | - | - |
| ETL | echo train length | cv: `opetl` |
| Readout type | (1) FSE, (2) SPGR, (3) bSSFP | cv: `echo_mode` |

##### FSE mode
| Variable | Description | Controlled by |
| - | - | - |
| TE | echo time (us) | cv: `opte` |
| α | flip angle (deg) | cv: `opflip` |
| Crush factor | crusher area (as % of kmax) | cv: `crushfac` |
| RF spoiling | option to using 117deg phase cycling | cv: `rfspoil_flag` |
