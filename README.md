# umvsasl
Written by David Frey for the fMRI lab at the University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our umvsasl project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL (Arterial Spin Labeling).

## Table of Contents
1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
2. [Usage](#usage)
    - [Pulse Sequence](#pulse-sequence)
        - [List of User Control Variables](#list-of-user-control-variables)
        - [Scanner Setup](#scanner-setup)
            - [Prerequisites](#prerequisites)
            - [Selecting the Pulse Sequence](#selecting-the-pulse-sequence)
    - [MATLAB Reconstruction](#matlab-reconstruction)
        - [Requirements](#requirements)
        - [Basic NUFFT Reconstruction Example](#basic-nufft-reconstruction-example)
        - [SENSE Reconstruction Example](#sense-reconstruction-example)
    - [EPIC Source Code](#epic-source-code)
        - [Pulse Sequence Code Structure](#pulse-sequence-code-structure)
        - [EPIC Compilation](#epic-compilation)
3. [License](#license)
4. [Contact](#contact)

## Introduction

This repository contains the implementation of an end-to-end pipeline for developing pulse sequences and reconstructing images for velocity-selective ASL. The project is designed to be used within the context of the University of Michigan fMRI Lab.

### Directory Structure
The repository is organized into several subdirectories:
- `/psdsrc`: Contains GE EPIC pulse sequence source code. This code should be used only on the machine dedicated to EPIC compilation.
- `/recon`: Contains MATLAB scripts for reconstructing MRI data. This code should be used only on the machine dedicated to reconstruction.
- `/scanner`: Contains a library of bash scripts and pulse .txt files used for scanning. These files should be on the path of the MRI system host (usually `/usr/g/bin`).
- `/docs`: Contains files documenting validation testing of our project.
    - `/tests`: Contains notes from experimental validation testing of our sequence.
- `/assets`: Contains assets for documentation on GitHub.

## Usage

### Pulse Sequence

Our sequence uses several ASL preparation pulses prior to a spoiled gradient echo (SPGR) acquisition. The sequence acquires data using an interleaved variable density 3D stack of spirals k-space trajectory.

A top-down view of our pulse sequence is shown in the figure below:
![image](https://github.com/fmrifrey/umvsasl/assets/96143939/69787155-f705-46be-83b8-1e5875561683)

From the bottom level to the top level, we have multiple nested loops:

#### Echoes
At each "echo," a single spiral arm for a particular kz plane is acquired using a spoiled GRE. The RF excitation pulse includes a pre-spoiler to spoil any transverse magnetization from the previous readout. The RF pulse is also slab-selective, exciting the prescribed image volume. The number of echoes in an echo train is controlled by the user control variable (CV): `opetl`. The echo time can be controlled by the user CV: `opte`, which describes the time from the center of the RF pulse to the beginning of the spiral waveform.

A diagram of each echo is shown below with named waveform objects and timing variables:
![image](https://github.com/fmrifrey/umvsasl/assets/96143939/8445e86c-56c0-4116-bcf6-332b4710d7ad)

#### Frames
At each frame, a single shot of kz planes for a single spiral arm is acquired. The number of frames in the image is controlled by the user CV: `opnframes`. The reconstruction will handle each frame as a separate image volume and concatenate them into a time series. If prep pulses are being modulated, they will do so along the this dimension.

#### Shots & arms
At each shot, multiple frames for a single undersampled kz volume and one spiral arm is acquired. Each shot interleaves the stack of spirals trajectory along kz in a center-out fashion. Each arm represents a full stack of a single spiral arm at a specific rotation about the kz axis.

The kspace interleaving pattern is shown in the diagram below:

### List of User Control Variables
| CV | Description |
| --- | --- |
| SLEWMAX | maximum allowed slew rate in G/cm/s |
| GMAX | maximum allowed gradient amplitude in G/cm |
| RFMAX | maximum allowed RF amplitude in mG |
| opflip | flip angle of each rf tipdown |
| optr | time between each SPGR echo train |
| esp | time between each echo in SPGR echo train |
| opetl | number of echoes in each SPGR echo train |
| opnshots | number of kz interleaves in stack of spirals readout |
| narms | number of interleaved spiral arms in each kz partition |
| nframes | number of frames to acquire |
| ndisdaqtrains | number of disdaq echo trains to play at the beginning of entire scan |
| ndisdaqechoes | number of disdaq echoes to play at beginning of each SPGR echo train |
| fatsat_flag | option to turn on the fat saturation pulse directly before each SPGR echo train |
| rfspoil_flag | option to use RF phase cycling to spoil each excitation |
| flowcomp_flag | option to add additional flow compensation gradients on each kz blip |
| crushfac | crusher amplitude factor in cycles/vox (i.e. dk_crusher = crushfac * kzmax) |
| kill_grads | option to kill the gradients in the readout to acquire FIDS |
| nnav | number of navigator points at beginning of spiral readout |
| kz_acc | undersampling (acceleration) factor along kz |
| vds_acc0 | undersampling factor at center of variable density spiral |
| vds_acc1 | undersampling factor at edge of variable density spiral |
| presat_flag | option to include pre-saturation pulse at beginning of each TR |
| presat_delay | delay between pre-saturation pulse and next object in each TR |
| prep*_id | id number for selected pulse corresponding to directory in /aslprep/pulses (0 = no pulse/delay) |
| prep*_pld | time between end of selected prep pulse and next object |
| prep*_rfmax | rf amplitude of selected prep pulse |
| prep*_gmax | gradient amplitude of selected prep pulse |
| prep*_mod | modulation pattern for selected prep pulse (1 = L/C, 2 = C/L, 3 = always L, 4 = always C) |
| prep*_tbgs* | time from last object to selected background suppression pulse |

### Scanner Setup

#### Prerequisites
1. GE MRI system compatible with EPIC MR30.1
2. All executables and libraries on the path (usually `/usr/g/bin`):
    - Compiled sequences: `umvsasl` (host executable file) and `umvsasl.psd.ice` (target executable file)
    - Pulse/scripts/schedules library: `/scanner`

#### Selecting the Pulse Sequence
To select the pulse sequence, choose any basic sequence (e.g., 3-plane localizer) and then go to imaging options:
[Provide detailed steps or screenshots if possible]

### MATLAB Reconstruction

#### Requirements
[List system requirements and dependencies here]

#### Basic NUFFT Reconstruction Example
[Provide step-by-step instructions for a basic NUFFT reconstruction]

#### SENSE Reconstruction Example
[Provide step-by-step instructions for a SENSE reconstruction]

### EPIC Source Code

#### Pulse Sequence Code Structure
[Provide detailed explanation of the source code structure]

#### EPIC Compilation
[Provide step-by-step compilation instructions]

## License
[Include license information here]

## Contact
[Provide contact information here]
