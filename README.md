# umvsasl
written by David Frey for fMRI lab at University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our umvsasl project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL.

## Table of Contents
1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
2. [Usage](#usage)
    - [Running the Pulse Sequence](#scanner-operation)
        - [Scanner Setup](#scanner-setup)
        - [Pulse Sequence Structure](#pulse-sequence-structure)
        - [List of User Control Variables](#list-of-user-control-variables)
        - [Scanner Operation Tips](#scanner-operation-tips)
    - [MATLAB Reconstruction](#matlab-reconstruction)
        - [Requirements](#recon-requirements)
        - [Basic NUFFT Reconstruction Example](#recon-example)
        - [SENSE Reconstruction Example](#sense-example)
    - [Pulse Sequence Source Code](#pulse-sequence-source-code)
        - [Pulse Sequence Code Structure](#pulse-sequence-source-code-structure)
        - [EPIC Compilation](#epic-compilation)
3. [License](#license)
4. [Contact](#contact)

## Introduction

### Directory Structure
The repository is broken up into several subdirectories:
- `/psdsrc`: Contains GE EPIC pulse sequence source code. This code should be used only on the machine dedicated to EPIC compilation.
- `/recon`: Contains MATLAB scripts for reconstructing MRI data. This code should be used only on the machine dedicated to reconstruction.
- `/scanner`: Contains a library of bash scripts and pulse .txt files used for scanning. These files should be on the path of the MRI system host (usually /usr/g/bin).
- `/docs`: Contains files documenting validation testing of our project
    - `/tests`: Contains notes from experimental validation testing of our sequence
- `/assets`: Contains assets for documentation on GitHub

## Usage

### Running the Pulse Sequence

#### Scanner Setup

##### Prerequisites
1. GE MRI system compatible with EPIC MR30.1
2. All executables and libraries on path (usually /usr/g/bin):
    - Compiled sequences: umvsasl (host executable file) and umvsasl.psd.ice (tgt executable file)
    - Pulse/scripts/schedules library: /scanner
  
##### Selecting the pulse sequence:
To select the pulse sequence, choose any basic sequence (i.e. 3-plane localizer) and then go to imaging options:

#### Pulse Sequence Structure:
Our sequence uses several ASL preperation pulses prior to a spoiled gradient echo (SPGR) acquisition. Our sequence acquires data using an interleaved variable density 3D stack of spirals kspace trajecory.

A top-down view of our pulse sequence is shown in the figure below:
![pulsesequencediagram](https://github.com/fmrifrey/umvsasl/assets/96143939/4bc1e478-81b4-4090-836e-895f8e85dd65)

From bottom level to top level, we have multiple nested loops:

##### - Echoes
At each "echo", a single spiral arm for a particular kz plane is acquired using a spoiled GRE. The RF excitation pulse object includes a pre-spoiler to spoil any transverse magnetization. The RF pulse is also slab-selective, exciting the prescribed image volume. The number of echoes in an echo train is controlled by user CV: `opetl`. The echo time can be controlled by user CV: `opte`, which describes the time from the center of the RF pulse and the beginning of the stack of spirals gradients.

##### - Shot Loop
At each "shot", several kz planes are acquired for each echo

##### - Arm Loop

##### - Frame Loop
At each frame, all data for a single image volume is collected, including all kz-encoding steps and arms of the spiral required for the k-space sampling specifications. The number of frames in the image is controlled by the User CV: `opnframes`. The reconstruction will handle each frame as a seperate image volume and concatenate them into a timeseries.

##### Stack of Spirals Interleaving


