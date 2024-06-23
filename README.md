# umvsasl
Written by David Frey for the fMRI lab at the University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our umvsasl project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL (Arterial Spin Labeling).

## Table of Contents
1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
2. [Usage](#usage)
    - [Pulse Sequence](#pulse-sequence)
        - [Pulse Sequence Structure](#pulse-sequence-structure)
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
![pulsesequencediagram](https://github.com/fmrifrey/umvsasl/assets/96143939/4bc1e478-81b4-4090-836e-895f8e85dd65)

From the bottom level to the top level, we have multiple nested loops:

#### Echoes
At each "echo," a single spiral arm for a particular kz plane is acquired using a spoiled GRE. The RF excitation pulse includes a pre-spoiler to spoil any transverse magnetization. The RF pulse is also slab-selective, exciting the prescribed image volume. The number of echoes in an echo train is controlled by the user control variable (CV): `opetl`. The echo time can be controlled by the user CV: `opte`, which describes the time from the center of the RF pulse to the beginning of the stack of spiral gradients.

#### Shots
At each "shot," several kz planes are acquired for each echo.

#### Arms
[Provide more detailed explanation or diagram here if possible]

#### Frames
At each frame, all data for a single image volume is collected, including all kz-encoding steps and arms of the spiral required for the k-space sampling specifications. The number of frames in the image is controlled by the user CV: `opnframes`. The reconstruction will handle each frame as a separate image volume and concatenate them into a time series.

Each loop index at a given time also represents a specific transformation in 3D k-space of an initial 2D spiral sampling trajectory. Each "arm" represents a full stack of a single spiral arm at a specific rotation about the kz axis.

### List of User Control Variables
[Provide a detailed list and explanation of user control variables here]

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
