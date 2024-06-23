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

At the top-level, we loop through temporal frames. At each frame, all data for a single image volume is collected. Within each frame, we acquire multiple spiral arms of the 3D stack of spirals. Between each arm, the stack of spirals are rotated by uniform rotations, as shown below:


