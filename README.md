# asl3dflex
written by David Frey for fMRI lab at University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our asl3dflex project repository. This project implements an end-to-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL.

## Table of Contents

1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
2. [Usage](#usage)
    - [Scanner Operation](#scanner-operation)
        - [Scanner Setup](#scanner-setup)
        - [Running Basic Sequences](#running-basic-sequences)
        - [Running ASL Sequences](#running-asl-sequences)
        - [List of User Control Variables](#list-of-user-control-variables)
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
- `/scanner`: Contains a library of bash scripts and pulse .txt files used for scanning. These files should be on the path of the MRI system host (generally in /usr/g/bin).

## Usage

### Scanner Operation

#### Scanner Setup

Prerequisites:
1. GE MRI system compatible with EPIC MR30.1
2. All executables and libraries on path:
    - Compiled sequences: asl3dflex (host executable file) and asl3dflex.psd.ice (tgt executable file)
    - On-scanner directory: /scanner
