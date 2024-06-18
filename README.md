# asl3dflex
written by David Frey for fMRI lab at University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our asl3dflex project repository. This project contains an end-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL using a flexible spiral projection readout.

## Table of Contents

1. [Introduction](#introduction)
    - [Directory Structure](#directory-structure)
4. [Usage](#usage)
    - [Pulse Sequence Source Code](#pulse-sequence-source-code)
        - [Code Structure](#psd-code-structure)
        - [EPIC Compilation](#psd-compilation)
    - [Scanner Operation](#scanner-operation)
        - [Scanner Setup](#scanner-setup)
        - [Running Basic Sequences](#basic-sequences)
        - [Running ASL Sequences](#asl-sequences)
        - [List of Variables](#scanner-variables)
    - [MATLAB Reconstruction](#matlab-reconstruction)
6. [License](#license)
7. [Contact](#contact)

## Introduction

### Directory Structure
This repository contains several subpackages aimed at different aspects of MRI sequence development and data reconstruction. 

- `/psdsrc`: Contains GE pulse sequence source code.
- `/recon`: Contains MATLAB scripts for reconstructing MRI data.
- `/scanner`: Contains a library of bash scripts and pulse .txt files used for scanning.




