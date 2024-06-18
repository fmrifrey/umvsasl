# asl3dflex
by David Frey
written for fMRI lab at University of Michigan, PI: Luis Hernandez-Garcia

Welcome to our asl3dflex project repository. In this project, we create a end-end pipeline for pulse sequence development and image reconstruction for velocity-selective ASL using a flexible spiral projection readout

## Table of Contents

1. [Introduction](#introduction)
2. [Directory Structure](#directory-structure)
3. [Installation](#installation)
4. [Usage](#usage)
    - [Pulse Sequence Source Code](#pulse-sequence-source-code)
    - [MATLAB Reconstruction](#matlab-reconstruction)
    - [Scanner operation + Bash Scripts and Pulse Files](#bash-scripts-and-pulse-files)
5. [Contributing](#contributing)
6. [License](#license)
7. [Contact](#contact)

## Introduction

This repository contains several subpackages aimed at different aspects of MRI sequence development and data reconstruction. 

- `/psdsrc`: Contains GE pulse sequence source code.
- `/recon`: Contains MATLAB scripts for reconstructing MRI data.
- `/scanner`: Contains a library of bash scripts and pulse .txt files used for scanning.

## Directory Structure

