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
