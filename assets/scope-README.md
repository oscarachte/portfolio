# Oscilloscope Signal Digitization

**EPFL · Programmation pour ingénieurs · BA3 · 4 ECTS · Team of 3**

**MATLAB · C · LabVIEW · Image Processing · Signal Reconstruction**

## Overview

Developed a multi-language workflow to extract and reconstruct signal traces from oscilloscope images. The project combined **LabVIEW** for workflow integration, **C** for pixel-level data extraction, and **MATLAB** for geometric correction, signal fitting, visualization, and PDF export.

## Results

The workflow transforms a distorted oscilloscope graph into a corrected, reconstructed set of signal traces with estimated waveform parameters.

| Distorted input | Reconstructed signals |
| --- | --- |
| ![Distorted oscilloscope graph](images/distorted_input.png) | ![Corrected and fitted signal traces](images/reconstructed_signals.png) |

## Technical Approach

- **LabVIEW:** Coordinated image processing and the execution of external routines.
- **C:** Parsed binary pixmap data, identified trace colors and reference corners, and exported pixel coordinates.
- **MATLAB:** Ordered reference corners, calculated a projective transformation, mapped traces into corrected coordinates, and normalized the signal axes.
- **Signal analysis:** Compared sinusoidal and cubic-polynomial fits using RMSE; extracted waveform parameters when a sinusoidal fit was selected.
- **Output:** Reconstructed plots with fitted traces and exported them as PDFs.

## My Work

Primarily developed the **MATLAB signal reconstruction and analysis** components, and contributed to the **C-based extraction** code. The project was completed collaboratively by three students.

## Repository Structure

- `matlab/` — Signal transformation, fitting, plotting and export scripts.
- `c/` — Pixel extraction source code.
- `labview/` — Original LabVIEW virtual instruments (`.vi`).

## Notes

This repository preserves the original course project source files. The workflow depends on LabVIEW integration and expected input data; it is **not packaged as a standalone application**. The MATLAB fitting step uses the `fit` function and may require the Curve Fitting Toolbox. Some axis scales in the scripts are configured for the original exercise setup.
