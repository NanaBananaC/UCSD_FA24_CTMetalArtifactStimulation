# Metal Artifact Simulation and Reduction in CT

## Overview
This project investigates **metal artifacts in computed tomography (CT)** and implements computational methods to **simulate** and **reduce** these artifacts. Metal implants (e.g., dental fillings, hip replacements) can create severe streaking and distortion in CT images due to physical effects such as **beam hardening** and **photon starvation**.

This work was completed as part of the **BENG 280A / ECE 207 Midterm Project**.

---

## What Causes Metal Artifacts?

**Beam Hardening**
- X-ray beams are polychromatic (multiple energies).
- Dense metals absorb lower-energy photons more strongly, shifting the beam spectrum.
- Results in dark streaks and bright edges in reconstructed images.

**Photon Starvation**
- Metals absorb or scatter most photons.
- Insufficient photon counts at the detector create streaking and speckle noise.

---

## Methods

### 1) Simulation of Metal Artifacts

**Photon Starvation Simulation**
- Convert HU image to attenuation coefficients (μ).
- Generate a sinogram using the Radon transform.
- Add metal attenuation.
- Reconstruct the image using inverse Radon transform.

**Beam Hardening Simulation**
- Model multiple energy levels instead of a single monochromatic beam.
- Generate and combine multiple sinograms.
- Reconstruct a polychromatic CT image showing beam hardening artifacts.

---

### 2) Metal Artifact Reduction (MAR)

**A. Sinogram Interpolation MAR**
1. Segment metal regions via thresholding.
2. Forward project to obtain metal-only sinogram.
3. Interpolate corrupted regions in the original sinogram.
4. Reconstruct the corrected image.

**B. Simplified Iterative MAR**
1. Segment metal regions.
2. Simulate metal artifacts.
3. Subtract simulated artifacts from the original image.
4. Iterate for refinement.

---

## Results
Both methods reduced metal streak artifacts:
- Sinogram interpolation is simpler but may blur details near metal.
- Iterative MAR better preserves structure but requires tuning.

---

## Challenges
- Generalizing to different CT datasets.
- Choosing appropriate thresholds.
- Selecting optimal iteration numbers.
- Balancing artifact removal with anatomical fidelity.

---

## Future Work
- Quantitative evaluation using SD, CNR, and artifact noise metrics.
- Explore deep learning–based inpainting or MAR methods.
- Extend to 3D CT volumes.

---

## Team
- **Benjia Zhang**
- **Juo-Hsuan Chang**
- **Iris Zaretzki**
