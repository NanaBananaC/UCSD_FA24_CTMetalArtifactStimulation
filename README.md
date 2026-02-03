# Metal Artifact Simulation and Reduction in CT

## Overview
This repository contains MATLAB implementations for **simulating metal artifacts in CT images** and applying **metal artifact reduction (MAR)** techniques. The project was completed as part of the **BENG 280A / ECE 207 Midterm Project** at **UC San Diego (Fall 2024)**.

Metal artifacts in CT arise mainly from:
- **Beam hardening**
- **Photon starvation**

This repo explores both the simulation of these effects and simple computational approaches to reduce them.

---

## Files in this Repo

- **`simulate_v2.m`** — Simulates metal artifacts in CT (photon starvation and beam hardening effects).
- **`Sinogram Interpolation MAR.m`** — Implements sinogram-based metal artifact reduction via interpolation.
- **`Simplified Iteration MAR.m`** — Implements a simplified iterative metal artifact reduction method.
- **`BENG280A_ECE207 Midterm.pdf`** — Project report/slides with background, methods, and results.

---

## Methods (High Level)

We implemented two MAR approaches:
- **Sinogram Interpolation MAR**: Identifies metal regions in the sinogram and interpolates corrupted areas before reconstruction.
- **Simplified Iterative MAR**: Simulates metal artifacts and subtracts them iteratively from the image.

---

## Team & Contact

- **Benjia Zhang** - b9zhang@stanford.edu
- **Juo-Hsuan Chang** — juc077@ucsd.edu  
- **Iris Zaretzki** -  izaretzki@ucsd.edu

---

## Reference

See **`BENG280A_ECE207 Midterm.pdf`** for full background, theory, methods, and results.
