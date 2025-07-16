---
title: Extracting the magnetopause surface from the simulation
---

Although, the question of the validity of the tangent hypothesis has been tackled by previous studies, further characterization is still necessary in order to derive the limits of the tangent fitting approach. In this chapter we will study how this translates to the X-ray emissivity that has been simulated by the LaTeP model, taking into account the evolution of the orbit and attitude of the satellite. Unlike the single point approach of A. Read {cite}`andi_2024`, we will attempt extracting the full maximum intensity arc of the integrated images and compare this with the projection of the magnetopause tangent that has been extracted from the simulation cubes.


The Tangent Fitting Approach is based on the hypothesis that the maximum intensity arc of the SXI images, correspond to the tangent direction of the magnetopause surface. To verify this claim, we will:
1. Extract the magnetopause surface directly from the emissivity cubes.
2. Project its tangent points to the SXI imager's FOV.
3. Extract the maximum intensity arc and compare it to the projection.


(ch:TH)=
# Extracting the magnetopause surface from the simulation

The extraction of the magnetopause location from the simulation requires that we understand the information that the output cubes can provide, as well as the various ways the magnetopause can be defined and derived from this information. This can be done either utilizing the electromagnetic fields of the MHD input that was given to the LaTeP model, either the X-ray production itself in the MHD or the LaTeP case.
