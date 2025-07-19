# About this work

In this project we present the analysis of the synthetic images produced by the [LATMOS](https://latmos.ipsl.fr/fr/) Test Particle ([LaTeP](https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2024JA032687)) model, in preparation for the [SMILE](https://www.esa.int/Science_Exploration/Space_Science/Smile/Smile_factsheet2) mission. You can either clone the repository of the [source code](https://github.com/notanobis/SMILE-TFA), read the online [jupyter book](https://notanobis.github.io/SMILE-TFA/Abstract.html), or download the document of my Master's thesis directly from <a href="MasterThesis_PB.pdf" target="_blank">here</a>.

Here is a complete table of content with the corresponding Jupyter Notebooks you can find on the GitHub repository:

Introduction
1. [Abstract : Locating the Earth’s magnetopause through X-ray imaging](https://notanobis.github.io/SMILE-TFA/Abstract.html)
2. [Introduction](https://notanobis.github.io/SMILE-TFA/Introduction.html)
    1. [Earth’s Magnetopause](https://notanobis.github.io/SMILE-TFA/Magnetopause.html)
    2. [Detecting the magnetopause](https://notanobis.github.io/SMILE-TFA/Detecting.html)
    3. [Previous work: Fitting methods](https://notanobis.github.io/SMILE-TFA/FittingMethods.html)
Tangent Hypothesis
1. [Extracting the magnetopause surface from the simulation](https://notanobis.github.io/SMILE-TFA/Chapter02.html)
    1. [LaTeP model: Emissivity cube](https://notanobis.github.io/SMILE-TFA/LatepCube.html): `LatepCube.ipynb`
    2. [Shue model fitted to slices](https://notanobis.github.io/SMILE-TFA/ShueToSlice.html): `ShueToSlice.ipynb`
    3. [Radial extraction of full magnetopause](https://notanobis.github.io/SMILE-TFA/MagSurface.html): `MagSurface.ipynb`
    4. [Fitting models to surface](https://notanobis.github.io/SMILE-TFA/FittingModels.html): `FittingModels.ipynb`
    5. [Extracting the magnetopause from the MHD emissivity](https://notanobis.github.io/SMILE-TFA/MHDsurface.html): `MHDsurface.ipynb`
2. [Extracting the maximum intensity curve from the image](https://notanobis.github.io/SMILE-TFA/3DView.html#)
    1. [Image processing](https://notanobis.github.io/SMILE-TFA/MaxIntensityArc.html): `MaxIntensityArc.ipynb`
3. [Characterizing the tangent hypothesis](https://notanobis.github.io/SMILE-TFA/TH.html)
    1. [Projecting the Surface](https://notanobis.github.io/SMILE-TFA/Projection.html): `Pojection.ipynb`
    2. [Results: Characterizing the Tangent Hypothesis](https://notanobis.github.io/SMILE-TFA/TangentHypothesis.html): `TangentHypothesis.ipynb`
