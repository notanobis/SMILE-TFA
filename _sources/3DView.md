# Extracting the maximum intensity curve from the image

As suspected from \autoref{fig:lines} the extraction of the maximum intensity arc is not as trivial as in the smooth-gradient MHD case. This will not change when introducing the \ac{SXI} \ac{FOV} and the satellite's \ac{POV}, but it will improve since the instrument never points directly to the cusps, while image processing techniques can also be applied to extract continuous arcs.

## Constructing the images

To construct the synthetic images of the \ac{SXI} instrument, we use the \href{https://3dview.irap.omp.eu/}{3DView}  open-source tool developed by CDPP \cite{3DView}. Through this tool we can simulate the orbit and attitude of the satellite for a specified interval, as well as the orientation and characteristics of the imaging system in relation to the satellite. We can directly import the emissivity cubes of the simulations, introduce the satellite and the \ac{FOV} of the imager, and compute the integrated flux over each pixel's \ac{LOS}. We therefore get the integrated images from each \ac{POV} of the orbit, according to the \ac{FOV} and resolution of the imager, as shown in \autoref{3DView}. No simulated noise has been added to account for the detection system of the instrument or the background contribution. A typical synthetic image of the magnetopause, generated through this method, is shown in \autoref{fig:typical_image}.

```{figure} ./images/3DView.png
:name: 3DView
3D View CDPP software
```


Along with the images, we can export XML files that contain information about the position and orientation of the satellite, for each time-step. This will be particularly useful in \autoref{section: tangent_hypothesis}.