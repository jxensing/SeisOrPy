# SeisOrPy
Estimating horizontal seismometer sensor orientation using empirical Green's functions computed from ambient seismic noise. 

This is the code and data used in our article: Ensing, J. X., & van Wijk, K. (2019). Estimating the orientation of borehole seismometers 
from ambient seismic noise. Bulletin of the Seismological Society of America, 109(1), 424-432. The ImpulseResponses (or Green's Functions) included here, were prepared as outlined in this article. The notation/naming of the files may be confusing, because the traces are rotated into ZRT as if they were originally in ZNE convention. With the borehoel seismometers, however, the orientation of the N and E traces is unknown, and the orientation of the R and T traces are also unkown. The orientation of the sensors is discovered by running the program.

The method is slightly modified after Zha, Y., S. C. Webb, and W. Menke (2013). Determining the orientations of ocean bottom seismometers 
using ambient noise correlation, Geophys. Res. Lett. 40, no. 14, 3585–3590.

To try it out, download the SeisOrPy folder, and unzip ImpulseResponses.zip, and run SeisOr.py

It will require the following python packages obspy, numpy, and scipy.

