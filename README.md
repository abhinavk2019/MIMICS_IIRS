# MIMICS_IIRS
Since we are modelling a naturally occurring phenomena there are a huge number of input parameters to consider. The major one are listed below:

System parameters Radar Frequency: Measured in Ghz, this is the operating frequency of the radar, kept at 1.25 GHz Wavelength: Appropriate wavelength has to be selected in order for backscattering to occur, we utilize L band, C band and X band for this. Incident Angle: Backscattering is quite dependent on the incident angle of the wave.

Environment parameters: Canopy density: modelled as tree per square metre, for eg 0.5 in one of the simulation run Trunk height and diameter: Leaf density:modelled as leaves per cubic meter, for eg 48.46

Dielectric value of : Leaves Needles Branches Trunks Ground Surface

The output of the simulation is the value of backscattering coefficient for various values of angle of incidence for different wavelengths.

We have 4 different types of polarized output to consider, HH,VV, HV and VH.

An example : following parameters are considered

Canopy Density =  0.11 Trees per square meter, Crown Height =  2.00 meters, Trunk Height =  8.00 meters, Vegetation Temperature =  20.00 degrees C.
 Surface Dielectric =  5.680 -j 1.130, Soil RMS Roughness =  0.45 cm, Soil Correlation Length =  18.75 cm, Model Type = Physical Optics               
 Trunk Dielectric = 28.190 -j 7.830 Trunk Diameter = 24.000 cm
 Needle Gravimetric Moisture = 0.800, Needle Dry Density = 0.377 cm, Needle Density =  0.8300E+03 needles per cubic meter, Needle Diameter = 0.100 cm, Needle Length = 13.001 cm
 Branch Dielectric = 28.190 -j 7.830, Branch Density =   4.100 branches per cubic meter, Branch Diameter =  0.700 cm, Branch Length =   0.750 meters

To run the simulation, you need to run the 'run.sh' shell file in the code folder.
Input data can be changed in the by changing the input values in the related file in the Data folder.
