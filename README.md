# GCPeakDetection
Peak detection code for gas chromatography data. 

The main code of peak detection is based on Vivo-TGruyols et al. [1]. Some modifications were made to the published work, these were noted in the code. The parameters to use for different data sets are in the file "gcmsProperties.m"

Create the following subdirectories in the results folder: 
* 0fullResults
* mz%d, where %d is the value of m/z, create one for each m/z in the data set. 

Also include in this package is code "readQTOF.m" to convert from XML file to CSV file. For examples of XML files (and the .D files they were converted from) see https://doi.org/10.25919/5b5e699817220

MATLAB code from other sources included in this package: 
* baseCorrALS.m : correcting baseline drift of data using asymmetric least square. [2]
* sgolay.m and sgolayfilt.m : Savitzky-Golay smoothing and fit written by Paul Kienzle and Pascal Dupuis released under GNU GPLv3




[1] Vivo-Truyols et al. "Automatic program for peak detection and deconvolution of multi-overlapped chromatographic signals, Part I: Peak detection", Journal of Chromatography A, vol. 1096 (2005) pp 133-145

[2] P. H. Eilers, H. F. Boelens, Baseline correction with asymmetric least squares smoothing, Leiden University Medical Centre Report 1 (1) (2005) 5.