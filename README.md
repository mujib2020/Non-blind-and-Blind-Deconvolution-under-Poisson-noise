
  

======================================================================
REFERENCE:
======================================================================   
Chowdhury, M. R. and Qin, J. and Lou, Y.; Non-blind and Blind Deconvolution under Poisson Noise using Fractional-order Total Variation, Journal of Mathematical Imaging and Vision, 2020



======================================================================
SOFTWARE
======================================================================   

SOFTWARE REVISION DATE:

       SEPTEMBER 2020

SOFTWARE LANGUAGE:

       MATLAB R2019b

======================================================================
PACKAGE
======================================================================

The directory contains the following files

README         	  : This file
demo_NB.m         : example of how to use the FOTV deblurring method (non-blind)
demo_blind.m      : example of how to use the FOTV blind deconvolution 
demo_blind_vs_NB  : non-blind Vs blind as Fig.6 and Fig.7 in the paper

--------------------------------------------------------
data              : This folder contains two test images.

--------------------------------------------------------
utilities	  : This folder contains the following functions.
FOTVDeblur_NB.m   : function implementing the FOTV deblurring method (non-blind)
FOTV_deconv_blind : function implementing the FOTV deblurring method (blind)
EM_Blind_Deconv   : function implementing the EM deblurring method (blind)
defDDt.m	  : function of fractional derivatives
conv2fft          : function of convolution to have a valid boundary (by P. Favaro)
PSNR.m		  : function of peak signal-to-noise ratio (PSNR)


======================================================================

If you have any questions, feel free to email at mujib.chowdhury@utdallas.edu or mrc.firoj@gmail.com
