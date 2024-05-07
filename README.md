# AM_JointInversion

This is an inversion for 1D seismic structure using an adaptive meshing scheme from reeciver functions and surface waves - essentially the classic Julia et al., 2000 routine but redesign to go to much higher frequencies for sharp structure. This work is currently pending at NSF Geophysics (ie I wrote it for a proposal). 

Advances over existing routine include
 - Inverison of multi-taper deconvolved recevier functions, not time domain, which preserves higher frequencies.
 - Scaling the windows in extended time multitaper deconvolution to better preserve high frequency signals
 - Adaptive mesh to mimick the behavior of transdimensional inversions, except this is much faster. Run times are typically ~<6 hours with full errors, which keep in mind is for high-frequency models requires hundreds of layers in the model for the receiver function calculation.
 - Iteratively updated free-surface rotation based on the Vs you get at the top of the model
 - Vp/Vs variations solved for near the top of the model
 - Weighting of the RF errors by magnitude of the harmonic terms - to avoid mapping anisotropy in the isotropic model.
 - A dampening scheme on the second derivative that allows big contrasts to grow (ie, a Cauchy function log(1 + weight*Vs'') instead of a Gaussian function weight*Vs'')

Things that could be added
 - Ellipicity, trivial addition. 
 - Love waves and a radial anisotropy term
 - A full Vp/Vs model where the reverbs are also modeled (I have tried this, but it made things slow)
 - A dereverb filter for sediments. Right now, if you try to invert a lot of sites, the high-f terms are all reverbs. You can model these well, sometimes. Sometimes, you get garbage and more work needs to be done to handle these. Lots of good targets with terribly reverbs though!
