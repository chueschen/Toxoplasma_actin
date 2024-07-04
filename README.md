# Toxoplasma_actin
COMSOL Multiphysics and MATLAB code for Hueschen et al., "Emergent actin flows explain diverse parasite gliding modes."

For COMSOL Multiphysics simulations (.mph), see COMSOLGuide.pdf for an orientation to the files in this repo. For a more thorough introduction to getting starting on understanding and working with these and similar COMSOL simulations, we have put together a tutorial (HowToGuide.pdf at https://github.com/RPGroup-PBoC/wildebeest_herds) associated with Hueschen, Dunn, Phillips, 2023 (https://journals.aps.org/pre/abstract/10.1103/PhysRevE.108.024610). 

In that study, we present a formulation of the Toner-Tu flocking theory that uses the finite element method to solve the governing equations on arbitrary curved surfaces. That general theoretical study was inspired by the particular biological puzzle of self-organization actin filaments that is the focus of the present paper.


For MATLAB code (.m), the following scripts correspond to the following analyses and figures.
1. PatchGlidingDuration.m:calculating patch gliding duration and number of reversals (Supp. Fig. S1B)
2. ActinMyosin_uTrackAnalysis: analyzing uTrack (Jaqaman et al., Nature Methods 2008) tracks of actin and myosin speckles (Supp. Fig. S2B-D)
3. FastDirectionalActinTracking.m: analyzing fast, directional actin speckle tracks (Fig. 2D-F)
4. segmentAlignManualTracks.m: function called within FastDirectionalActinTracking.m
5. JaspProtrusionSpeeds.m: calculating and plotting Speeds of jasplakinolide-induced actin protrusions (Supp. Fig. S4)
6. MotilityModeFractions.m: analyzing frequency of occurence of different Toxoplasma motility modes for different concentrations of jasplakinolide (Supp. Fig. S8)
