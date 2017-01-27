# 4D MR Imaging - Self-navigated 4D Cartesian Imaging of Periodic Motion
Acquistion and reconstruction method for MR motion imaging and correction

## Acquisition
A Cartesian 3D random Gaussian ESPReSSo subsampling is proposed which acquires a random combination of k-space lines ky and kz in each repetition, as illustrated in the Figure below. Periodically the central k-space lines are acquired which serve as a self-navigation signal in the reconstruction.
- precompiled 4D MR imaging sequence: CS_Retro (Siemens, VB20P)

## Reconstruction
retrospective motion gating and reconstruction of subsampled data
- self-navigation signal extraction
- gating procedure
- for the CS reconstruction, please refer to https://github.com/thomaskuestner/CS_LAB
 
### Matlab
code included here

### Gadgetron
please refer to https://github.com/thomaskuestner/CS_LAB/reconstruction/gadgetron/CS_LAB_Gadget

--------------------------------------------------------
Please read LICENSE file for licensing details.

Detailed information are available at:
https://sites.google.com/site/kspaceastronauts/motion-correction/4d-mr-imaging
