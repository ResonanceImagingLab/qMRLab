# Adiabatic Inversion Pulses 

## What are adiabatic pulses? 
Adiabatic pulses are a special class of RF pulses that can excite, refocus, or **invert** magnetization vectors uniformly (Matt A. Bernstein et al., 2004a). In a RF frequency sweep from one side of resonance to the other, the net rotation of the magnetization vector **M** is highly insensitive to changes in B1 amplitude (Tannús & Garwood, 1997). So, adiabatic pulses rotate **$\overrightarrow{M}$** by a constant 180&deg; flip angle, even when B1 is extremely inhomogeneous (Tannús & Garwood, 1997). Adiabatic pulses also operate under the **adiabatic passage principle** or **adiabatic condition** which states that **$\overrightarrow{M}$** of a spin system follows the direction of **Beff** such that the direction of **$\overrightarrow{Beff}$** does not change much during one period of precession (Matt A. Bernstein et al., 2004a). 

## Adiabatic vs Non-Adiabatic 
### In an adiabatic pulse... 
- ${Θ \not{=} γ\int0^T{B1(τ)dτ}}$
- long duration (10-20ms)
- higher B1 amplitude (>10μT)
- generally NOT multi-purpose (inversion can't be used for refocus, etc.)
  
### In a regular RF pulse... 
- ${Θ = γ\int_0^T{B1(τ)dτ}}$
- AM, constant carrier frequency
- short duration (0.3-1ms)
- lower B1 amplitude (10μT)
- generally multi-purpose

## What are adiabatic inversion pulses? 
Adiabatic inversion pulses are a special class of RF pulses that rotate **M** from the +z to -z axis (Matt A. Bernstein et al., 2004a). They will uniformly invert **$\overrightarrow{M}$** across an imaged object even when the B1 field is spatially non-uniform (Matt A. Bernstein et al., 2004a). These pulses also operate under the **adiabatic condition** and can be displayed visually in the following figure. 

![image](https://github.com/ResonanceImagingLab/qMRLab/assets/154541326/b36dc143-b5d3-4070-8deb-a428cc18debc)

Adiabatic inversion pulses can be displayed as a frequency modulated pulse or a phase modulated pulse: 
- **Frequency Modulated pulse:**
  $B1(t) = A(t)\mathrm{e}^{-i\int{w1(t')dt'}}$
- **Phase Modulated pulse:**
  $B1(t) = A(t)\mathrm{e}^{-i\phi(t)}$ where $\phi(t)=\int{w1(t)dt}$
  
A(t) is defined as the envelop/max amplitude and w1(t) is the frequency sweep 

These functions are designed do display the six different adiabatic inversion pulses including hyperbolic secant (Hs1), lorentz, gaussian, hanning, Hsn (n=2-8) and Sin40. The following section will take you step-by-step how do work this development. 

## Software Requirements 
MATLAB_R2023b

## Step-by-step tutorial 
1. Download all Matlab functions and code within all folders as well as the adiabaticexample.m file.
2. Add all function and code to your path on MATLAB
3. Open adiabaticexample.m

## How does adiabaticexample.m work? 
### Within adiabaticexample.m ...
- The beginning section gives a brief overview on what is included within the code and where you can find certain information
- Each of the six pulses are listed twice: the first set is for a 1 pool case and the last set is for the 2 pool case
- This code is designed to run section by section so all pulses are separated within their own section
    - You can also run the entire code itself however your computer will be over run by figures
- The parameters currently in place are the default parameters but this code is designed so you can change around the parameters as you please
    - The default parameters are commented out on the side so if you forget you can always go back to the original
    - These parameters are also defined in the defaultpulseparams.m functions for each respective pulse
- Three different plotting options are available to you:
  1. Plot adiabatic pulse including amplitude and frequency modulation functions
  2. Plot adiabatic inversion pulses by calling the Bloch sim results for 1 pool or 2 pool
  3. Plot RF pulse by removing frequency (omega1)
- Option 2 is currently the only option uncommented as it is the main goal of this learning tool but uncomment the other options as you please
- Still more to be added to this readme 
