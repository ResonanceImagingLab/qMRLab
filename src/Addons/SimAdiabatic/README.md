# Adiabatic Inversion Pulses 

## Background: What are adiabatic pulses? 
Adiabatic pulses are a special class of RF pulses that can excite, refocus, or **invert** magnetization vectors uniformly (Matt A. Bernstein et al., 2004a). In a RF frequency sweep from one side of resonance to the other, the net rotation of the magnetization vector **$\overrightarrow{M}$** is highly insensitive to changes in B<sub>1</sub> amplitude (Tannús & Garwood, 1997). So, adiabatic pulses rotate **$\overrightarrow{M}$** by a constant 180&deg; flip angle, even when B<sub>1</sub> is extremely inhomogeneous (Tannús & Garwood, 1997). Adiabatic pulses also operate under the **adiabatic passage principle** or **adiabatic condition** which states that **$\overrightarrow{M}$** of a spin system follows the direction of **$\overrightarrow{B_{eff}}$** such that the direction of **$\overrightarrow{B_{eff}}$** does not change much during one period of precession (Matt A. Bernstein et al., 2004a). 

## Adiabatic vs Non-Adiabatic 

### In a standard RF pulse: 
- ${θ = γ\int_0^T{B_1(τ)dτ}}$
- constant carrier frequency
- short duration (0.3-1ms)
- lower B<sub>1</sub> amplitude 
- generally multi-purpose

### In an adiabatic pulse: 
- ${θ≠γ\int_0^T{B_1(τ)dτ}}$
- long duration (10-20ms)
- higher B<sub>1</sub> amplitude (>10μT)
- generally NOT multi-purpose (inversion can't be used for refocus, etc.)

## What are adiabatic inversion pulses? 
Adiabatic inversion pulses are a special class of RF pulses that rotate **$\overrightarrow{M}$** from the +z to -z axis (Matt A. Bernstein et al., 2004a). They will uniformly invert **$\overrightarrow{M}$** across an imaged object even when the B<sub>1</sub> field is spatially non-uniform (Matt A. Bernstein et al., 2004a). These pulses also operate under the **adiabatic condition** and can be displayed visually in the following figure. Note: Figure a-e denotes beginning to end of an adiabatic inversion pulse and Figure f denotes the trajectory of a non-adiabatic pulse.

 ![image](https://github.com/ResonanceImagingLab/qMRLab/assets/154541326/b36dc143-b5d3-4070-8deb-a428cc18debc)

Adiabatic inversion pulses can be displayed as a frequency modulated pulse or a phase modulated pulse: 
- **Frequency Modulated pulse:**
  $B_1(t) = A(t)\mathrm{e}^{-i\int{ω_1(t')dt'}}$
- **Phase Modulated pulse:**
  $B_1(t) = A(t)\mathrm{e}^{-i\phi(t)}$ where $\phi(t)=\int{ω_1(t)dt}$
  
A(t) is defined as the envelop/max amplitude and ω<sub>1</sub>(t) is the frequency sweep 

The functions included in this package are designed to display six different adiabatic inversion pulses including: 
 - hyperbolic secant (Hs1)
 - Lorentz
 - Gaussian
 - Hanning
 - hsn (n=2-8)
 - sin40

The following section will take you step-by-step how this document works. 

## Software Requirements 
Tested on MATLAB_R2022a. May work on earlier versions  

## Step-by-step tutorial 
1. Download qMRLab
2. Add qMRlab to your MATLAB path
   - addpath("c:\matlab\MyFolder")
3. Run `open adiabaticExample.m` in your command window on MATLAB
     
### Within adiabaticExample.m ...
- The beginning section gives a brief overview on what is included within the code and where you can find certain information
- Each of the six pulses are listed twice: the first set is for a 1 pool case and the last set is for the 2 pool case
- This code is designed to run section by section so all pulses are separated within their own section
    - You can also run the entire code itself however your computer will be over run by figures
- The parameters currently in place are the default parameters but this code is designed so you can change around the parameters as you please
    - The default parameters are commented out on the side so if you forget you can always go back to the original
    - These parameters are also defined in the defaultPulseParams.m functions for each respective pulse
- Three different plotting options are available to you:
  1. Plot adiabatic pulse including amplitude and frequency modulation functions
  2. Plot adiabatic inversion pulses by calling the Bloch sim results for 1 pool or 2 pool
  3. Plot RF pulse by removing frequency (ω<sub>1</sub>)
- Option 2 is currently the only option uncommented as it is the main goal of this learning tool but uncomment the other options as you please
- All references used to aid in your understanding of adiabatic inversion pulses and this learning tool are listed at the end of getAdiabaticPulse.m.
    - Under each reference is a list of the important information the selected papers or book chapters cover so you do not need to read through every single reference to find the information you are looking for

### Step-by-step example to run through adiabaticexample.m
#### Example for hyperbolic secant 1 pool 
1. Clear your workspace `clear all`
    - The code is set to clear all parameters at the beginning of each pulse
2. Run the first set of initial parameters from Params.B0 = 3 to Params.NumPools = 1 by highlighting the specified values, right click and select "Evaluate section in command window" 
    -  If you want to see the Default tissue params: `open defaultCortexTissueParams` into command window
3. Run the specified pulse parameters from Params.Inv.PulseOpt.beta = 550 to Params.Inv.shape = 'Hs1' by repeating step 2 
4. Apply the inversion pulse by calling `[inv_pulse, omega1, A_t, ~] = getAdiabaticPulse(Params.Inv.Trf, Params.Inv.shape, Params.Inv`
    - You can run this step by highlighting --> Right click --> "Evaluate section in command window"
    - or copy and paste this line into the command window
    - `open getAdiabaticPulse`
5. Define your time array: `t = linspace(0, Params.Inv.Trf, Params.Inv.nSamples)`
6. Select the desired plotting option you wish to see
    - `plotAdiabaticPulse(t, inv_pulse, Params)`
    - `blochSimCallFunction(inv_pulse, t, A_t, omega1, Params)`
    - `blochSimCallFunction(abs(inv_pulse), t, A_t, omega1, Params)`
7. To open any function: `open` and then the name of the desired function (refer to step 2) OR highlight pulse name-->right click-->select **Open "function name"** 
8. Expected plots for each plotting function are displayed below in order they are listed in step 6
   ##### plotAdiabaticPulse
    <img width="1197" alt="Screenshot 2024-04-14 at 4 07 41 PM" src="https://github.com/ResonanceImagingLab/qMRLab/assets/154541326/3e52fffc-fe2d-45f0-868c-23dc66ade43f">

   ##### blochSimCallFunction Adiabatic Inv
   <img width="1197" alt="Screenshot 2024-04-14 at 4 10 55 PM" src="https://github.com/ResonanceImagingLab/qMRLab/assets/154541326/7ba7cf61-93b2-4173-ae05-7fad2802d0a8">

   ##### blochSimCallFunction RF pulse
   <img width="1197" alt="Screenshot 2024-04-14 at 4 12 04 PM" src="https://github.com/ResonanceImagingLab/qMRLab/assets/154541326/05095348-93ea-4d09-ac46-9d7a7a1af7e3">

## Questions to ask yourself 
&#x2610;  What happens if I increase or decrease beta? 

&#x2610;  What happens if I increase or decrease A0?

&#x2610;  What happens if I increase or decrease Q? 

&#x2610;  Does anything happen if nSamples changes from 512 to 256?

&#x2610;  What happens if Trf increases or decreases? 

&#x2610;  What would happen if the Hsn pulse changes from n=8?

&#x2610;  What if a higher Tesla magnet was used? (ex. 7T)

&#x2610;  What if a lower Tesla magnet was used? (ex. 1.5T) 
