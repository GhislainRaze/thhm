# Two-harmonic homotopy method

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14775562.svg)](https://doi.org/10.5281/zenodo.14775562)

## Description

This repository contains Matlab code for minimal working examples using the two-harmonic homotopy method (THHM) to uncover isolated secondary resonances [1].

A basic pseudoarclength continuation procedure is given in `simpleContinuation.m` to solve the equations of motion with a harmonic balance formalism [2]. It may require manual tuning for other cases. A more robust continuation code/toolbox is recommended for a more advanced use.

The method can be applied by running the following examples:
* `duffingOscillator_ConstantForce.m`: constant-force THHM applied to a Duffing oscillator (1:3 subharmonic resonance by default)
* `duffingOscillator_ConstantAmplitude.m`: constant-amplitude THHM applied to a Duffing oscillator (1:3 subharmonic resonance by default)
* `twoDof.m`: constant-amplitude THHM applied to a two-degree-of-freedom system (1:2 subharmonic resonance by default)
* `beam.m`: constant-amplitude THHM applied to a cantilever beam with clearance contact (1:3 and 1:5 subharmonic resonances by default)

## References
[1] Raze G, Kerschen G. 2024 A two-harmonic homotopy method to connect a primary resonance to its secondary resonances. Proceedings of the Royal Society A 481, 20240875. doi: [10.1098/rspa.2024.0875](https://doi.org/10.1098/rspa.2024.0875). Also available [here](https://hdl.handle.net/2268/330717).

[2] Detroux T, Renson L, Masset L, Kerschen G. 2015 The harmonic balance method for bifurcation analysis of large-scale nonlinear mechanical systems. Computer Methods in Applied Mechanics and Engineering 296, 18–38. doi: [10.1016/j.cma.2015.07.017](https://doi.org/10.1016/j.cma.2015.07.017).

## Disclaimer
This code comes with no warranty.

## License
<p xmlns:cc="http://creativecommons.org/ns#" >This work is licensed under <a href="http://creativecommons.org/licenses/by/4.0/?ref=chooser-v1" target="_blank" rel="license noopener noreferrer" style="display:inline-block;">CC BY 4.0<img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/cc.svg?ref=chooser-v1"><img style="height:22px!important;margin-left:3px;vertical-align:text-bottom;" src="https://mirrors.creativecommons.org/presskit/icons/by.svg?ref=chooser-v1"></a></p>

