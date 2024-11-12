# Two-harmonic homotopy method

## Description

This repository contains Matlab code for minimal working examples using the two-harmonic homotopy method (THHM) to uncover isolated secondary resonances [1].

A basic pseudoarclength continuation procedure is given in `simpleContinuation.m` to solve the equations of motion with a harmonic balance formalism [2]. It may require fine tuning, and for a more advanced use, a more robust continuation code/toolbox is recommended.

The method can be applied by running the following examples:
* `duffingOscillator_ConstantForce`: constant-force THHM applied to a Duffing oscillator
* `duffingOscillator_ConstantAmplitude`: constant-amplitude THHM applied to a Duffing oscillator
* `twoDof`: constant-amplitude THHM applied to a two-degree-of-freedom system
* `beam`: constant-amplitude THHM applied to a cantilever beam with clearance contact

## References

[1] Raze G, Kerschen G. 2024 A two-harmonic homotopy method to connect a primary resonance to its secondary resonances. Preprint.

[2] Detroux T, Renson L, Masset L, Kerschen G. 2015 The harmonic balance method for bifurcation analysis of large-scale nonlinear mechanical systems. Computer Methods in Applied Mechanics and Engineering 296, 18–38. doi: [10.1016/j.cma.2015.07.017](https://doi.org/10.1016/j.cma.2015.07.017).

## Disclaimer
This code comes with no warranty.

## License
CC BY 4.0
