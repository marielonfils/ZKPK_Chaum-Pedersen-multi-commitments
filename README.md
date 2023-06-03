# ZKPK_Chaum-Pedersen-multi-commitments

### Supervisors: Olivier PEREIRA, Thomas PETERS
### Authors: Marie LONFILS, Alexis VUILLE
### UCL - Academic year 2022–2023


This repository contains the implementation of the protocols developed for the Master thesis "*_ Risk-Limiting Audit Optimization with ElectionGuard - Zero Knowledge Arguments on Chaum-Pedersen multi-commitments_*.
## Requirements
The implementation requires 
- python3.8+
- numpy
- gmpy2

## Structure
The repository is composed as follows :
- _utils_ folder:
    - group.py : provides random group exponents, exponentiation with precomputation and proof size computation
    - constants.py : provides group parameters and loading of generators
    - generators.gmz : contains 33024 generators
    - hash2.py: provides a hash function
- logarithmic.py : provides an implementation of the **Logarithmic 0-1** proof
- logarithmic_precomputation.py : provides an implementation of the **Logarithmic 0-1** proof using precomputed tables for group exponentiations
- minmax_multibatching.py: provides an implementation of the **K-selection** proof
- minmax_multibatching_precomputation.py: provides an implementation of the **K-selection** proof using precomputed tables for group exponentiations
- partial_opening.py: provides an implementation of the **Partial Opening** protocol
- partial_opening_precomputation.py: provides an implementation of the **Partial Opening** protocol using precomputed tables for group exponentiations

## Running the code
For running the code please replace the following values by values of your choice
- l : integer representing the number of choices, should be a power of 2
- m : integer representing the number of voters; should be a power of 2 
- n : integer such that the selection limit is 2^n-1.
- vs: array of size lxm containing the votes as 0/1 values (by default, the votes are chosen randomly)
When running the code, the prover's function will create a proof and the verifier's function will verify it.
