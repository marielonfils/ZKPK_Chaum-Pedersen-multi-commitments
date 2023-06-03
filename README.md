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
The repository is composed as follows inside the folder _zkpk_python_:
- _utils_ folder :
    - _group.py_ : provides random group exponents, exponentiation with precomputation and proof size computation
    - _constants.py_ : provides group parameters and loading of generators
    - _constants.npz_ : contains 33024 generators
    - _hash2.py_ : provides a hash function
- _bulletproof.py_ : provides an implementation of the **Bulletproofs inner-product** argument
- _extended_schnorr.py_ : provides an implementation of the **Extended Schnorr** argument
- _logarithmic.py_ : provides an implementation of the **Logarithmic 0-1** proof
- _logarithmic_precomputation.py_ : provides an implementation of the **Logarithmic 0-1** proof using precomputed tables for group exponentiations
- _minmax_multibatching.py_ : provides an implementation of the **K-selection** proof
- _minmax_multibatching_precomputation.py_ : provides an implementation of the **K-selection** proof using precomputed tables for group exponentiations
- _partial_opening.py_ : provides an implementation of the **Partial Opening** protocol
- _partial_opening_precomputation.py_ : provides an implementation of the **Partial Opening** protocol using precomputed tables for group exponentiations

## Running the code
For running the code please replace the following values by values of your choice
- _l_ : integer representing the number of choices, should be a power of 2
- _m_ : integer representing the number of voters; should be a power of 2 
- _n_ : integer such that the selection limit is 2^n-1.
- _vs_ : array of size lxm containing the votes as 0/1 values (by default, the votes are chosen randomly)
When running the code, the prover's function will create a proof and the verifier's function will verify it.
