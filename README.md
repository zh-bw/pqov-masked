# Masked UOV Implementation

This is the C implementation of masked UOV (Unbalanced Oil and Vinegar), a post-quantum cryptographic signature scheme with masking countermeasures against side-channel attacks.

It is the artifact of the paper "Masked Solving of Linear Equations System and Application to UOV Signatures" (IACR CHES 2026)

The unmasked part is built upon the [UOV reference implementation](https://github.com/pqov/pqov).



### Compilation
```bash
make all
```
### Benchmarking
For correctness test, run
```bash
./sign_test
```
For performance test, run
```bash
python3 run_benchmark.py
```
and the results appear in bench_res.txt