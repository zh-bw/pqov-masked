# Masked UOV Implementation

This is the C implementation of masked UOV (Unbalanced Oil and Vinegar) signature, and we implemented its pkc variant.

It is the artifact of the [paper](https://eprint.iacr.org/2026/048.pdf) "Masked Solving of Linear Equations System and Application to UOV Signatures" (IACR TCHES 2026)

The unmasked part is built upon the [UOV reference implementation](https://github.com/pqov/pqov).
### Requirements
To successfully run this artifact, you will need the following
- Python3
- GNU Make
- OpenSSL
### Compilation
```bash
make all
```
### Benchmarking
All data presented in the paper is obtained on an Intel Xeon Gold 6128 @ 3.4 GHz machine, compiled using gcc 13.3.

To replicate Tables presented in the paper, run
```bash
python3 run_benchmark.py
```
All results appear in bench_res.txt

