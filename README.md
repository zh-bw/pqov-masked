# Masked UOV Implementation

This is the C implementation of masked UOV (Unbalanced Oil and Vinegar) signature, and we implemented its pkc variant.

It is the artifact of the [paper](https://eprint.iacr.org/2026/048.pdf) "Masked Solving of Linear Equations System and Application to UOV Signatures" (IACR TCHES 2026)

The unmasked part is built upon the [UOV reference implementation](https://github.com/pqov/pqov).

### Compilation
```bash
make all
```
### Benchmarking

To replicate Table 2, enable #define COUNT in 'src/params.h', and uncomment
``` C
test_rand_usage()
```
in 'test/bench_test.c'.

Then run
```bash
python3 run_benchmark.py
```
The results appear in bench_res.txt

To replicate Table 3 and 4, uncomment
``` C
timing_gadgets();
timing_masked_ov();
```
in 'test/bench_test.c', and run

```bash
python3 run_benchmark.py
```
The results appear in bench_res.txt

The PRNG is the slow one by default, to select the xorshift PRNG (fast one), enable
#define RNGXOR in 'src/params.h'.
