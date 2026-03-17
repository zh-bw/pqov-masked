from os import system, popen
import sys

PATH = "./bench_res.txt"
MAX_ORDER = 6
VERBOSE_COMPILE = False
REDIRECT = "> /dev/null" if not VERBOSE_COMPILE else ""

def run_experiment(param, order, mode):
    """
    mode 1: Standard Timing (Default)
    mode 2: Optimized Timing (RNG=XOR)
    mode 3: Randomness Count (COUNT=1)
    """
    cmd_make = f"make clean {REDIRECT} && make bench_test PARAM={param} ORDER={order}"
    
    if mode == 2:
        cmd_make += " RNG=XOR"
        label = "Fast-PRNG-Timing"
    elif mode == 3:
        cmd_make += " COUNT=1"
        label = "Randomness-usage-test"
    else:
        label = "slow-PRNG-Timing"

    system(cmd_make + REDIRECT)
    output = popen("./bench_test").read()
    return label, output

with open(PATH, 'w') as f:
    for PARAM in range(3, 6):
        variant_name = {3: "UOV-I", 4: "UOV-III", 5: "UOV-V"}[PARAM]
        f.write(f"=== Parameter Set: {variant_name} ===\n")
        print(f"\nTesting {variant_name}...")

        for i in range(1, MAX_ORDER + 1):
            print(f"  Order {i}: ", end='', flush=True)
            
            # 1. Performance Run (Standard)
            _, res_std = run_experiment(PARAM, i, mode=1)
            
            # 2. Performance Run (Xorshift)
            # This triggers timing_gadgets() specifically due to the RNGXOR macro
            _, res_xor = run_experiment(PARAM, i, mode=2)
            
            # 3. Randomness usage test Run
            # This triggers test_rand_usage() specifically due to the COUNT macro
            _, res_count = run_experiment(PARAM, i, mode=3)
            
            combined_res = (
                f"--- Order {i} ---\n"
                f"[Results using Slow PRNG]\n {res_std.strip()}\n"
                f"[Results using Fast PRNG]\n {res_xor.strip()}\n"
                f"[Randomness usage]\n {res_count.strip()}\n\n"
            )
            
            f.write(combined_res)
            print("Done.")

print(f"\nAll tests completed. Results saved to {PATH}")