#include "cpucycles.h"
#include <stdint.h>
#include <time.h>

#if defined(__x86_64__) || defined(__i386__)  // Check if compiling on x86

#include <x86intrin.h>

int64_t cpucycles(void) {
    return __rdtsc();  // Use Clang-supported intrinsic
}

#elif defined(__aarch64__)  // Check if compiling on ARM (Apple Silicon)

int64_t cpucycles(void) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (int64_t)ts.tv_sec * 1000000000LL + ts.tv_nsec;  // Return nanoseconds
}

#else
#error "Unsupported architecture"
#endif