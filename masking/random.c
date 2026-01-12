#include "random.h"
// #include "utils/nistkat/rng.h"
#include "utils/utils_randombytes.h"

#if defined(RNGXOR)
static unsigned x=123456789, y=362436069, z=521288629;
#endif

#ifdef COUNT
uint64_t count_rand = 0;
#endif

uint32_t rand32(){

#if defined(RNGXOR)
  unsigned t;

  x ^= x << 16;
  x ^= x >> 5;
  x ^= x << 1;

  t = x;
  x = y;
  y = z;
  z = t ^ x ^ y; 
  return z;
#else
uint32_t r;
randombytes(&r, 4);
#ifdef COUNT
    count_rand+=4;
#endif
return r;
#endif
}

uint8_t rand8(){
  #if defined(RNGXOR)
  return (uint8_t)(rand32() & 0xFF);
  #else
  #ifdef COUNT
    count_rand++;
  #endif
  uint8_t r;
  randombytes(&r, 1);
  return r;
  #endif
}