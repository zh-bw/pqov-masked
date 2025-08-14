#ifndef RANDOM_H
#define RANDOM_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "params.h"

#ifndef RNG_MODE
#define RNG_MODE 2
#endif


#ifdef COUNT
extern uint64_t count_rand;
#endif

uint8_t  rand8();
uint16_t rand16();
uint32_t rand32();
uint64_t rand64();
#endif
