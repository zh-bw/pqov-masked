#ifndef RANDOM_H
#define RANDOM_H
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <time.h>
#include "params.h"

#ifdef COUNT
extern uint64_t count_rand;
#endif

uint8_t  rand8();
uint32_t rand32();
#endif
