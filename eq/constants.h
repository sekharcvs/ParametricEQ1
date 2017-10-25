#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "options.h"

#define FRAME_SIZE_MS 20 /* Processing frame size in ms */
#define OLA_SIZE_MS 10 /* Overlap add size in ms for DFT */
#define WIN_SIZE_MS FRAME_SIZE_MS + OLA_SIZE_MS

#define N_BANDS 14

#define MS2NSAMPLES(fs, ms) (short)( (((long)fs) * ms) / 1000L)
