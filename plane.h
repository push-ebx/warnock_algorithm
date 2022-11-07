#include <iostream>
#include <math.h>
#include <limits>
#include <algorithm>
#include <array>
#include "graphics.h"
#include "Matrix.h"

#define deg2rad(deg) (deg * M_PI / 180.0)

#define width 720
#define height 720
#define delay 1
#define angle 7
#define inc_coef 1.1
#define speed 0.05

#define vertex 6

#define UP_KEY 72
#define LEFT_KEY 75
#define DOWN_KEY 80
#define RIGHT_KEY 77
#define MINUS_KEY 45
#define PLUS_KEY 61
#define W_KEY 119
#define S_KEY 115
#define A_KEY 97
#define D_KEY 100
#define Z_KEY 122
#define X_KEY 120
#define ESC_KEY 27

using namespace std;