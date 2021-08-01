#pragma once
#include "image.h"
#include <complex>
#include <math.h>
#include <vector>
#include <utility>
#include <iostream>
#include <stdio.h>
#include <random>
#include <chrono>

struct HSV {
    uint16_t h;
    double s, v;

    HSV();
    HSV(uint16_t h, double s, double v);
    ~HSV();
};

int inMandelbrotQuick(double x, double y, std::pair<double, double>*& seq, int& n, int maxIterations = 1000);
int inJuliaQuick(double x, double y, double cx, double cy, std::pair<double, double>*& seq, int& n, int maxIterations = 1000);
double inMandelbrot(double x, double y, int maxIterations = 1000);
Image mandelbrot(uint16_t width, uint16_t height, std::vector<std::pair<double, Colour>> controlPoints, 
                 int maxIterations=1000, double xMin = -2.1, double xMax = 0.6, double yMin = -1.2, double yMax = 1.2);
Colour doubleToColour(double d, double x, double y);
Colour* generatePalette(std::vector<std::pair<double, Colour>> colourSet, int maxIterations);
uint8_t* hermiteInterpolation(std::vector<std::pair<double,uint8_t>> d, int sampleCount);
double interpolation(double d, double y1, double y2, double g1, double g2, double t);
void showGradient(Colour* c, int size, const char* filename="grad.bmp");
Image buddhabrot(uint16_t width, uint16_t height, 
                 int maxIterations=1000, long nSamples=20000000, double xMin=-2.1, double xMax=0.6, double yMin=-1.2, double yMax=1.2);
Image juliaBuddhabrot(uint16_t width, uint16_t height, double cx, double cy, int maxIterations=1000, long nSamples = 100000000, 
                      double xMin=-1.5, double xMax=1.5, double yMin=-1.2, double yMax=1.2);
