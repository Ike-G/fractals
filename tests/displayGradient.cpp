#include "Image.h"
#include <functional>
#include <cmath>
#include <iostream>

struct HSV {
    uint16_t h = 0;
    double s = 0;
    double v = 0;

    HSV(uint16_t H, double S, double V) : h(H), s(S), v(V) {}
};

void showGradient(std::function<Colour(double, double)> f, const char* filename) {
    Image image(1000,1000);
    for (uint16_t x = 0; x < 1000; x++) {
        for (uint16_t y = 0; y < 1000; y++) { 
            image.setColour(f((double)x/1000,(double)y/1000), x, y);
        }
    }
    image.Export(filename);
}

Colour a(double x, double y) {
    return Colour(x*255,(1-x)*255,y*255);
}
Colour hsvToRgb(HSV z) {
    double a, b, c, d, e;
    long i;
    if (z.s <= 0.0)
        return Colour(z.v, z.v, z.v);
    a = z.h % 360;
    a /= 60;
    i = (long)a;
    b = a-i;
    c = z.v * (1. - z.s);
    d = z.v * (1. - (z.s*b));
    e = z.v * (1. - (z.s*(1. - b)));
    switch (i) {
        case 0:
            return Colour(z.v*255, e*255, c*255);
        case 1:
            return Colour(d*255, z.v*255, c*255);
        case 2:
            return Colour(c*255, z.v*255, e*255);
        case 3:
            return Colour(c*255, d*255, z.v*255);
        case 4:
            return Colour(e*255, c*255, z.v*255);
        case 5:
        default:
            return Colour(z.v*255, c*255, d*255);
    }
}

Colour b(double x, double y) {
    return hsvToRgb(HSV(2*360*(int)(x*(1-x)+y*(1-y)), 1., 1.));
}

Colour c(double x, double y) {
    return hsvToRgb(HSV(360*x*y, 1., 1.));
}

Colour d(double x, double y) {
    return hsvToRgb(HSV(180*(pow(x,3)-pow(x,2)*y+pow(y,2)*x-pow(y,3)), 1., 1.));
}

Colour e(double x, double y) {
    return hsvToRgb(HSV(90*(pow(x,3)+pow(x,2)*y+pow(y,2)*x+pow(y,3)+0.5), 1., 1.));
}

Colour f(double x, double y) {
    return hsvToRgb(HSV(95+90*(-(x+y-1)/(1+(2*x-1)*(2*y-1))+1), 1., 1.));
}

Colour g(double x, double y) {
    return hsvToRgb(HSV(120*(pow(2*x-1,3)-3*(2*x-1)*pow(2*y-1,2)+1), 1., 1.));
}

Colour h(double x, double y) {
    return hsvToRgb(HSV(360*(-0.5*(x+y-1)/(1+(x-0.5)*(y-0.5))+0.5), 1., 1.));
}


int main() { 
    showGradient(h, "hsv5.bmp");
    return 0;
}