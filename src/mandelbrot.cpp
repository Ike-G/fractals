#include "mandelbrot.h"

HSV::HSV() : h(0), s(0), v(0) {}
HSV::HSV(uint16_t H, double S, double V) : h(H), s(S), v(V) {}
HSV::~HSV() {}

int inMandelbrotQuick(double x0, double y0, std::pair<double, double>*& seq, int& n, int maxIterations) {
    double x, y;
    double x2 = 0;
    double y2 = 0;
    double w = 0;
    int i = 0;
    n = 0;
    while (x2 + y2 <= 4 && i < maxIterations) {
        x = x2-y2 + x0; // x^2-y^2+c0 - Algebraic form of Re((x+yi)^2+c)
        y = w - x2 - y2 + y0; // 2xy + y0 - Algebraic form of Im((x+yi)^2)+c)
        seq[i+1].first = x;
        seq[i+1].second = y;
        x2 = x*x;
        y2 = y*y;
        w = (x+y)*(x+y);
        ++i; 
    }
    if (i == maxIterations)
        return 1;
    n = i+1;
    return 0;
}

int inJuliaQuick(double x, double y, double cx, double cy, std::pair<double, double>*& seq, int& n, int maxIterations) {
    double x2 = x*x;
    double y2 = y*y;
    double w = (x+y)*(x+y);
    int i = 0;
    n = 0;
    while (x2 + y2 <= 3 && i < maxIterations) {
        x = x2-y2 + cx;
        y = w - x2 - y2 + cy;
        seq[i+1].first = x;
        seq[i+1].second = y;
        x2 = x*x;
        y2 = y*y;
        w = (x+y)*(x+y);
        ++i;
    }
    if (i == maxIterations)
        return 1;
    n = i+1;
    return 0;
}

double inMandelbrot(double x0, double y0, int maxIterations) {
    double x = 0;
    double y = 0;
    double xt;
    double log_zn;
    double nu;
    double i = 0;

    // Taking 2**8 as the bailout radius
    while (x*x+y*y <= (1 << 16) && (int)i < maxIterations) {
        xt = x*x - y*y + x0;
        y = 2*x*y + y0;
        x = xt;
        i++;
    }

    if (i < maxIterations) {
        log_zn = log(x*x+y*y)/2; // phi(z) * log(N)/log(2) 
        nu = log(log_zn / log(2)) / log(2); // n-nu(z)
        i++;
        i -= nu; 
    }
    return i;
}

double interpolation(double d, double y1, double y2, double g1, double g2, double t) {
    double h00, h10, h01, h11;
    h00 = 2*pow(t,3)-3*pow(t,2)+1;
    h10 = pow(t,3)-2*pow(t,2)+t;
    h01 = -2*pow(t,3)+3*pow(t,2);
    h11 = pow(t,3)-pow(t,2);
    return y1*h00 + d*g1*h10 + y2*h01 + d*g2*h11;
}

uint8_t* hermiteInterpolation(std::vector<std::pair<double, uint8_t>> d, int sampleCount) {
    double* alpha = new double[d.size()-1];
    double* beta = new double[d.size()-1];
    double* gamma = new double[d.size()];
    double* delta = new double[d.size()-1];
    double* epsilon = new double[d.size()-1];
    double* zeta = new double[d.size()-1];
    double* dx = new double[d.size()-1];
    double* dy = new double[d.size()-1];

    // Consecutive difference - delta_k = (y_{k+1}-y_k)/(x_{k+1}-x_k)
    for (int i = 0; i < d.size()-1; ++i) {
        dy[i] = (double)(d[i+1].second-d[i].second);
        dx[i] = (d[i+1].first-d[i].first);
        delta[i] = dy[i]/dx[i];
    } 

    // Degree-1 coefficients - gamma_k = (delta_{k-1}+delta_k)/2
    gamma[0] = delta[0];
    gamma[d.size()-1] = delta[d.size()-2];
    // double a, b;
    for (int i = 1; i < d.size()-1; ++i) {
        // Enforce flatness for extremum
        if (delta[i-1]*delta[i] <= 0) 
            gamma[i] = 0;
        else 
            gamma[i] = 3*(dx[i-1]+dx[i])/((dx[i-1]+2*dx[i])/delta[i-1] + (2*dx[i-1]+dx[i])/delta[i]);
    }

    // Degree-2 and 3 coefficients
    for (int i = 0; i < d.size()-1; ++i) {
        epsilon[i] = (3*delta[i]-2*gamma[i]-gamma[i+1])/dx[i];
        zeta[i] = (gamma[i]+gamma[i+1]-2*delta[i])/pow(dx[i],2);
    }

    double tau;
    for (int i = 0; i < d.size()-1; ++i) {
        if (delta[i] != 0) {
            alpha[i] = gamma[i]/delta[i];
            beta[i] = gamma[i+1]/delta[i];       
            if (alpha[i] - pow(2*alpha[i]+beta[i]-3, 2)/(3*(alpha[i]+beta[i]-2)) <= 0) {
                std::cout << "Strict monotonicity failed." << std::endl;
                if (alpha[i]+2*beta[i]-3 <= 0 || 2*alpha[i]+beta[i]-3 <= 0) {
                    std::cout << "Weak monotonicity achieved." << std::endl;
                } else {
                    std::cout << "Weak monotonicity failed." << std::endl;
                }
            }        
        }
    }
    delete[] alpha, beta;

    uint8_t* y = new uint8_t[sampleCount];
    double x, xt, diff;
    int j, l, m, h;
    bool found(false);
    for (int i = 0; i < sampleCount; ++i) {
        x = (double)i/(double)sampleCount;
        l = 0;
        h = d.size()-1;
        while (l <= h) {
            m = (h+l)/2;
            xt = d[m].first;
            if (xt < x) {
                l = m+1; 
            } else if (xt > x) {
                h = m-1;
            } else {
                y[i] = d[m].second;
                found = true;
                break;
            }
        }
        if (!found) {
            j = h > 0 ? h : 0;
            diff = x-d[j].first;
            y[i] = static_cast<uint8_t>(d[j].second+gamma[j]*diff+epsilon[j]*pow(diff,2)+zeta[j]*pow(diff,3));
        }
        found = false;
    }
    delete[] gamma, delta, epsilon, zeta, dx, dy;
    return y;
}

void showGradient(Colour* c, int size, const char* filename) {
    Image image(size,256);
    std::cout << "Initialised gradient image" << std::endl;
    for (int x = 0; x < size; ++x) {
        for (int y = 0; y < 256; ++y) {
            image.setColour(c[x], x, y);
        }
        image.setColour(Colour(255,0,0), x, c[x].r);
        image.setColour(Colour(0,255,0), x, c[x].g);
        image.setColour(Colour(0,0,255), x, c[x].b);
    }
    std::cout << "Finished creating gradient image - now exporting." << std::endl;
    image.Export(filename);
}

Colour* generatePalette(std::vector<std::pair<double, Colour>> controlPoints, int maxIterations) {
    Colour* palette = new Colour[maxIterations];
    std::vector<std::pair<double, uint8_t>> R, G, B;
    uint8_t *Rs, *Gs, *Bs;
    controlPoints.push_back(std::pair<double, Colour>(controlPoints[controlPoints.size()-1].first-1.2, controlPoints[controlPoints.size()-1].second));
    for (std::pair<double, Colour>& x : controlPoints) {
        R.push_back(std::pair<double, uint8_t>(x.first, x.second.r));
        G.push_back(std::pair<double, uint8_t>(x.first, x.second.g));
        B.push_back(std::pair<double, uint8_t>(x.first, x.second.b));
    }
    controlPoints.push_back(std::pair<double, Colour>(controlPoints[0].first+1.2, controlPoints[0].second));
    Rs = hermiteInterpolation(R, maxIterations);
    Gs = hermiteInterpolation(G, maxIterations);
    Bs = hermiteInterpolation(B, maxIterations);
    for (int i = 0; i < maxIterations; ++i) {
        palette[i] = Colour(Rs[i], Gs[i], Bs[i]);
        // printf("%d: %d %d %d\n", i, palette[i].r, palette[i].g, palette[i].b);
    }
    delete[] Rs, Gs, Bs;
    return palette;
}

//Colour doubleToColour(double d, double x, double y) {
    // double r, g, b;
    // r = x;
    // g = 1-x;
    // b = y;
    // return Colour(x*d*d, (1-x)*d*d, y*d*d);
//}

Image mandelbrot(uint16_t width, uint16_t height, std::vector<std::pair<double, Colour>> controlPoints, 
                 int maxIterations, double xMin, double xMax, double yMin, double yMax) {
    Image image(width, height);
    Colour* palette = generatePalette(controlPoints, maxIterations+1);
    // showGradient(palette, maxIterations+1);
    Colour c, c1, c2;
    double hr = yMax-yMin;
    double wr = xMax-xMin;
    double h = (double)height;
    double w = (double)width;
    double xs, ys, d;
    for (uint16_t x = 0; x < width; ++x) {
        for (uint16_t y = 0; y < height; ++y) {
            xs = xMin+(double)x*wr/w;
            ys = yMin+(double)y*hr/h;
            double iteration = inMandelbrot(xs, ys, maxIterations);
            if (floor(iteration) == iteration) { 
                image.setColour(palette[(int)iteration], x, y);
            } else {
                c1 = palette[(int)floor(iteration)];
                c2 = palette[(int)floor(iteration)+1];
                d = iteration-floor(iteration);
                image.setColour(Colour((uint8_t)((double)(c2.r-c1.r)*d+(double)c1.r), 
                                (uint8_t)((double)(c2.g-c1.g)*d+(double)c1.g), 
                                (uint8_t)((double)(c2.b-c1.b)*d+(double)c1.b)), x, y); 
            }
        }
    }
    delete[] palette;
    return image;
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

double colourIntensity(double x) {
    return x <= 0.25 ? 3*x : x/3 + 2/3;
}

double spatialColour(double x, double y) {
    return 1.25*(-0.5*(x+y-1)/(1+(x-0.5)*(y-0.5))+0.4);
}

Image buddhabrot(uint16_t width, uint16_t height, 
                 int maxIterations, long nSamples, double xMin, double xMax, double yMin, double yMax) {
    uint16_t* im = new uint16_t[width*height];
    std::mt19937 rng;
    std::uniform_real_distribution<double> xDistribution(xMin, xMax);
    std::uniform_real_distribution<double> yDistribution(yMin, yMax);
    rng.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    int n, xs, ys, biggest(0);
    double x, y;
    for (long i = 0; i < nSamples; ++i) {
        x = xDistribution(rng);
        y = yDistribution(rng);
        std::pair<double, double>* sequence = new std::pair<double, double>[maxIterations+1];
        sequence[0] = std::pair<double, double>(x, y);
        if (!inMandelbrotQuick(x, y, sequence, n)) {
            for (int i = 0; i < n; ++i) {
                if (sequence[i].first < xMax && sequence[i].second < yMax
                    && sequence[i].first >= xMin && sequence[i].second >= yMin) {
                    xs = static_cast<uint16_t>(width * (sequence[i].first - xMin) / (xMax - xMin));
                    ys = static_cast<uint16_t>(height * (sequence[i].second - yMin) / (yMax - yMin));
                    ++im[xs + ys * width];
                    if (im[xs + ys * width] > biggest)
                        biggest = im[xs + ys * width];
                }
            }
        };
        delete[] sequence;
    }
    std::cout << biggest << std::endl;
    int smallest(biggest);
    for (int i = 0; i < width*height; ++i) {
        if (im[i] < smallest)
            smallest = im[i];
    }
    // Colour* palette = generatePalette(controlPoints, maxIterations);
    std::cout << "Finished calculations." << std::endl;
    double v, w((double)width), h((double)height);
    // Colour c1, c2;
    Image image(width, height);
    std::cout << "Initialised image." << std::endl;
    for (uint16_t i = 0; i < width; ++i) {
        for (uint16_t j = 0; j < height; ++j) {
            v = colourIntensity(static_cast<double>((im[i+j*width]-smallest)/(double)(biggest-smallest)));
            // p = spatialColour((double)i/w, (double)j/h);
            image.setColour(Colour((uint8_t)(((double)i/w)*v*255), (uint8_t)((1-((double)i/w))*v*255), (uint8_t)(((double)j/h)*v*255)), i, j);
            // v = static_cast<double>(maxIterations*(im[i+j*width]-smallest))/(double)(biggest-smallest);
            // if (floor(p) == p) { 
            //     image.setColour(palette[static_cast<int>(p)]*v, x, y);
            // } else {
            //    c1 = palette[static_cast<int>(floor(p))]*v;
            //    c2 = palette[static_cast<int>(floor(p)+1)]*v;
            //    d = p-floor(p);
            //    image.setColour(Colour((uint8_t)((double)(c2.r-c1.r)*d+(double)c1.r), 
            //                    (uint8_t)((double)(c2.g-c1.g)*d+(double)c1.g), 
            //                    (uint8_t)((double)(c2.b-c1.b)*d+(double)c1.b)), i, j); 
            //}
        }
    }
    delete[] im;
    std::cout << "Created image." << std::endl;
    return image;
}

Image juliaBuddhabrot(uint16_t width, uint16_t height, double cx, double cy, int maxIterations, 
                      long nSamples, double xMin, double xMax, double yMin, double yMax) {
    uint32_t* im = new uint32_t[width*height];
    std::mt19937 rng;
    std::uniform_real_distribution<double> xDistribution(xMin, xMax);
    std::uniform_real_distribution<double> yDistribution(yMin, yMax);
    rng.seed(std::chrono::high_resolution_clock::now().time_since_epoch().count());
    uint32_t biggest(0);
    int n, xs, ys;
    double x, y;
    // R2 = sqrt(cx*cx+cy*cy)+1;
    for (long i = 0; i < nSamples; ++i) {
        x = xDistribution(rng);
        y = yDistribution(rng);
        std::pair<double, double>* sequence = new std::pair<double, double>[maxIterations+1];
        sequence[0] = std::pair<double, double>(x, y);
        if (!inJuliaQuick(x, y, cx, cy, sequence, n)) {
            for (int i = 0; i < n; ++i) {
                if (sequence[i].first < xMax && sequence[i].second < yMax
                    && sequence[i].first >= xMin && sequence[i].second >= yMin) {
                    xs = static_cast<uint16_t>(width * (sequence[i].first - xMin) / (xMax - xMin));
                    ys = static_cast<uint16_t>(height * (sequence[i].second - yMin) / (yMax - yMin));
                    ++im[xs + ys * width];
                    if (im[xs + ys * width] > biggest)
                        biggest = im[xs + ys * width];
                }
            }
        };
        delete[] sequence;
    }
    std::cout << biggest << std::endl;
    uint32_t smallest(biggest);
    for (int i = 0; i < width*height; ++i) {
        if (im[i] < smallest)
            smallest = im[i];
    }
    std::cout << smallest << std::endl;
    // Colour* palette = generatePalette(controlPoints, maxIterations);
    std::cout << "Finished calculations." << std::endl;
    double v, w((double)width), h((double)height);
    // Colour c1, c2;
    Image image(width, height);
    std::cout << "Initialised image." << std::endl;
    for (uint16_t i = 0; i < width; ++i) {
        for (uint16_t j = 0; j < height; ++j) {
            v = colourIntensity(static_cast<double>((long double)(im[i+j*width]-smallest)/(long double)(biggest-smallest)));
            // p = spatialColour((double)i/w, (double)j/h);
            image.setColour(Colour((uint8_t)(((double)i/w)*v*255), (uint8_t)((1-((double)i/w))*v*255), (uint8_t)(((double)j/h)*v*255)), i, j);
            // v = static_cast<double>(maxIterations*(im[i+j*width]-smallest))/(double)(biggest-smallest);
            // if (floor(p) == p) { 
            //     image.setColour(palette[static_cast<int>(p)]*v, x, y);
            // } else {
            //    c1 = palette[static_cast<int>(floor(p))]*v;
            //    c2 = palette[static_cast<int>(floor(p)+1)]*v;
            //    d = p-floor(p);
            //    image.setColour(Colour((uint8_t)((double)(c2.r-c1.r)*d+(double)c1.r), 
            //                    (uint8_t)((double)(c2.g-c1.g)*d+(double)c1.g), 
            //                    (uint8_t)((double)(c2.b-c1.b)*d+(double)c1.b)), i, j); 
            //}
        }
    }
    delete[] im;
    std::cout << "Created image." << std::endl;
    return image;
}