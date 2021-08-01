#include "mandelbrot.h"

#define PHI 1.618033988749895

void scales(Image image, uint16_t width, uint16_t height, std::string fileName) {
    for (int x = width/10.; x < width; x += width/10) {            
        for (int y = 0; y < height; y++) {
            image.setColour(Colour(255,255,255), x, y);
        } 
    }
    for (int y = height/10; y < height; y += height/10) {
        for (int x = 0; x < width; x++) {
            image.setColour(Colour(255,255,255), x, y);
        }
    }
    image.Export(fileName.c_str());
}

void mandelbrotZoom(uint16_t width, uint16_t height, std::vector<std::pair<double, Colour>> controlPoints, int maxIterations,
    double xMin, double xMax, double yMin, double yMax) {
    std::string input, xi, yi, xh, yh;
    double x(xMin), y(yMin), dx((double)xMax-xMin), dy((double)yMax-yMin);
    while (input != "q") {
        Image currentImage = mandelbrot(width, height, controlPoints, maxIterations, x, x+dx, y, y+dy);
        currentImage.Export(("b"+xh+"-"+yh+".bmp").c_str());
        scales(currentImage, width, height, "b"+xh+"-"+yh+"-Scales.bmp");
        dx /= 10;
        dy /= 10;
        std::cout << "\nEnter new x [0-9]: ";
        std::cin >> xi;
        x = dx*static_cast<double>(xi.c_str()[0] - '0')+x;
        std::cout << "\nEnter new y [0-9]: ";
        std::cin  >> yi;
        y = dy*static_cast<double>(yi.c_str()[0] - '0')+y;
        xh += xi[0];
        yh += yi[0];
        std::cout << "\nQuit? ";
        std::cin >> input;
    }
}

void buddhabrotZoom(uint16_t width, uint16_t height, int maxIterations,
    long nSamples, double xMin, double xMax, double yMin, double yMax) {
    std::string input, xi, yi, si, xh, yh;
    double x(xMin), y(yMin), dx((double)xMax-xMin), dy((double)yMax-yMin);
    while (input != "q") {
        Image currentImage = buddhabrot(width, height, maxIterations, nSamples, x, x+dx, y, y+dy);
        currentImage.Export(("b"+xh+"-"+yh+".bmp").c_str());
        scales(currentImage, width, height, "b"+xh+"-"+yh+"-Scales.bmp");
        dx /= 10;
        dy /= 10;
        std::cout << "\nEnter new x [0-9]: ";
        std::cin >> xi;
        x = dx*static_cast<double>(xi.c_str()[0] - '0')+x;
        std::cout << "\nEnter new y [0-9]: ";
        std::cin  >> yi;
        y = dy*static_cast<double>(yi.c_str()[0] - '0')+y;
        std::cout << "\nEnter sample number: ";
        std::cin >> si;
        nSamples = stoi(si);
        xh += xi[0];
        yh += yi[0];
        std::cout << "\nQuit? ";
        std::cin >> input;
    }
}

void juliaIteration(uint16_t width, uint16_t height, double cxMin=-0.8, double cxMax=0.8, 
                    double cyMin=-0.2, double cyMax=0.2, int cxSamples=5, int cySamples=5) {
    double dcx = (cxMax-cxMin)/cxSamples;
    double dcy = (cyMax-cyMin)/cySamples;
    int i(0), j;

    for (double cx = cxMin; cx <= cxMax; cx+=dcx) {
        j = 0;
        for (double cy = cyMin; cy <= cyMax; cy+=dcy) {
            Image current = juliaBuddhabrot(width, height, cx, cy);
            current.Export(("j-"+std::to_string(i)+"-"+std::to_string(j)+".bmp").c_str());
            j++;
        }
        i++;
    }
}

int main() {
    const uint16_t width = 1080;
    const uint16_t height = 1080;

    Colour c1(0,7,100), c2(32,107,203), c3(237,255,255), c4(255,170,0), c5(0,2,0);
    std::vector<std::pair<double, Colour>> blueWhiteOrangeBlack {
        std::pair<double, Colour>(0.0, c1),
        std::pair<double, Colour>(0.16, c2),
        std::pair<double, Colour>(0.42, c3),
        std::pair<double, Colour>(0.6425, c4),
        std::pair<double, Colour>(0.8575, c5),
    };

    std::vector<std::pair<double, Colour>> purpleToGreen { 
        std::pair<double, Colour>(0.0, Colour(116,0,184)),
        std::pair<double, Colour>(0.5, Colour(94,96,206)),
        std::pair<double, Colour>(0.6, Colour(83,144,217)),
        std::pair<double, Colour>(0.7, Colour(72,191,227)),
        std::pair<double, Colour>(0.8, Colour(114,239,221)),
        std::pair<double, Colour>(0.9, Colour(128,255,219))
    };


    std::vector<std::pair<double, Colour>> blueDarkToLight { 
        std::pair<double, Colour>(0.0, Colour(3,4,94)),
        std::pair<double, Colour>(0.3, Colour(2,62,138)),
        std::pair<double, Colour>(0.5, Colour(0,119,182)),
        std::pair<double, Colour>(0.6, Colour(0,180,216)),
        std::pair<double, Colour>(0.8, Colour(72,202,228))
    };

    std::cout << "began program" << std::endl;
    //buddhabrotZoom(width, height, blueDarkToLight, 1000, 5000000, -2.1, 0.6, -1.2, 1.2);
    // blueDarkToLight.push_back(std::pair<double, Colour>(1.0, Colour(0,0,0)));
    // mandelbrotZoom(width, height, blueDarkToLight, 1000, -2.1, 0.6, -1.2, 1.2);
    // juliaIteration(width, height);
    std::pair<double, double> points[] = {
        // std::pair<double, double>(1-PHI, 0),
        std::pair<double, double>(PHI-2, PHI-1),
        std::pair<double, double>(0.285, 0),
        std::pair<double, double>(0.285, 0.01),
        std::pair<double, double>(0.45, 0.1428),
        std::pair<double, double>(-0.70176, -0.3842),
        std::pair<double, double>(-0.835, -0.2321),
        std::pair<double, double>(-0.8, 0.156),
        std::pair<double, double>(-0.7269, 0.1889),
        std::pair<double, double>(0, -0.8)
    };

    // int i = 1;
    // for (std::pair<double, double>& point : points) {
    //     Image current = juliaBuddhabrot(width, height, point.first, point.second);
    //     current.Export(("j2-"+std::to_string(i)+".bmp").c_str());
    //     i++;
    // }
    juliaIteration(width, height);
    // showGradient(generatePalette(purpleToGreen, 1001), 1001, "purpleGrad.bmp");
    // showGradient(generatePalette(blueDarkToLight, 1001), 1001, "blueDarkToLight.bmp");
    // Image b = buddhabrot(width, height, blueDarkToLight, 1000);
    // b.Export("buddhabrotColour.bmp");
    // Image m = mandelbrot(width, height, blueDarkToLight);
    // m.Export("mandelbrotColour3.bmp");

    // Image m = mandelbrot(width, height, controlPoints, 1000);
    // Image m = mandelbrot(width, height, controlPoints, 1000, -0.7500000000001, -0.740000000009, 0.09, 0.11);
    // m.Export("mandelbrotColour2.bmp");
    return 0;
} 
   