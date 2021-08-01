#pragma once
#include <vector>
#include <cstdint>

struct Colour {
    uint8_t r, g, b;

    Colour();
    Colour(uint8_t r, uint8_t g, uint8_t b);
    ~Colour();
    Colour operator*(double x) const;
};

class Image {
public:
    Image(uint16_t width, uint16_t height);
    ~Image();

    Colour getColour(uint16_t x, uint16_t y);
    void setColour(const Colour& colour, uint16_t x, uint16_t y);
    void Export(const char* path);

private:
    uint16_t m_width;
    uint16_t m_height;
    std::vector<Colour> m_colours;
};
