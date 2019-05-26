#pragma once

#include <vector>

class Dimension : protected std::vector<unsigned int> {
  public:
    Dimension() : std::vector<unsigned int>() {}
    Dimension(std::initializer_list<unsigned int> value) : std::vector<unsigned int>(value) {}
};

class Dimension2D : protected Dimension {
  public:
    unsigned int getRow() const { return this->at(0); }
    unsigned int getColumn() const { return this->at(1); }
    Dimension2D() : Dimension() {}
    Dimension2D(std::initializer_list<unsigned int> value) : Dimension(value) {}
};
