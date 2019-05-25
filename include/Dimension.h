#pragma once

#include <vector>

class Dimension : protected std::vector<unsigned int> {};

class Dimension2D : protected Dimension {
  public:
    unsigned int getRow() const { return this->at(0); }
    unsigned int getColumn() const { return this->at(1); }
};
