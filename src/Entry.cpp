// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Entry.cpp
// This file is part of the Garamon Generator.
// Authors: Vincent Nozick and Stephane Breuils
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program


#include "Entry.hpp"


// Entry::Entry(const std::string &equation, const std::string &objectName, const std::vector<unsigned int> &color)
//    : _equation(equation), _objectName(objectName), _color(color)
//    {
//    }


Entry::Entry(const std::string &equation, const std::string &objectName, const int &red, const int &green, const int &blue)
   : _equation(equation), _objectName(objectName)
   {
      _color = std::vector<int>(3);
      _color[0] = red;
      _color[1] = green;
      _color[2] = blue;
   }



Entry::~Entry()
{}

void Entry::display() const
{
    std::cout << "Entry " << std::endl;
    std::cout << "  object  : " << _objectName << std::endl;
    std::cout << "  equation: " << _equation << std::endl;
    std::cout << "  color   : (";
    for(int i=0; i<(int)_color.size()-1; ++i)
        std::cout << std::to_string(_color[i]) << ", ";
    std::cout << std::to_string(_color[_color.size()-1]) << ")" << std::endl;
    std::cout << std::endl;
}



