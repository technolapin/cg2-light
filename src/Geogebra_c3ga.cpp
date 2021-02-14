// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// Geogebra_c3ga.cpp
// Authors: Vincent Nozick and Stephane Breuils 
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file main.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief Convert garamon data to geogebra for visualization


#include <iostream>

#include "Entry.hpp"

#include "Geogebra_c3ga.hpp"


Viewer_c3ga::Viewer_c3ga() {}

Viewer_c3ga::~Viewer_c3ga() {}

void Viewer_c3ga::display() const
{
    for(const auto& e : entries)
        e.display();
}


void Viewer_c3ga::removeNameDoublons(){

    // for each entry, find if there is multiple occurences of the same name
    auto e = entries.begin();
    while(e != entries.end()){

        // check all the remaining entries
        int count = 2;
        auto e2 = e;
        ++e2;
        while(e2 != entries.end()){
            if((*e)._objectName == (*e2)._objectName){
                (*e2)._objectName += std::to_string(count);
                ++count;
            }
            ++e2;
        }
        ++e;
    }
}



void Viewer_c3ga::render(const std::string &filename) {

    // open the template file
    std::string data = readFile("../data/geogebra_c3ga_template.html");

    // remove name doublons
    removeNameDoublons();

    // build output commands
    std::string commands = "";

    // options (can be removed)
    commands += "   api.evalCommand(\"ShowAxes(false)\");\n";
    commands += "   api.evalCommand(\"ShowGrid(false)\");\n";
     
    // for each entry, extract the command
    for(const auto& e : entries){

        // equation
//        commands += "   api.evalCommand(\"" + e._objectName + " : " +  e._equation + "\");\n";
        commands += "   api.evalCommand(\"";
        commands += (e._objectName == "") ? "" : e._objectName + " : ";
        commands +=  e._equation + "\");\n";

        // color
        if(e._color[0] != -1)
            commands += "   api.evalCommand(\"SetColor(" + e._objectName + ", " 
                                + std::to_string(e._color[0]) + ", " 
                                + std::to_string(e._color[1]) + ", " 
                                + std::to_string(e._color[2]) + ")\");\n"; 
    }

    // substitute command on template file
    substitute(data,"__GARAMON_INSTER_DATA_FLAG__", commands);

    // save file
    writeFile(data, filename);
}

