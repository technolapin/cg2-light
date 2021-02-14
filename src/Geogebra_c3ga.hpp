// Copyright (c) 2019 by University Paris-Est Marne-la-Vallee
// Geogebra_c3ga.hpp
// Authors: Vincent Nozick and Stephane Breuils 
// Contact: vincent.nozick@u-pem.fr
//
// Licence MIT
// A a copy of the MIT License is given along with this program

/// \file Geogebra_c3ga.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief define an entry for a convertor from Garamon to Geogebra files


#ifndef GA_GEOGEBRA_C3GA_HPP__
#define GA_GEOGEBRA_C3GA_HPP__


#include <list>
#include <vector>
#include <iostream>
#include <string>
#include <algorithm> 
#include <cstdlib>
#include <limits>

#include <c3ga/Mvec.hpp>

#include "c3gaTools.hpp"
#include "Entry.hpp"
#include "Directory.hpp"



class Viewer_c3ga
{
public:

    // list of entries
    std::list<Entry> entries;


public:

    Viewer_c3ga();

    ~Viewer_c3ga();

    void display() const;

    void removeNameDoublons();

    void render(const std::string &filename);

    template<typename T>
    inline void pushPolygon(std::list<c3ga::Mvec<T> > &mvList, const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){ 

        // init equation
        std::string equation = "Polygon( ";

        // add all points for polygon
        for(auto& e : mvList){

            // homogeneous coordinates
            if(fabs(e[c3ga::E0]) > std::numeric_limits<T>::epsilon())
                e = e / e[c3ga::E0];

            // add point to equation
            equation += "Point({" + std::to_string(e[c3ga::E1]) + "," + std::to_string(e[c3ga::E2]) + "," + std::to_string(e[c3ga::E3]) + "})";
        
            // last element ?
            if(&e != &mvList.back())
                equation += ",";
        }

        // close equation
        equation += ")";

        // push the entry
        entries.push_back(Entry(equation, "poly", red, green, blue));


        // add all points (just points)
        for(auto& e : mvList){
#if 0
            // homogeneous coordinates
            if(fabs(e[c3ga::E0]) > std::numeric_limits<T>::epsilon())
                e = e / e[c3ga::E0];

            // add point to equation
            equation = "Point({" + std::to_string(e[c3ga::E1]) + "," + std::to_string(e[c3ga::E2]) + "," + std::to_string(e[c3ga::E3]) + "})";
        
            // push the entry
            entries.push_back(Entry(equation, "pt", red*0.5, green*0.5, blue*0.5));
#else
            // add point to equation
            push(e, red*0.5, green*0.5, blue*0.5);
#endif
        }
    }   


    template<typename T>
    inline int push(c3ga::Mvec<T> mv, const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){    
        return push(mv,"",red,green,blue);
    }

    template<typename T>
    inline int push(c3ga::Mvec<T> mv, std::string objectName = "", const unsigned int &red = -1, const unsigned int &green = -1, const unsigned int &blue = -1 ){

        // exrtact multivector type (point / sphere / ....)
        std::string mvType = c3ga::whoAmI(mv);
        std::cout << "type = " <<  mvType << std::endl;

        // remove space in the name
        objectName.erase(std::remove(objectName.begin(), objectName.end(), ' '), objectName.end());

        // final equation
        std::string equation;

        // grade 1 ///////////////////////////////////

        // Euclidean vector
        if(mvType == "Euclidean vector") {
            equation = "Vector((" + std::to_string(mv[c3ga::E1]) + "," + std::to_string(mv[c3ga::E2]) + "," + std::to_string(mv[c3ga::E3]) + "))";

            // put a default name
            if(objectName == "")
                objectName = "vec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // points
        if(mvType == "point (dual tangent trivector)"){

            // homogeneous coordinate to 1
            mv /= mv[c3ga::E0];
            equation = " Point({" + std::to_string(mv[c3ga::E1]) + "," + std::to_string(mv[c3ga::E2]) + "," + std::to_string(mv[c3ga::E3]) + "})";

            // put a default name
            if(objectName == "")
                objectName = "pt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;        }

        // dual sphere
        if(mvType == "dual sphere" || mvType == "imaginary dual sphere"){
            // back to correct scale
            mv /= mv[c3ga::E0]; 

            // extract signedradius
            double squaredRadius = mv | mv;
            double radius = sqrt(abs(squaredRadius));

            // extract center
            c3ga::Mvec<double> center = mv;
            center /= center[c3ga::E0];

            // equation
            equation = " Sphere((" + std::to_string(mv[c3ga::E1]) + "," + std::to_string(mv[c3ga::E2]) + "," + std::to_string(mv[c3ga::E3]) + "), " + std::to_string(fabs(radius)) +  ")";

            // put a default name
            if(objectName == "")
                objectName = "sphere";          

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // dual plane
        if(mvType == "dual plane"){
            equation = std::to_string(mv[c3ga::E1]) + " x ";
            if(mv[c3ga::E2] >= 0) equation += " + "; 
            equation += std::to_string(mv[c3ga::E2]) + " y ";
            if(mv[c3ga::E3] >= 0) equation += " + "; 
            equation += std::to_string(mv[c3ga::E3]) + " z ";
            if(-mv[c3ga::Ei] >= 0) equation += " + "; 
            equation += std::to_string(-mv[c3ga::Ei]) + " = 0";

            // put a default name
            if(objectName == "")
                objectName = "plane";           

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }



        // grade 2 ///////////////////////////////////

        // tangent vector (dual tangent bivector)
        if(mvType == "tangent vector (dual tangent bivector)") {

            // position and orientation
            c3ga::Mvec<T> pos,dir;
            c3ga::extractTangentVector(mv,pos,dir);
            equation = " Vector( Point({" + std::to_string(pos[c3ga::E1]) + ","  + std::to_string(pos[c3ga::E2]) + ","  + std::to_string(pos[c3ga::E3]) + "}), Point({"
                                          + std::to_string(pos[c3ga::E1]+dir[c3ga::E1]) + ","  + std::to_string(pos[c3ga::E2]+dir[c3ga::E2]) + ","  + std::to_string(pos[c3ga::E3]+dir[c3ga::E3]) + "}) )";

            // put a default name
            if(objectName == "")
                objectName = "tangentVec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // add tangent bivector
            c3ga::Mvec<T> pt1 = pos - 1.0e-5 * dir;
            c3ga::Mvec<T> pt2 = pos + 1.0e-5 * dir;

            equation = "Cylinder(Point({" + std::to_string(pt1[c3ga::E1]) + ","  + std::to_string(pt1[c3ga::E2]) + ","  + std::to_string(pt1[c3ga::E3]) + "}), Point({"
                        + std::to_string(pt2[c3ga::E1]) + ","  + std::to_string(pt2[c3ga::E2]) + ","  + std::to_string(pt2[c3ga::E3]) + "}),0.3)";
            objectName = "tanBiv";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // pair point (imaginary dual circle)
        if(mvType == "pair point (imaginary dual circle)" || mvType == "imaginary pair point (dual circle)") {

            // extract the 2 points
            c3ga::Mvec<T> pt1,pt2;
            c3ga::extractPairPoint(mv,pt1,pt2);

            // equation of point 1
            equation = " Point({" + std::to_string(pt1[c3ga::E1]) + "," + std::to_string(pt1[c3ga::E2]) + "," + std::to_string(pt1[c3ga::E3]) + "})";

            // put a default name
            if(objectName == "")
                objectName = "pp";

            // push the first point
            entries.push_back(Entry(equation, objectName + std::to_string(1), red, green, blue));

            // equation of point 2
            equation = "Point({" + std::to_string(pt2[c3ga::E1]) + "," + std::to_string(pt2[c3ga::E2]) + "," + std::to_string(pt2[c3ga::E3]) + "})";

            // push the entry
            entries.push_back(Entry(equation, objectName + std::to_string(2), red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // flat point
        if(mvType == "flat point") {

            // extract the point
            c3ga::Mvec<T> pt;
            c3ga::extractFlatPoint(mv,pt);

            // equation of point
            equation = " Point({" + std::to_string(pt[c3ga::E1]) + "," + std::to_string(pt[c3ga::E2]) + "," + std::to_string(pt[c3ga::E3]) + "})";

            // put a default name
            if(objectName == "")
                objectName = "flatPt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // dual line
        if(mvType == "dual line"){

            // extract a point and a direction
            c3ga::Mvec<double> direction;
            direction[c3ga::E1] =  mv[c3ga::E23];
            direction[c3ga::E2] = -mv[c3ga::E13];
            direction[c3ga::E3] =  mv[c3ga::E12];
            c3ga::Mvec<double> point;
            point[c3ga::E1] =  mv[c3ga::E1i];
            point[c3ga::E2] =  mv[c3ga::E2i];
            point[c3ga::E3] =  mv[c3ga::E3i];
            point = (direction ^ point) / (direction | direction); // can prabably find better than this

            equation = " Line((" + std::to_string(point[c3ga::E23]) + "," + std::to_string(-point[c3ga::E13]) + "," + std::to_string(point[c3ga::E12]) + "), Vector(("
                                 + std::to_string(direction[c3ga::E1]) + "," + std::to_string(direction[c3ga::E2]) + "," + std::to_string(direction[c3ga::E3]) + ")))";

            // put a default name
            if(objectName == "")
                objectName = "l";           

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }





        // grade 3 ///////////////////////////////////

        // tangent vector (dual tangent bivector)
        if(mvType == "tangent bivector (dual tangent vector)") {

            // position and orientation
            c3ga::Mvec<T> pos,dir;
            c3ga::extractTangentBivector(mv,pos,dir);
            equation = " Vector( Point({" + std::to_string(pos[c3ga::E1]) + ","  + std::to_string(pos[c3ga::E2]) + ","  + std::to_string(pos[c3ga::E3]) + "}), Point({"
                                          + std::to_string(pos[c3ga::E1]+dir[c3ga::E1]) + ","  + std::to_string(pos[c3ga::E2]+dir[c3ga::E2]) + ","  + std::to_string(pos[c3ga::E3]+dir[c3ga::E3]) + "}) )";

            // put a default name
            if(objectName == "")
                objectName = "tangentVec";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // add tangent bivector
            c3ga::Mvec<T> pt1 = pos - 1.0e-5 * dir;
            c3ga::Mvec<T> pt2 = pos + 1.0e-5 * dir;

            equation = "Cylinder(Point({" + std::to_string(pt1[c3ga::E1]) + ","  + std::to_string(pt1[c3ga::E2]) + ","  + std::to_string(pt1[c3ga::E3]) + "}), Point({"
                        + std::to_string(pt2[c3ga::E1]) + ","  + std::to_string(pt2[c3ga::E2]) + ","  + std::to_string(pt2[c3ga::E3]) + "}),0.3)";
            objectName = "tanBiv";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // circle
        if(mvType == "circle (imaginary dual pair point)" || mvType == "imaginary circle (dual pair point)") {

            T radius;
            c3ga::Mvec<T> center, direction;
            extractDualCircle(mv.dual(), radius,  center, direction);
            equation = "Circle( Point({" + std::to_string(center[c3ga::E1]) + ","  + std::to_string(center[c3ga::E2]) + ","  + std::to_string(center[c3ga::E3]) + "}), "
                        + std::to_string(fabs(radius))
                        + ", Vector((" + std::to_string(direction[c3ga::E1]) + ","  + std::to_string(direction[c3ga::E2]) + ","  + std::to_string(direction[c3ga::E3]) + ")))";

            // put a default name
            if(objectName == "")
                objectName = "circle";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }


        // line
        if(mvType == "line"){
 
            // dualize
            mv = mv.dual();

            // extract a point and a direction
            c3ga::Mvec<double> direction;
            direction[c3ga::E1] =  mv[c3ga::E23];
            direction[c3ga::E2] = -mv[c3ga::E13];
            direction[c3ga::E3] =  mv[c3ga::E12];
            c3ga::Mvec<double> point;
            point[c3ga::E1] =  mv[c3ga::E1i];
            point[c3ga::E2] =  mv[c3ga::E2i];
            point[c3ga::E3] =  mv[c3ga::E3i];
            point = (direction ^ point) / (direction | direction);

            equation = " Line((" + std::to_string(point[c3ga::E23]) + "," + std::to_string(-point[c3ga::E13]) + "," + std::to_string(point[c3ga::E12]) + "), Vector(("
                                 + std::to_string(direction[c3ga::E1]) + "," + std::to_string(direction[c3ga::E2]) + "," + std::to_string(direction[c3ga::E3]) + ")))";

            // put a default name
            if(objectName == "")
                objectName = "l";           

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // dual flat point
        if(mvType == "dual flat point") {

            // extract the point
            c3ga::Mvec<T> pt;
            c3ga::extractFlatPoint(mv.dual(),pt);

            // equation of point
            equation = " Point({" + std::to_string(pt[c3ga::E1]) + "," + std::to_string(pt[c3ga::E2]) + "," + std::to_string(pt[c3ga::E3]) + "})";

            // put a default name
            if(objectName == "")
                objectName = "dualFlatPt";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;        }




        // grade 4 ///////////////////////////////////

        // plane
        if(mvType == "plane"){
            equation = std::to_string(-mv[c3ga::E023i]) + " x ";
            if(mv[c3ga::E013i] >= 0) equation += " + "; 
            equation += std::to_string(mv[c3ga::E013i]) + " y ";
            if(-mv[c3ga::E012i] >= 0) equation += " + "; 
            equation += std::to_string(-mv[c3ga::E012i]) + " z ";
            if(mv[c3ga::E123i] >= 0) equation += " + "; 
            equation += std::to_string(mv[c3ga::E123i]) + " = 0";

            // put a default name
            if(objectName == "")
                objectName = "plane";

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;           
        }

        // sphere
        if(mvType == "sphere" || mvType == "imaginary sphere"){

            // dualize
            mv = mv.dual();

            // extract parameters
            mv /= mv[c3ga::E0]; // back to correct scale
            double squareRadius = mv | mv;
            double radius = sqrt(squareRadius);
            c3ga::Mvec<double> center = mv;
            center /= center[c3ga::E0];
            equation = " Sphere((" + std::to_string(mv[c3ga::E1]) + "," + std::to_string(mv[c3ga::E2]) + "," + std::to_string(mv[c3ga::E3]) + "), " + std::to_string(fabs(radius)) +  ")";

            // put a default name
            if(objectName == "")
                objectName = "sphere";          

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        // tangent trivector
        if(mvType == "tangent trivector (dual point)"){

            // dualize
            mv = mv.dual();

            // extract parameters
            mv /= mv[c3ga::E0]; // back to correct scale
            double radius = 0.5;
            equation = " Sphere((" + std::to_string(mv[c3ga::E1]) + "," + std::to_string(mv[c3ga::E2]) + "," + std::to_string(mv[c3ga::E3]) + "), " + std::to_string(fabs(radius)) +  ")";

            // put a default name
            if(objectName == "")
                objectName = "tanTrivec";          

            // push the entry
            entries.push_back(Entry(equation, objectName, red, green, blue));

            // finish
            return EXIT_SUCCESS;
        }

        return EXIT_FAILURE;
    }

};




#endif // GA_GEOGEBRA_C3GA_HPP__
