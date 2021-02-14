// Copyright (c) 2018 by University Paris-Est Marne-la-Vallee
// cgaTools.hpp
// Authors: Stephane Breuils and Vincent Nozick
// Contact: vincent.nozick@u-pem.fr

/// \file c3gaTools.hpp
/// \author Stephane Breuils, Vincent Nozick
/// \brief some useful functions when using conformal geometric algebra of R^3. Use this file if you generated the lib using the "c3ga.conf" configuration file.


// Anti-doublon
#ifndef C2GA_TOOLS_HPP__
#define C2GA_TOOLS_HPP__
#pragma once

// Internal Includes
#include <c2ga/Mvec.hpp>


/// \namespace grouping the multivectors object
namespace c2ga{

    /// \brief build a point from a vector
    /// \param x vector component related to e1
    /// \param y vector component related to e2
    /// \param z vector component related to e3
    /// \return a multivector corresponding to a point p = e0 +  x e1 + y e2 + z e3 + 0.5 || (x e1 + y e2 + z e3) ||^2 einf
    template<typename T>
    c2ga::Mvec<T> point(const T &x, const T &y){

        c2ga::Mvec<T> mv;
        mv[c2ga::E1] = x;
        mv[c2ga::E2] = y;
        mv[c2ga::Ei] = 0.5 * mv.quadraticNorm();
        mv[c2ga::E0] = 1.0;

        return mv;
    }

    /// \brief build a point from a vector
    /// \param vec is should be a multivector of grade 1 vec = v1 e1 + v2 e2 + v3 e3. If vec has other components, they will be ignored during the execution.
    /// \return a multivector corresponding to a point p = e0 + v1 e1 + v2 e2 + v3 e3 + 0.5 || vec ||^2 einf
    template<typename T>
    c2ga::Mvec<T> point(const c2ga::Mvec<T> &vec){
        return point(vec[c2ga::E1], vec[c2ga::E2]);
    }

    /// \brief build a dual sphere from a center and a radius
    /// \param centerX dual sphere center component related to e1
    /// \param centerY dual sphere center component related to e2
    /// \param centerZ dual sphere center component related to e3
    /// \param radius of the sphere
    /// \return a multivector corresponding to a dual sphere s = center - 0.5 radius ei, with center being e0 +  x e1 + y e2 + z e3 + 0.5 || (x e1 + y e2 + z e3) ||^2 einf.
    template<typename T>
    c2ga::Mvec<T> dualSphere(const T &centerX, const T &centerY, const T &radius){
        c2ga::Mvec<T> dualSphere = point(centerX,centerY);
        dualSphere[c2ga::Ei] -= 0.5*radius;
        return dualSphere;
    }

    /// \brief extract the center and radius of a dual sphere
    /// \param dualSphere the dual sphere
    /// \param radius positive number for real spheres and negative for imaginary spheres.
    /// \param center center of the dual sphere
    template<typename T>
    void radiusAndCenterFromDualSphere(const c2ga::Mvec<T> &dualSphere, T &radius, c2ga::Mvec<T> &center){
        radius = dualSphere | dualSphere;
        center = dualSphere;
        dualSphere[c2ga::Ei] += 0.5*radius*radius;
    }

    /// \brief interpret the nature of the geometric object (line, circle, pair point, ...)
    /// \param multivector: the multivector to be studied
    /// \todo ... todo :)
    template<typename T>
    void whoAmI(const c2ga::Mvec<T> &mv){

        std::vector<unsigned int> grades_ = mv.grades();
        if(grades_.size() == 0) {
            std::cout << "null vector" << std::endl;
            return;
        }

        if(grades_.size() == 1)
            std::cout << grades_[0] << "-vector (homogeneous)" << std::endl;
        else
            std::cout << "non-homogeous multivector" << std::endl;
    }


    /// \brief extract 2 points pt1 and pt2 from a pair point p = pt1 ^ pt2
    /// \param pairPoint implicitly contains 2 points
    /// \param epsilon is the minimum threshold to specify if 2 points are disjoint
    /// \return a list of 2 points (if they are disjoint) or a single point.
    template<typename T>
    std::vector<c2ga::Mvec<T>> extractPairPoint(const c2ga::Mvec<T> &pairPoint, const T &epsilon = 1.0e-7){

        std::vector<c2ga::Mvec<T>> points;
        T innerSqrt = sqrt(pairPoint | pairPoint);
        if(innerSqrt < epsilon)
            points.push_back(pairPoint / pairPoint[c2ga::E0]);
        else {
            points.push_back( (pairPoint+innerSqrt)/ pairPoint[c2ga::E0]);
            points.push_back( (pairPoint-innerSqrt)/ pairPoint[c2ga::E0]);
        }
        return points;
    }


	/// \brief compute the normal of a surface on a point.
	/// \param surface is a quad vector (sphere of plane).
	/// \param point is a normalized point (e0 = 1) lying on the surface where the normal is estimated.
	/// \return a normal vector (e1,e2,e3) with L2 norm = 1 
    template<typename T>
    c2ga::Mvec<T> surfaceNormal(c2ga::Mvec<T> &surface, c2ga::Mvec<T> &point){
        std::cout << "not sure if this works at all" << std::endl;
        c2ga::Mvec<T> normal;
        normal[c2ga::E1] = - point[c2ga::E1] * surface[c2ga::E012] / point[c2ga::E0] + surface[c2ga::E02i];  
        normal[c2ga::E2] = - point[c2ga::E2] * surface[c2ga::E012] / point[c2ga::E0] - surface[c2ga::E01i]; 
        
        normal = normal / (double) sqrt(normal[c2ga::E1]*normal[c2ga::E1] + normal[c2ga::E2]*normal[c2ga::E2]);
        
        return normal;
    }


} // namespace

#endif // projection_inclusion_guard
