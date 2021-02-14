#include <iostream>
#include <cstdlib>
#include <list>
#include <c3ga/Mvec.hpp>

#include "c3gaTools.hpp"
#include "Geogebra_c3ga.hpp"

#include "Entry.hpp"
#include <glm/glm.hpp>
#include <GL/glew.h>
#include <glimac/SDLWindowManager.hpp>
#include <glimac/Program.hpp>

using namespace glimac;

///////////////////////////////////////////////////////////
void sample1()
{
    glm::vec3 v(0.);
    
    Viewer_c3ga viewer;

    auto exercice = 6;
    
    if (exercice == 0)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(1,2,0.5);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        c3ga::Mvec<double> pt4 = c3ga::point<double>(2,2,2.0);
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);
        viewer.push(pt4, "pt4", 200,0,0);

        // a line
        c3ga::Mvec<double> line = pt3 ^ pt2 ^ c3ga::ei<double>();
        viewer.push(line, "line", 0,200,0);

        // a plane
        c3ga::Mvec<double> plane = pt1 ^ pt2 ^ pt3 ^ c3ga::ei<double>();
        viewer.push(plane, "plane", 150,200,0);

        // a sphere
        c3ga::Mvec<double> sphere = pt1 ^ pt2 ^ pt3 ^ pt4;
        viewer.push(sphere, "sphere", 0,0,200);
        std::cout << "sphere " << sphere << std::endl;

        // a circle (sphere-plane intersection)
        c3ga::Mvec<double> circle = !( !plane ^ !sphere );
        // c3ga::Mvec<double> circle = pt1 ^ pt2 ^ pt3;
        viewer.push(circle, "circle", 0,200,0);
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 1)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(1,2,0.5);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        c3ga::Mvec<double> pt4 = c3ga::point<double>(2,2,2.0);
        
        c3ga::Mvec<double> pt5 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt6 = c3ga::point<double>(5,2,0.5);
        c3ga::Mvec<double> pt7 = c3ga::point<double>(4,3,0.5);
        c3ga::Mvec<double> pt8 = c3ga::point<double>(4,2,2.0);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);
        viewer.push(pt4, "pt4", 200,0,0);
        viewer.push(pt5, "pt5", 200,0,0);
        viewer.push(pt6, "pt6", 200,0,0);
        viewer.push(pt7, "pt7", 200,0,0);
        viewer.push(pt8, "pt8", 200,0,0);

        // a sphere
        c3ga::Mvec<double> sphere1 = pt1 ^ pt2 ^ pt3 ^ pt4;
        viewer.push(sphere1, "sphere1", 0,0,200);
        std::cout << "sphere1 " << sphere1 << std::endl;
        c3ga::Mvec<double> sphere2 = pt5 ^ pt6 ^ pt7 ^ pt8;
        viewer.push(sphere2, "sphere2", 0,0,200);
        std::cout << "sphere2 " << sphere2 << std::endl;

        c3ga::Mvec<double> intersection = !(!sphere1 ^ !sphere2);
        viewer.push(intersection, "intersection", 0,200,200);
        std::cout << "intersection " << intersection << std::endl;

        
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 2)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(0,0,0);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);

        c3ga::Mvec<double> cercle = pt1 ^ pt2 ^ pt3;
        viewer.push(cercle, "cercle", 0,0,200);
        std::cout << "cercle " << cercle << std::endl;

        std::vector<c3ga::Mvec<double>> lines;
        for (auto i = 0; i < 10; ++i)
        {
            c3ga::Mvec<double> line = c3ga::randomPoint<double>() ^ c3ga::randomPoint<double>() ^ c3ga::ei<double>();

            auto dual_sph = (!line) | cercle;
            double radius;
            c3ga::Mvec<double> center;
            //auto radius = (dual_sph |dual_sph)/dual_sph[c3ga::E0];
            radiusAndCenterFromDualSphere(dual_sph, radius, center);
            radius = dual_sph | dual_sph;
            if (radius >= 0)
            {
                viewer.push(line, "line", 0,200,0);
            }
            else
            {
                viewer.push(line, "line", 200,0,0);
            }
            lines.push_back(line);
        }
        
        
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 3)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(0,0,0);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);

        c3ga::Mvec<double> cercle = pt1 ^ pt2 ^ pt3;
        viewer.push(cercle, "cercle", 0,0,200);
        std::cout << "cercle " << cercle << std::endl;

        const double t = 0.5; // translation amplitude
        c3ga::Mvec<double> translator = 1 - 0.5*t*c3ga::e1<double>()*c3ga::ei<double>();
        for (auto i = 0; i < 10; ++i)
        {
            cercle = translator * cercle * translator.inv();
            viewer.push(cercle, "cercle", 0,0,200);
            std::cout << "cercle " << cercle << std::endl;
        }
        
        
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 4)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(0,0,0);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);

        c3ga::Mvec<double> cercle = pt1 ^ pt2 ^ pt3;
        viewer.push(cercle, "cercle", 0,0,200);
        std::cout << "cercle " << cercle << std::endl;
        auto porth1 = c3ga::e1<double>();
        auto porth2 = c3ga::e3<double>();
        
        auto axe = porth1 ^ porth2;
        axe /= axe.norm();
        
        auto alpha = 0.1;
        auto rotor = cos(0.5*alpha) - axe*sin(0.5*alpha);
        
        for (auto i = 0; i < 10; ++i)
        {
            cercle = rotor * cercle * rotor.inv();
            viewer.push(cercle, "cercle", 0,0,200);
            std::cout << "cercle " << cercle << std::endl;
        }
        
        
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 5)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(0,0,0);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);

        c3ga::Mvec<double> cercle = pt1 ^ pt2 ^ pt3;
        viewer.push(cercle, "cercle", 0,0,200);
        std::cout << "cercle " << cercle << std::endl;

        auto scale = 1.2;
        auto dilator = 1 - (1-scale)/(1+scale)*c3ga::e0i<double>();
        
        for (auto i = 0; i < 10; ++i)
        {
            cercle = dilator * cercle * dilator.inv();
            viewer.push(cercle, "cercle", 0,0,200);
            std::cout << "cercle " << cercle << std::endl;
        }
        
        
        viewer.display();
        viewer.render("output.html");
    }
    else if (exercice == 6)
    {
        
        // 4  points
        c3ga::Mvec<double> pt1 = c3ga::point<double>(0,0,0);
        c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.5);
        c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.5);
        
        viewer.push(pt1, "pt1", 200,0,0);
        viewer.push(pt2, "pt2", 200,0,0);
        viewer.push(pt3, "pt3", 200,0,0);

        c3ga::Mvec<double> cercle = pt1 ^ pt2 ^ pt3;
        viewer.push(cercle, "cercle", 0,0,200);
        std::cout << "cercle " << cercle << std::endl;

        auto nbIter = 50;
        auto nbTours = 5.;
        auto size = 5.;

        
        auto axe = c3ga::e12<double>();;
        auto alpha = 2.0*3.14*nbTours/( (double) nbIter );
        auto rotor = cos(0.5*alpha) - axe*sin(0.5*alpha);

        const double t = c3ga::e3<double>()*size/((double) nbIter);
        c3ga::Mvec<double> translator = 1 - 0.5*t*c3ga::e1<double>()*c3ga::ei<double>();

        auto scale = pow(10., -1./((double) nbIter));
        auto dilator = 1 - (1-scale)/(1+scale)*c3ga::e0i<double>();

        auto transform = scale*rotor*translator;
        
        for (auto i = 0; i < nbIter; ++i)
        {
            cercle = transform * cercle * transform.inv();
            viewer.push(cercle, "cercle", 0,0,200);
            std::cout << "cercle " << cercle << std::endl;
        }
        
        
        viewer.display();
        viewer.render("output.html");
    }

}


///////////////////////////////////////////////////////////
void sample2(){

    Viewer_c3ga viewer;

    // plane from 3  points (plane z=0)
    c3ga::Mvec<double> pt1 = c3ga::point<double>(1,2,0.0);
    c3ga::Mvec<double> pt2 = c3ga::point<double>(3,2,0.0);
    c3ga::Mvec<double> pt3 = c3ga::point<double>(2,3,0.0);
    c3ga::Mvec<double> plane = pt1 ^ pt2 ^ pt3 ^ c3ga::ei<double>();
    //viewer.push(plane, 50,50,50);

    // sphere
    c3ga::Mvec<double> sphere = c3ga::point<double>(0,0,2.0) 
                              ^ c3ga::point<double>(2,0,2.0) 
                              ^ c3ga::point<double>(1,1,2.0) 
                              ^ c3ga::point<double>(1,0,3.0);
    viewer.push(sphere, 200,200,00);

    // circle
    c3ga::Mvec<double> circle = c3ga::point<double>(0,-0.75,0.5) 
                              ^ c3ga::point<double>(2,-0.75,0.5) 
                              ^ c3ga::point<double>(1,-0.75,1.5);
    viewer.push(circle, 200, 0, 0);

    // line 
    c3ga::Mvec<double> line = c3ga::point<double>(0,-0.75,2.3) 
                              ^ c3ga::point<double>(2,-0.75,2.3) 
                              ^ c3ga::ei<double>();
    viewer.push(line, 200, 0, 0);

    // circle plane intersection
    c3ga::Mvec<double> pp1 = !circle ^ !plane;
    viewer.push(!pp1, 0, 200, 0);

    // circle sphere intersection
    c3ga::Mvec<double> pp2 = !circle ^ !sphere;
    viewer.push(!pp2, 0, 200, 0);

    // line sphere intersection
    c3ga::Mvec<double> pp3 = !line ^ !sphere;
    viewer.push(!pp3, 0, 200, 0);


    viewer.display();
    viewer.render("output.html"); 
}


///////////////////////////////////////////////////////////
void sample_polygons()
{
    Viewer_c3ga viewer;

    // a set of points
    c3ga::Mvec<double> pt1 = c3ga::point<double>(1  , 2, 0.5);
    c3ga::Mvec<double> pt2 = c3ga::point<double>(2  , 3, 0.5);
    c3ga::Mvec<double> pt3 = c3ga::point<double>(3  , 3, 0.5);
    c3ga::Mvec<double> pt4 = c3ga::point<double>(4  , 2, 0.5);
    c3ga::Mvec<double> pt5 = c3ga::point<double>(2.5,-2, 0.5);

    // put points on a list
    std::list<c3ga::Mvec<double>> myList;
    myList.push_back(pt1);
    myList.push_back(pt2);
    myList.push_back(pt3);
    myList.push_back(pt4);
    myList.push_back(pt5);

    viewer.pushPolygon(myList, 0,200,0);

    viewer.display();
    viewer.render("output.html"); 
}


// just a sphere for now
struct Object
{
    c3ga::Mvec<double> point;
};


struct World
{
    std::vector<Object> objects;
};
    





int main(int argc, char** argv)
{
        
    int width = 1000;
    int height = 1000;
    // Initialize SDL and open a window
    SDLWindowManager windowManager(width, height, "GLImac");
    glm::vec2 dims((float) width, (float) height);
    // Initialize glew for OpenGL3+ support
    GLenum glewInitError = glewInit();
    if(GLEW_OK != glewInitError) {
        std::cerr << glewGetErrorString(glewInitError) << std::endl;
        return EXIT_FAILURE;
    }

    std::cout << "OpenGL Version : " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLEW Version : " << glewGetString(GLEW_VERSION) << std::endl;

    /*********************************
     * HERE SHOULD COME THE INITIALIZATION CODE
     *********************************/
    
    FilePath applicationPath(argv[0]);


    // Application loop:
    bool done = false;
    while(!done) {
        // Event loop:
        SDL_Event e;
        while(windowManager.pollEvent(e)) {
            if(e.type == SDL_QUIT) {
                done = true; // Leave the loop after this iteration
            }
        }

        if (windowManager.isKeyPressed(SDLK_ESCAPE))
        {
           done = true;
        }

        // *********************************
        // * HERE SHOULD COME THE RENDERING CODE
        // *********************************


        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        
        // Update the display
        windowManager.swapBuffers();
    }
      
    return EXIT_SUCCESS;
}
