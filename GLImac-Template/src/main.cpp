#include "declarations.hpp"
#include "GLObject.hpp"




GLObject
circle(int n_points)
{
    std::vector<Vertex> verts;

    for (int i = 0; i < n_points; ++i)
    {
        float angle = ((float) i)/((float) n_points) * M_PI*2.0;
        float x = cos(angle);
        float y = sin(angle);

        Vertex v = {glm::vec2(x, y)};
        
        verts.push_back(v);
    }
    GLObject obj(verts);
    return obj;
}



int main(int argc, char** argv)
{

    // plane from 3  points (plane z=0)
    c2ga::Mvec<double> pt1 = c2ga::point<double>(1,2);
    c2ga::Mvec<double> pt2 = c2ga::point<double>(3,2);
    c2ga::Mvec<double> pt3 = c2ga::point<double>(2,3);
    c2ga::Mvec<double> plane = pt1 ^ pt2 ^ pt3 ^ c2ga::ei<double>();

    std::cout << pt1 << std::endl;
    
        
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


    auto obj = circle(64);


    

    glEnable(GL_DEPTH_TEST);
    

    
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


        // test draw
        obj.draw(GL_LINE_LOOP);

        
        // Update the display
        windowManager.swapBuffers();
    }

    obj.free_gpu();
    
    return EXIT_SUCCESS;
}

