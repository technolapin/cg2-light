#include "declarations.hpp"




void
build_rays(Scene & scene,
           const glm::vec2 & source_coords,
           const glm::vec3 & source_color,
           const std::vector<OpticObject> & objects,
           const float pas)
{
    auto pt_source = CGA::point(source_coords);
    std::vector<BinTree<Ray>> rays;
    for (auto alpha = 0.; alpha < 2.0*M_PI; alpha+=pas)
    {
        auto trans = CGA::translator(CGA::vector(0.1, 0.0));
        auto motor = CGA::motor(alpha, pt_source);
        auto versor = motor*trans;
        auto point2 = CGA::apply_versor(versor, pt_source);

        auto pos = CGA::point_to_glm(point2);

        auto dir = CGA::unit_vector(alpha);

        auto ray = Ray(point2, dir, {1.0, 1.0, 0.1});

        BinTree tree(ray);
        rays.push_back(tree);
        
    } 

    for (auto & tree: rays)
    {
        std::vector<int> stack = {0};

        while (stack.size() != 0)
        {
            int mother = stack.back();
            stack.pop_back();
            const auto & ray = tree.nodes[mother];

            // loop over segments (one for now)
            for (const auto & obj: objects)
            {
                const auto intern = ray.intersect(obj);
                if (intern)
                {
                    const auto & inter = intern.value()[0];
                    const auto & norm = intern.value()[1];
                    const auto refs = ray.reflect_on(obj, inter, norm);
                    const auto refl = refs[0];
                    const auto refr = refs[1];
 
                    if (refl)
                    {
                        int left = tree.push_left(refl.value(), mother); 
                        if (refl.value().power2() >= 0.1)
                        {
                            stack.push_back(left);
                        }
                    }
                    
                    if (refr)
                    {
                        int right = tree.push_right(refr.value(), mother);
                        if (refr.value().power2() >= 0.1)
                        {
                            stack.push_back(right);
                        }
                    }
                }
            }
        }
    }

    for (auto & tree: rays)
    {
        for (auto parent = 0; parent < tree.nodes.size(); ++parent)
        {
            auto & ray = tree.nodes[parent];
            auto & pt0 = ray.source;
            auto & color = ray.color;

            auto left = tree.child_left(parent); 
            auto right = tree.child_right(parent);

            auto ptl = (left >= 0)?tree.nodes[left].source : ray.infinity();

            auto seg_refl = pt0^ptl;

            scene.push(CGA::segment_object(seg_refl), Instance().colored(color));

        }
    }

    
    
}


void
ray_tracer3(Scene & scene, const glm::vec2 & source_coords)
{
    auto pas = 0.1;
    auto cir = disk(10);
    auto cir_empty = circle(100);
    Instance inst;

    auto segment = Circle({0., 0}, 0.2, {0.7, 0.9, 0.9}, 2.0);
    auto outer_circle = CGA::circle(CGA::e0(), 2.0);
    scene.push(cir_empty, Instance()
               .translate(CGA::point_to_glm(segment.center))
               .scale(segment.radius)
               .colored(segment.color));
    
    
    auto pt_source = CGA::point(source_coords);
    {
        auto pos =  CGA::point_to_glm(pt_source);
        Instance inst;
        auto inst2 = inst.translate(pos).colored(1., 1., 0.).scale(0.01);
        scene.push(cir, inst2);
    }
    
    std::vector<BinTree<Ray>> rays;

    for (auto alpha = 0.; alpha < 2.0*M_PI; alpha+=pas)
    {
        auto trans = CGA::translator(CGA::vector(0.1, 0.0));
        auto motor = CGA::motor(alpha, pt_source);
        auto versor = motor*trans;
        auto point2 = CGA::apply_versor(versor, pt_source);

        auto pos = CGA::point_to_glm(point2);

        scene.push(cir, inst.translate(pos)
                   .colored(1., 1., 1.)
                   .scale(0.002));

        auto dir = CGA::unit_vector(alpha);

        auto ray = Ray(point2, dir, {1.0, 1.0, 0.1});

        BinTree tree(ray);
        rays.push_back(tree);
        
    } 

    for (auto & tree: rays)
    {
        std::vector<int> stack = {0};

        while (stack.size() != 0)
        {
            int mother = stack.back();
            stack.pop_back();
            const auto & ray = tree.nodes[mother];

            // loop over segments (one for now)
            {
                const auto intern = ray.intersect(segment);
                if (intern)
                {
                    const auto & inter = intern.value()[0];
                    const auto & norm = intern.value()[1];
                    const auto refs = ray.reflect_on(segment, inter, norm);
                    const auto refl = refs[0];
                    const auto refr = refs[1];
 
                    if (refl)
                    {
                        int left = tree.push_left(refl.value(), mother); 
                        if (refl.value().power2() >= 0.1)
                        {
                            stack.push_back(left);
                        }
                    }
                    
                    if (refr)
                    {
                        int right = tree.push_right(refr.value(), mother);
                        if (refr.value().power2() >= 0.1)
                        {
                            stack.push_back(right);
                        }
                    }
                }
            }
        }
    }



    for (auto & tree: rays)
    {
        for (auto parent = 0; parent < tree.nodes.size(); ++parent)
        {
            auto & ray = tree.nodes[parent];
            auto & pt0 = ray.source;
            auto & color = ray.color;

            auto left = tree.child_left(parent); 
            auto right = tree.child_right(parent);

            auto ptl = (left >= 0)?tree.nodes[left].source : ray.infinity();

            auto seg_refl = pt0^ptl;

            scene.push(CGA::segment_object(seg_refl), inst.colored(color));

            
            //auto ptr = (left >= 0)?tree.nodes[left].source : ray.infinity();
            //auto seg_refr = pt0^ptr;  
            //scene.push(CGA::segment_object(seg_refr), inst.colored(color));
         
        }
    }

    
}




int main(int argc, char** argv)
{
    

    const int width = 1000;
    const int height = 1000;
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

    auto path_vs = applicationPath.dirPath() + "shaders" + "2D_instancied.vs.glsl";
    auto path_fs = applicationPath.dirPath() + "shaders" + "2D.fs.glsl";

    std::cout << path_vs << std::endl;
    Renderer rdr(
        path_vs,
        path_fs,
        GL_LINE_LOOP
        );


    auto pgrm = loadProgram(path_vs, path_fs);

    pgrm.use();



    glLineWidth(3.0);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_DST_ALPHA);  
    glBlendEquation(GL_MAX);
    
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

        const auto mouse_pixel_pos = windowManager.getMousePosition();
        const glm::vec2 mouse_pos = {((float) mouse_pixel_pos[0])/((float) width)*2.0f -1.0f,
            -((float) mouse_pixel_pos[1])/((float) height)*2.0f +1.0f};
        
        Scene scene;

        ray_tracer3(scene, mouse_pos);

        
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // test draw
        //rdr.setup();
        rdr.draw(scene);
        
        // Update the display
        windowManager.swapBuffers();
    }

    
    return EXIT_SUCCESS;
}

