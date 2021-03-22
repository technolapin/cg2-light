#include "declarations.hpp"
#include "GLObject.hpp"

using Mvec = c2ga::Mvec<float>;


const float EPSILON = 0.000001;

struct Instance
{
    glm::mat4 transform = glm::mat4(1);
    glm::vec3 color = glm::vec3(1);

    void
    send(const Program & pgrm) const
    {
        const GLuint pgrm_id = pgrm.getGLId();
        GLuint trans_loc = glGetUniformLocation(pgrm_id,
                                                "trans");
        GLuint color_loc = glGetUniformLocation(pgrm_id,
                                                "color");
        glUniformMatrix4fv(trans_loc,
                           1,
                           GL_FALSE,
                           glm::value_ptr(transform));
        glUniform3fv(color_loc,
                     1,
                     glm::value_ptr(color));
    }

    Instance
    translate(glm::vec2 pos) const
    {
        Instance copy = *this;
        copy.transform = glm::translate(transform, glm::vec3(pos, 0));
        return copy;
    }

    Instance
    translate(float x, float y) const
    {
        return translate({x, y});
    }


    Instance
    scale(float s) const
    {
        Instance copy = *this;
        copy.transform = glm::scale(transform, glm::vec3(s, s, s));
        return copy;
    }


    Instance
    colored(glm::vec3 col) const
    {
        Instance copy = *this;
        copy.color = col;
        return copy;
    }

    Instance
    colored(float r, float g, float b) const
    {
        return colored({r, g, b});
    }
};





struct Scene
{
    std::unordered_map<GLObject, std::vector<Instance>> objects;

    //Scene() = default;

    void
    push(GLObject obj, Instance inst)
    {
        if (objects.find(obj) == objects.end())
        {
            objects.insert({obj, {}});
        }
        objects[obj].push_back(inst);
    }

    void
    draw()
    {
        for (const auto & pairs: objects)
        {
            auto obj = pairs.first;
            for (const auto & inst: pairs.second)
            {
//                inst.send();
                obj.draw();
            }
        }
    }
};


struct Renderer
{ 
    Program pgrm;
    GLenum draw_mode = GL_TRIANGLES;

    void
    setup()
    {
        pgrm.use();
    }

    Renderer(std::string path_vs,
             std::string path_fs,
             GLenum mode = GL_TRIANGLES):
        pgrm(loadProgram(path_vs, path_fs)),
        draw_mode(mode)
   {}

    GLuint
    getGLId()
    {
        return pgrm.getGLId();
    }

    void
    draw(const Scene & scene)
    {
        for (const auto & pair: scene.objects)
        {
            draw_instances(pair.first, pair.second);
        }
    }
    
    void
    draw_instances(const GLObject & obj, const std::vector<Instance> & instances)
    {
        for (const auto & inst: instances)
        {
            inst.send(pgrm);
            obj.draw();
        }
        
    }
    
};



    


/// creates a circle with a given number of vertices
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
    GLObject obj(verts, GL_LINE_LOOP);
    return obj;
}

GLObject
disk(int n_points)
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
    std::vector<Vertex> verts2;

    Vertex v0 = {glm::vec2(0, 0)};

    for (int i = 0; i < n_points; ++i)
    {
        int j = (i+1) % n_points;
        verts2.push_back(verts[i]);
        verts2.push_back(verts[j]);
        verts2.push_back(v0);

    }
    
    GLObject obj(verts2, GL_TRIANGLES);
    return obj;
}

namespace CGA
{
    Mvec e0() { return c2ga::e0<float>(); }
    Mvec e1() { return c2ga::e1<float>(); }
    Mvec e2() { return c2ga::e2<float>(); }
    Mvec ei() { return c2ga::ei<float>(); }

    Mvec
    vector(float x, float y)
    {
        return x*CGA::e1() + y*CGA::e2();
    }

    Mvec
    unit_vector(float angle)
    {
        return CGA::vector(std::cos(angle), std::sin(angle));
    }
    
    /** Creates a point of given center and radius */
    Mvec
    circle(Mvec center, float radius)
    {
        return !(center - radius*radius/2.0*c2ga::ei<float>());
    }

    /** Creates a line passing through the two given points */
    Mvec
     line(Mvec p1, Mvec p2)
    {
        return p1 ^ p2 ^ c2ga::ei<float>();
    }

    /**
       Creates a point of given x and y coordinates
    */
    Mvec
    point(float x, float y)
    {
        return c2ga::point<float>(x, y);
    }

    /**
       Creates a point of given glm coordinates
    */
    Mvec
    point(glm::vec2 pos)
    {
        return CGA::point(pos[0], pos[1]);
    }

    Mvec
    translator(const Mvec & vector)
    {
        return 1 - 0.5*vector*CGA::ei();
    }

    Mvec
    rotor(const float angle)
    {
        auto a = angle/2.0;
        return std::cos(a) - std::sin(a)*(CGA::e1()^CGA::e2());
    }

    Mvec
    vector_from_point(const Mvec & point)
    {
        Mvec e0i = CGA::e0()^CGA::ei();

        return (point ^ e0i | e0i);
    }
    
    Mvec
    motor(const float angle, const Mvec & point_axis)
    {
        Mvec pos = vector_from_point(point_axis);
        Mvec translate_to_0 = translator(-pos);
        Mvec translate_from_0 = translator(pos);
        Mvec rot = rotor(angle);
        Mvec versor = translate_from_0 * rot * translate_to_0;
        return versor;
        
    }
    
    Mvec
    conjugate_versor(const Mvec & versor)
    {
        auto n = versor.norm();
        auto inv = versor.reverse();
        if (versor.grade() % 2 == 0)
        {
            return inv;
        }
        else
        {
            return -inv;
        }
    }

    Mvec
    apply_versor(const Mvec & versor, const Mvec & obj)
    {
        auto conj = CGA::conjugate_versor(versor);
        std::cout << "V  = " << versor << std::endl; 
        std::cout << "V~ = " << conj << std::endl; 
        std::cout << "VP " << versor*obj << std::endl;
        std::cout << "VPV~ " << versor*obj*conj << std::endl;
        return versor*obj*conj;
    }

    /**
       Creates a random point
    */
    Mvec
    rand_point()
    {
        std::random_device rd;  //Will be used to obtain a seed for the random number engine
        std::mt19937 gen(rd()); //Standard mersenne_twister_engine seeded with rd()
        std::uniform_real_distribution<> dis(-1.0, 1.0);

        return CGA::point(dis(gen), dis(gen));
    }
 
    /**
       Creates a random line
    */
    Mvec
    rand_line()
    {
        return CGA::line(CGA::rand_point(),
                         CGA::rand_point());
    }


    /**
       Returns the intersection of two objects.
       Made to intersect lines together, not studied for other cases
     */
    Mvec
    intersect(Mvec obj1, Mvec obj2)
    {
        return !((!obj1) ^ (!obj2));
    }


    /** Create a GL-renderable object from the given Mvec, supposed to be a point */
    GLObject
    point_object(Mvec point)
    {
        auto xy = glm::vec2(point[c2ga::E1], point[c2ga::E2]);
        auto size = 0.01;

        auto p00 = xy - glm::vec2(size, size);
        auto p01 = xy - glm::vec2(size, -size);
        auto p11 = xy + glm::vec2(size, size);
        auto p10 = xy + glm::vec2(size, -size);
        
        std::vector<Vertex> vertices =
            {
                {p00}, {p11},
                {xy},
                {p01}, {p10}
            };
        return GLObject(vertices);
        
    }

    Mvec
    extractFlatPoint(const Mvec &flatPoint){
        // Euclidean vector
        auto pt = - ( c2ga::e0i<float>() | (c2ga::e0<float>() ^ flatPoint) )/flatPoint[c2ga::E0i];
        // remove numerical error (nearly zero remaining parts)
        pt.roundZero();

        // make a point from this Euclidean vector
        pt = point(pt[c2ga::E1],pt[c2ga::E2]);  // e0 + x e1 + ye2 + 1/2 ...
        return pt;
    }

    std::array<Mvec, 2>
    extractPairPoint(const Mvec & pair)
    {
        auto d = sqrt(pair|pair);
        auto div = -c2ga::ei<float>() | pair;
        auto p1 = (pair - d)/div;
        auto p2 = (pair + d)/div;
        return {p1, p2};
    }

    Mvec
    intersect_lines(const Mvec & l1, const Mvec & l2)
    {
        auto inter = intersect(l1, l2);
        return extractFlatPoint(inter);
    }

    glm::vec2
    point_to_glm(const Mvec & point)
    {
        return glm::vec2(point[c2ga::E1],
                         point[c2ga::E2]);
    }

    /**
       We intersect the line with a circle bigger than the screen
       so that the vertices are never too far away.
       Can be optimized
    */
    GLObject
    line_object(Mvec line)
    {
        auto cir = CGA::circle(CGA::point(0, 0), 2);
        auto inters = extractPairPoint(intersect(cir, line));
        glm::vec2 p1 = point_to_glm(inters[0]);
        glm::vec2 p2 = point_to_glm(inters[1]);


        std::vector<Vertex> vertices = {{p1}, {p2}};
        return GLObject(vertices, GL_LINES);
        
        
    }
    GLObject
    segment_object(const Mvec & pair)
    {
        auto pts = extractPairPoint(pair);
        glm::vec2 p1 = point_to_glm(pts[0]);
        glm::vec2 p2 = point_to_glm(pts[1]);
        std::vector<Vertex> vertices = {{p1}, {p2}};
        return GLObject(vertices, GL_LINES);
    }


    /** A segment can be seen as a point pair */
    Mvec
    intersect_segment(const Mvec & pair, const Mvec & line)
    {
        auto inter = CGA::extractFlatPoint(CGA::intersect(pair^CGA::ei(), line));
        if (inter < EPSILON)
            return 0.0;

        auto div = -CGA::ei() | pair;
        auto center = pair / div;
        auto d = inter - center;
        auto dist2 = d|d;
        float div2 = div | div; 
        auto radius2 = (pair | pair)/div2;

        if (dist2 <= radius2)
        {
            return inter;
        }
        else
        {
            return 0.;
        }
        
    }

    Mvec
    line_director_vec(const Mvec & line)
    {
        auto u = CGA::e1()*line[c2ga::E01i] + CGA::e2()*line[c2ga::E02i];

        return u / u.norm();
    }

    Mvec
    line_normal(const Mvec & line)
    {
        auto n = CGA::e1()*line[c2ga::E02i] - CGA::e2()*line[c2ga::E01i];
        return n/n.norm();
    }
    
    Mvec
    reflect(const Mvec & ray_dir, const Mvec & normal, const Mvec & inter)
    {
        auto reflector = normal;

        auto reflected = reflector * ray_dir * CGA::conjugate_versor(reflector);
        //auto refl_dir = CGA::line_director_vec(reflected);
        //auto reflected2 = inter ^ reflected ^ CGA::ei();
        return reflected/reflected.norm();
    }

    /** 
        indice_ratio = n2/n1
        n1 sin(angle in) = n2 sin (angle out)
        (angles relative to the normal line)
        angle out = Arcsin(1/indice_ratio * sin(angle in))
        
        source for the CGA formula: https://pure.uva.nl/ws/files/4375498/52687_fontijne.pdf

     */
    Mvec
    refract(const Mvec & l, const Mvec & n, float r)
    {
        float c = n | l;
        float sign = ((c >= 0.0)?1.0:-1.0);
        c = c * sign;
        auto b = r*c - sqrt(1-r*r*(1-c*c));
        auto refrac = (r*l - sign*b*n);
        auto normalized = refrac/refrac.norm();
        
        std::cout << "      sign: " << sign << std::endl;
        std::cout << "        l : " << l << std::endl;
        std::cout << "        n : " << n << std::endl;
        std::cout << "        c : " << c << std::endl;
        std::cout << "        r : " << r << std::endl;
        std::cout << "        b : " << b << std::endl;
        std::cout << "  b² + r² = " << b*b+r*r << std::endl;
        std::cout << "Refraction: " << refrac << std::endl;
        std::cout << "normalized: " << normalized << std::endl;

        return refrac/refrac.norm();
    }

    float
    radius(const Mvec & round_thing)
    {
        auto conj = !round_thing;
        auto r = conj/(conj | CGA::ei()); 
        return r|r;
    }

    // returns director vectors, not lines
    std::array<Mvec, 2>
    ray_on_segment(const Mvec & ray, const Mvec & seg)
    {
        
        auto inter = CGA::intersect_segment(seg, ray);

        auto ray_dir = CGA::line_director_vec(ray);
        
        auto normal = line_normal(seg^CGA::ei());

        // a line that is parallel to the reflexion
        auto reflected = reflect(ray_dir, normal, inter);

        // ratio = n1/n2
        auto refracted = refract(ray_dir, normal, 0.5);
        
        std::cout << "#### ray   " << ray << std::endl;
        std::cout << "#### seg   " << seg << std::endl;
//        std::cout << "#### flat  " << inter_flat << std::endl;
        std::cout << "#### inter " << inter << std::endl;
        std::cout << "#### norma " << normal << std::endl;
        std::cout << "#### refle " << reflected << std::endl;
        std::cout << "#### refract " << refracted << std::endl;

        
        return {reflected, refracted};
    }
    
    void
    ray_on_sphere(const Mvec & ray, const Mvec & sphere)
    {
        auto inter = intersect(ray, sphere);
        if (inter > EPSILON)
        {
            auto intersections = extractPairPoint(inter);
            auto int1 = intersections[0];
            auto int2 = intersections[1];

            auto u = CGA::line_director_vec(ray);
        }
        else
        {
            std::cout << "no intersection" << std::endl;
        }
    }

    
    
};

struct Ray
{
    Mvec ray;
    Mvec dir;
    Mvec source;
    glm::vec3 color;

    Ray(Mvec pt, Mvec dir, glm::vec3 col):
        ray(pt ^ dir ^ CGA::ei()),
        dir(dir/dir.norm()),
        source(pt),
        color(col)
    {}

    template<class Obj>
    std::optional<std::array<Ray, 2>>
    reflect(const Obj & object) const;

    GLObject
    to_gl_object() const
    {
        auto outer_cir = CGA::circle(CGA::e0(), 2.0);
        auto far = CGA::extractPairPoint(CGA::intersect(ray, outer_cir))[1];
        return CGA::segment_object(source^far);
    }
};

struct Segment
{
    Mvec segment;
    Mvec normal;
    Mvec tangent;
    glm::vec3 color;
    float optic_ratio = 1.0;
    
    Segment(glm::vec2 vert1, glm::vec2 vert2, glm::vec3 col, float ratio):
        segment(CGA::point(vert1)^CGA::point(vert2)),
        normal(CGA::line_normal(CGA::point(vert1)^CGA::point(vert2)^CGA::ei())),
        tangent(CGA::line_director_vec(CGA::point(vert1)^CGA::point(vert2)^CGA::ei())),
        color(col),
        optic_ratio(ratio)
    {}
};


template<>
std::optional<std::array<Ray, 2>>
Ray::reflect(const Segment & segment) const
{
    auto inter = CGA::intersect_segment(segment.segment, ray);
    auto v = (inter - source);
    float sign = v | dir;

    if (sign > EPSILON && inter.grade() > 0)
    {
        // a line that is parallel to the reflexion
        auto reflected = CGA::reflect(dir, segment.normal, inter);

        // ratio = n1/n2
        auto refracted = CGA::refract(dir, segment.normal, segment.optic_ratio);

        /*
          Ireflect + Irefract + Iabsorb = Iinput

         */
        const float cos_theta_in = std::abs(dir|segment.tangent);
        const float cos_theta_tr = std::abs(refracted|segment.tangent);
        const float nu = segment.optic_ratio;

        const float x = cos_theta_in/cos_theta_tr;

        // power reflected (angular part)
        const float R = 1.0 - 2.0*(std::pow(1+nu*x, -2) + std::pow(x+nu, -2));
        const float T = 1.0 - R;

        auto absorbed = color*segment.color;
        
        auto non_absorbed = (glm::vec3(1.0) - absorbed);

        auto col_reflected = non_absorbed*R; 
        auto col_transmited = non_absorbed*T;
        
        auto ray_reflected = Ray(inter, reflected, col_reflected);
        auto ray_transmited = Ray(inter, refracted, col_transmited);
        std::array<Ray, 2> result = {ray_reflected, ray_transmited};
        return result;
    }
    else
    {
        return {};
    }
}




void
ray_tracer(Scene & scene)
{
    auto pas = 0.1;
    auto cir = disk(10);
    Instance inst;
    auto segment = CGA::point(0.5, -1.0) ^ CGA::point(0.5, 1.0);
    scene.push(CGA::segment_object(segment), inst.colored(1., 1., 0.));

    auto pt_source = CGA::point(-0.5, 0.5);
    {
        auto pos =  CGA::point_to_glm(pt_source);
        auto inst2 = inst.translate(pos).colored(1., 1., 0.).scale(0.01);
        scene.push(cir, inst2);
    }

    std::vector<Mvec> rays;

    {
        Mvec e0i = CGA::e0()^CGA::ei();
        auto vec = (pt_source^e0i | e0i);
        auto trans = CGA::translator(vec);
        auto point = CGA::point(0.1, 0.0);
        auto pos0 =  CGA::point_to_glm(pt_source);

        for (auto alpha = 0.; alpha < 2.0*M_PI; alpha+=pas)
        {
            Instance inst;
            auto rotor = CGA::rotor(alpha);
            auto point2 = CGA::apply_versor(trans*rotor, point);

            auto pos = CGA::point_to_glm(point2);

            scene.push(cir, inst.translate(pos)
                       .colored(1., 1., 1.)
                       .scale(0.002));

            auto ray = pt_source ^ point2 ^ CGA::ei();

            auto outer_cir = CGA::circle(CGA::e0(), 2.0);

            auto far = CGA::extractPairPoint(CGA::intersect(ray, outer_cir))[1];

            auto obj = CGA::segment_object(point2^far);
            {
                Instance inst;
                scene.push(obj, inst.colored(0., 1., 0.));
            }

            rays.push_back(ray);
        
       } 

    }


    for (const auto & ray: rays)
    {
        auto inter = CGA::intersect_segment(segment, ray);
        auto u = CGA::line_director_vec(ray);
        auto v = (inter - pt_source);
        float sign = v | u;

        if (sign > EPSILON && inter.grade() > 0)
        {
            {
                Instance inst;
            
                scene.push(cir, inst
                           .translate(CGA::point_to_glm(inter))
                           .colored(1, 0, 0)
                           .scale(0.01));
            }

            {
                auto refs = CGA::ray_on_segment(ray, segment);

                auto reflexion = refs[0];
                auto refraxion = refs[1];

                float t = 0.1;
                auto inter_v = CGA::point_to_glm(inter);
                auto refle_v = inter_v + glm::normalize(CGA::point_to_glm(reflexion))*t;
                auto refra_v = inter_v + glm::normalize(CGA::point_to_glm(refraxion))*t;

                std::vector<Vertex> vert_refle = {{inter_v}, {refle_v}};
                std::vector<Vertex> vert_refra = {{inter_v}, {refra_v}};
            
                GLObject obj_refle(vert_refle, GL_LINES);
                GLObject obj_refra(vert_refra, GL_LINES);
            
                Instance inst2;
            
                scene.push(obj_refle, inst2.colored(0., 0.1, 1.));
                scene.push(obj_refra, inst2.colored(1., 0.1, 0.));

                auto n = !(segment^CGA::ei());
                Instance inst3;
                scene.push(CGA::line_object(inter^n^CGA::ei()), inst3.colored(0.5, 0.5, 0.5));
                
                
            }            
            
        }
    }
    
}

void
ray_tracer2(Scene & scene)
{
    auto pas = 0.1;
    auto cir = disk(10);
    Instance inst;
    auto segment = Segment({0.5, -0.8}, {0.7, 0.9}, {0.1, 0.2, 0.2}, 0.5);

    scene.push(CGA::segment_object(segment.segment), inst.colored(segment.color));
    
    
    auto pt_source = CGA::point(-0.5, 0.5);
    {
        auto pos =  CGA::point_to_glm(pt_source);
        Instance inst;
        auto inst2 = inst.translate(pos).colored(1., 1., 0.).scale(0.01);
        scene.push(cir, inst2);
    }
    
    std::vector<Ray> rays;

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

        auto ray = Ray(point2, dir, {0.9, 0.8, 0.2});

        scene.push(ray.to_gl_object(), inst.colored(ray.color));
            
        rays.push_back(ray);
        
    } 


    for (const auto & ray: rays)
    {
        // loop over segments (one for now)
        {
            const auto inter = ray.reflect(segment);
            if (inter)
            {
                const auto refl = inter.value()[0];
                const auto refr = inter.value()[1];

                scene.push(refl.to_gl_object(), inst.colored(refl.color));
                scene.push(refr.to_gl_object(), inst.colored(refr.color));
                
            }
        }
    }

    
}

    

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

    auto path_vs = applicationPath.dirPath() + "shaders" + "2D_instancied.vs.glsl";
    auto path_fs = applicationPath.dirPath() + "shaders" + "2D.fs.glsl";

    std::cout << path_vs << std::endl;
    Renderer rdr(
        path_vs,
        path_fs,
        GL_LINE_LOOP
        );

    Scene scene;

    if (false)
    {
        std::vector<Mvec> lines;
        for (auto i = 0; i < 10; ++i)
        {
            auto line = CGA::rand_line();
            lines.push_back(line);
            scene.push(CGA::line_object(line), {glm::mat4(1)});
        }

        auto segment = CGA::point(0.2, 0.2) ^ CGA::point(0.3, 0.4);
        Instance seg_inst;
        scene.push(CGA::segment_object(segment), seg_inst.colored(1., 1., 0.));
        auto cir = disk(10);
        auto scale_f = 0.01;
        auto scale = glm::scale(glm::mat4(1), glm::vec3(scale_f, scale_f, scale_f));
        for (auto & l1: lines)
        {
            for (auto & l2: lines)
            {
                if (l1 != l2)
                {
                    auto inter = CGA::intersect(l1, l2);
                    auto non_flat = CGA::extractFlatPoint(inter);
                    auto pos =  glm::vec2(non_flat[c2ga::E1], non_flat[c2ga::E2]);

                    Instance inst;
                    inst.translate(pos)
                        .scale(0.01)
                        .colored(1.0, 0.2, 0.2);

                    scene.push(cir, inst);
                
                
                }
            }
        }

        // intersect with segment
        for (auto & l: lines)
        {
            auto inter = CGA::intersect_segment(segment, l);
            if (inter|inter > EPSILON)
            {
                std::cout << "SEG INTERSECTION?  " << inter << std::endl;
                auto pos =  CGA::point_to_glm(inter);

                Instance inst;
                inst.translate(pos)
                    .scale(0.01)
                    .colored(0.2, 1.0, 0.2);

                scene.push(cir, inst);

                auto tang = CGA::line_director_vec(segment ^ CGA::ei());

                std::cout << "real inter  " << inter <<std::endl;
                auto refs = CGA::ray_on_segment(l, segment);

                auto reflexion = refs[0];
                auto refraxion = refs[1];

                float t = 0.1;
                auto inter_v = CGA::point_to_glm(inter);
                auto refle_v = inter_v + glm::normalize(CGA::point_to_glm(reflexion))*t;
                auto refra_v = inter_v + glm::normalize(CGA::point_to_glm(refraxion))*t;

                std::vector<Vertex> vert_refle = {{inter_v}, {refle_v}};
                std::vector<Vertex> vert_refra = {{inter_v}, {refra_v}};
            
                GLObject obj_refle(vert_refle, GL_LINES);
                GLObject obj_refra(vert_refra, GL_LINES);
            
                Instance inst2;
            
                scene.push(obj_refle, inst2.colored(0., 0.1, 1.));
                scene.push(obj_refra, inst2.colored(1., 0.1, 0.));

                //reflect(const Mvec & ray, const Mvec & normal, const Mvec & inter)
            
            }
        }
    }
    else
    {
        ray_tracer2(scene);
    }
    // test
    {
        auto p1 = CGA::rand_point();
        auto p2 = CGA::rand_point();
        auto seg1 = p1^p2;
        auto seg2 = p2^p1;

        auto div1 = - CGA::ei()|seg1;
        auto div2 = - CGA::ei()|seg2;
        
        auto c1 =  seg1/div1;
        auto c2 =  seg1/div2;

        auto norm1 = seg1 | seg1;
        auto norm2 = seg2 | seg2;
        
        auto u1 = sqrt(norm1)/div1;
        auto u2 = sqrt(norm2)/div2;

        std::cout << "div1: " << div1 << std::endl;
        std::cout << "div2: " << div2 << std::endl;
        std::cout << "norm1: " << norm1 << std::endl;
        std::cout << "norm2: " << norm2 << std::endl;
        std::cout << "u1: " << u1 << std::endl;
        std::cout << "u2: " << u2 << std::endl;
        
    }
    /// are lines oriented?
    {
        auto p1 = CGA::rand_point();
        auto p2 = CGA::rand_point();
        auto line1 = p1^p2^CGA::ei();
        auto line2 = p2^p1^CGA::ei();


        auto u0 = p1 - p2;

        auto dx = (p1[c2ga::E1] - p2[c2ga::E1]);
        auto dy = (p1[c2ga::E2] - p2[c2ga::E2]);
        auto d2 = dx*dx + dy*dy;
        auto d = sqrt(d2);

        auto u1 = line1 | (CGA::e0()^CGA::ei());
        auto u2 = line2 | (CGA::e0()^CGA::ei());

        auto ll11 = p1 ^ u1 ^ CGA::ei();
        auto ll12 = p1 ^ u2 ^ CGA::ei();

        auto n1 = CGA::e1()*line1[c2ga::E02i] - CGA::e2()*line1[c2ga::E01i];
        
        std::cout << "line1:  " << line1 << std::endl;
        std::cout << "line2:  " << line2 << std::endl;
        std::cout << "line1*: " << !line1 << std::endl;
        std::cout << "line2*: " << !line2 << std::endl;
        std::cout << "u1:  " << u1 << std::endl;
        std::cout << "n1:  " << n1 << std::endl;
        std::cout << "n1:  " << (CGA::ei()|n1) << std::endl;
        std::cout << "u2:  " << u2 << std::endl;
        std::cout << "u1*: " << !u1 << std::endl;
        std::cout << "u2*: " << !u2 << std::endl;
        std::cout << "u  : " << u0 << std::endl;
        std::cout << "test  : " << (line1 | (CGA::e0()^CGA::ei())) << std::endl;
        std::cout << "d = " << d << "     d2 = " << d2 << std::endl;
        std::cout << "dx = " << dx << "     dy = " << dy << std::endl;
        std::cout << "ll11*: " << ll11 << std::endl;
        std::cout << "ll12*: " << ll12 << std::endl;

        // ll11 = line2
        // ll12 = line1
        // p1^u2^ei = line1
        // p1^()^ei = p2^p1^ei
    }    
    auto pgrm = loadProgram(path_vs, path_fs);

    pgrm.use();

    glEnable(GL_DEPTH_TEST);
    



    for (auto pair: scene.objects)
    {
        std::cout << "OBJ: " << pair.first << " instances = " << pair.second.size() << std::endl;
    }
    {
        std::cout << "TESTS\n\n";

        auto e0 = CGA::e0();
        auto e1 = CGA::e1();
        auto e2 = CGA::e2();
        auto ei = CGA::ei();
        
        std::cout << "e0" << std::endl;
        std::cout << e0*e0*e0 << std::endl;
        std::cout << e1*e0*e1 << std::endl;
        std::cout << e2*e0*e2 << std::endl;
        std::cout << ei*e0*ei << std::endl;
        std::cout << "e2" << std::endl;
        std::cout << e0*e1*e0 << std::endl;
        std::cout << e1*e1*e1 << std::endl;
        std::cout << e2*e1*e2 << std::endl;
        std::cout << ei*e1*ei << std::endl;
        std::cout << "e2" << std::endl;
        std::cout << e0*e2*e0 << std::endl;
        std::cout << e1*e2*e1 << std::endl;
        std::cout << e2*e2*e2 << std::endl;
        std::cout << ei*e2*ei << std::endl;
        std::cout << "ei" << std::endl;
        std::cout << e0*ei*e0 << std::endl;
        std::cout << e1*ei*e1 << std::endl;
        std::cout << e2*ei*e2 << std::endl;
        std::cout << ei*ei*ei << std::endl;

        std::cout << std::endl;
        std::cout << std::endl;

        auto e0i = e0^ei;
        
        std::cout << (CGA::point(3.0, 5.0)^e0i | e0i)<< std::endl;
        

        Mvec v = 1.0;
        auto rev = v.reverse();
        auto n2 = v.quadraticNorm();
        auto n = v.norm();
        auto n22 = rev | v;
        
        std::cout << "v  = " << v << std::endl;
        std::cout << "n2 = " << n2 << std::endl;
        std::cout << "n  = " << n << std::endl;
        std::cout << "rev  = " << rev << std::endl;
        std::cout << "n22  = " << n22 << std::endl;

        auto e01i = CGA::e0()^CGA::e1()^CGA::ei();
        // auto e0i = CGA::e0()^CGA::ei();
        std::cout << (e01i | e0i) << std::endl;
        
    }
    // testing flat point extraction
    {
        for (auto i = 0; i < 1000; i++)
        {
            auto p = CGA::rand_point();

            auto flat = p^CGA::ei();

            auto extracted = CGA::extractFlatPoint(flat);

            auto diff = (p - extracted).norm();
            if (diff > EPSILON)
            {
                std::cout << "DIFF " << diff << std::endl;
                std::cout << "original " << p << std::endl;
                std::cout << "flat     " << flat << std::endl;
                std::cout << "deflated " << extracted << std::endl;
            }
            
        }
    }
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
        //rdr.setup();
        rdr.draw(scene);
        
        // Update the display
        windowManager.swapBuffers();
    }

    
    return EXIT_SUCCESS;
}

