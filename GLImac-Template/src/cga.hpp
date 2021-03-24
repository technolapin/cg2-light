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
    unit_vector(float x, float y)
    {
        const auto u = CGA::vector(x, y);
        return u/u.norm();
    }

    Mvec
    unit_vector(float angle)
    {
        return CGA::vector(std::cos(angle), std::sin(angle));
    }

    Mvec
    unit_vector(const Mvec & pt1, const Mvec & pt2)
    {
        return CGA::unit_vector(pt2[c2ga::E1] - pt1[c2ga::E1],
                                pt2[c2ga::E2] - pt1[c2ga::E2]);
    }

    /** Creates a point of given center and radius */
    Mvec
    circle(const Mvec center, const float radius)
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

    float
    dist2(Mvec p1, Mvec p2)
    {
        return - 2.0 * p1|p2;
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
    reflect(const Mvec & ray_dir, const Mvec & normal)
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
        auto reflected = reflect(ray_dir, normal);

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

    
    
}

