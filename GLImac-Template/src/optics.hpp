struct Ray
{
    Mvec ray;
    Mvec dir;
    Mvec source;
    
    glm::vec3 color;

    Ray(const Mvec pt, const Mvec dir, const glm::vec3 col):
        ray(pt ^ dir ^ CGA::ei()),
        dir(dir/dir.norm()),
        source(pt),
        color(col)
    {}




    std::array<std::optional<Ray>, 2>
    reflect(const Mvec & valid_inter,
            const Mvec & normal,
            const Mvec & tangent,
            const glm::vec3 & obj_color,
            const float optic_ratio) const;

    
    template<class Obj>
    std::array<std::optional<Ray>, 2>
    reflect_on(const Obj & object, const Mvec & valid_inter, const Mvec & normal) const;

    template<class Obj>
    std::optional<std::array<Mvec, 2>>
    intersect(const Obj & object) const;
    
    GLObject
    to_gl_object() const
    {
        auto outer_cir = CGA::circle(CGA::e0(), 2.0);
        auto far = CGA::extractPairPoint(CGA::intersect(ray, outer_cir))[1];
        return CGA::segment_object(source^far);
    }

    float
    power2() const
    {
        return glm::dot(color, color);
    }

    Mvec
    infinity() const
    {
        auto outer_cir = CGA::circle(CGA::e0(), 2.0);
        auto far = CGA::extractPairPoint(CGA::intersect(ray, outer_cir))[1];
        return far;
    }
    
};


struct Segment
{
    Mvec segment;
    Mvec normal;
    Mvec tangent;
    
    Segment(const glm::vec2 vert1,
            const glm::vec2 vert2):
        segment(CGA::point(vert1)^CGA::point(vert2)),
        normal(CGA::line_normal(CGA::point(vert1)^CGA::point(vert2)^CGA::ei())),
        tangent(CGA::line_director_vec(CGA::point(vert1)^CGA::point(vert2)^CGA::ei()))
    {}

};






struct Circle
{
    Mvec circle;
    Mvec center;
    float radius;
    
    Circle(const glm::vec2 center,
           const float radius):
        circle(CGA::circle(CGA::point(center), radius)),
        center(CGA::point(center)),
        radius(radius)
    {}

    bool
    is_inside(const Mvec & pt) const
    {
        return CGA::dist2(pt, center) <= radius*radius;
    }
    
    Mvec
    normal(const Mvec & pt) const
    {
        const auto sign = (is_inside(pt))? -1.0:1.0;
        const auto u = sign*(pt - center); // degenerate
        return CGA::unit_vector(u[c2ga::E1], u[c2ga::E2]);
    }

    Mvec
    tangent(const Mvec & pt) const
    {
        const auto n = normal(pt);
        return -n[c2ga::E2]*CGA::e1() + n[c2ga::E1]*CGA::e2();
    }
};


struct SegMesh
{
    std::vector<Mvec> segments;
    // for rought intersection check (can be made smaller)
    Mvec circle;
    
    SegMesh(const std::vector<glm::vec2> & vertices)
    {
        auto center = glm::vec2(0);
        auto radius = 0.0f;
        for (const auto & v: vertices)
        {
            center += v;
        }
        center /= (float) vertices.size();
        for (const auto & v: vertices)
        {
            radius = std::max(radius, glm::distance(center, v));
        }

        circle = CGA::circle(CGA::point(center), radius);

        const auto n = vertices.size();
        for (auto i = 0; i < n; ++i)
        {
            const auto v0 = vertices[i];
            const auto v1 = vertices[(i+1) % n];

            const auto p0 = CGA::point(v0);
            const auto p1 = CGA::point(v1);
            const auto seg = p0 ^ p1;

            segments.push_back(seg);
        }
    }
};
    



enum OpticObjT
{
    SegmentT,
    CircleT,
    SegMeshT
};


struct OpticObject
{
    OpticObjT type;
    glm::vec3 color;
    float optic_ratio;

    union
    {
        Segment _segment;
        Circle _circle;
        SegMesh _mesh;

    };

    OpticObject():
        type(SegmentT),
        color(glm::vec3(0)),
        optic_ratio(1.0),
        _segment(Segment(glm::vec2(0), glm::vec2(0)))
    {
        /*
        {
                SegmentT,
                glm::vec3(0),
                1.0,
                Segment(glm::vec2(0), glm::vec2(0))

            };
        */
    }
    
    OpticObject(const OpticObject & obj):
        type(obj.type),
        color(obj.color),
        optic_ratio(obj.optic_ratio)
    {
        switch (obj.type)
        {
            case SegmentT:
                _segment = obj._segment;
                break;
            case CircleT:
                _circle = obj._circle;
                break;
            case SegMeshT:
                _mesh = obj._mesh;
                break;
        }
    }
    
    ~OpticObject()
    {
        switch (type)
        {
            case SegmentT:
                _segment.~Segment();
                break;
            case CircleT:
                _circle.~Circle();
                break;
            case SegMeshT:
                _mesh.~SegMesh();
                break;
        }
    }

    

    static OpticObject 
    segment(const glm::vec2 & vert1,
            const glm::vec2 & vert2,
            const glm::vec3 & col,
            const float ratio)
    {
        OpticObject obj;
        obj.type = SegmentT;
        obj.color = col;
        obj.optic_ratio = ratio;
        obj._segment = Segment(vert1, vert2);

        return obj;
    }

    static OpticObject 
    circle(const glm::vec2 & center,
           const float radius,
           const glm::vec3 & col,
           const float ratio)
    {
        OpticObject obj;
        obj.type = CircleT;
        obj.color = col;
        obj.optic_ratio = ratio;
        obj._circle = Circle(center, radius);

        return obj;
    }

    static OpticObject 
    mesh(const std::vector<glm::vec2> & vertices,
         const glm::vec3 & col,
         const float ratio)
    {
        OpticObject obj;
        obj.type = SegmentT;
        obj.color = col;
        obj.optic_ratio = ratio;
        obj._mesh = SegMesh(vertices);

        return obj;
    }



    
};
