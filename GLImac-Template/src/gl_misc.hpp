
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
