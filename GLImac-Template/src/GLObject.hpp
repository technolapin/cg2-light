#pragma once

#include "declarations.hpp"

struct GLObject
{
private:
    GLuint _vao;
    GLuint _vbo;
    GLuint _size;
    GLuint _mode = GL_TRIANGLES;
    
public:
    GLObject()
    {
    }
    
    GLObject(std::vector<Vertex> & obj, GLuint mode = GL_TRIANGLES):
        _mode(mode)
    {
        _size = obj.size();
      
        glGenBuffers(1, &_vbo);
        glGenVertexArrays(1, &_vao);

        glBindBuffer(GL_ARRAY_BUFFER, _vbo);
        glBindVertexArray(_vao);
                                           
      
        glBufferData(GL_ARRAY_BUFFER,
                     sizeof(Vertex) * obj.size(),
                     obj.data(),
                     GL_STATIC_DRAW);

        glEnableVertexAttribArray(VERTEX_ATTR_POSITION);

        glVertexAttribPointer
            (
                VERTEX_ATTR_POSITION,
                2,
                GL_FLOAT,
                GL_FALSE,
                sizeof(Vertex),
                (const GLvoid*) offsetof(Vertex, pos)
                );
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    void free_gpu()
    {
        glDeleteVertexArrays(1, &_vao);
        glDeleteBuffers(1, &_vbo);
    }

   

    void
    draw() const
    {
        glBindVertexArray(_vao);
        glDrawArrays(_mode, 0, _size);
        glBindVertexArray(0);
    }

    GLuint
    vbo() const
    {
        return _vbo;
    }
    GLuint
    vao() const
    {
        return _vao;
    }
    GLuint
    size() const
    {
        return _size;
    }
    GLuint
    mode() const
    {
        return _mode;
    }
    
    bool operator==(const GLObject &other) const
    {
        return this->vbo() == other.vbo()
            && this->vao() == other.vao()
            && this->size() == other.size()
            && this->mode() == other.mode();
    }
};


namespace std {
  template<> struct hash<GLObject> {
    size_t operator()(GLObject const& obj) const {
        return obj.vbo();
    }
  };
}

std::ostream & operator << (std::ostream &out, const GLObject & obj)
{
    out << "Object [vbo:" << obj.vbo()
        << ", vao:" << obj.vao()
        << ", size:" << obj.size() << "]";
    return out;
}
