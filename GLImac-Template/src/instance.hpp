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

        // std::cout << "GL COLOR " << color[0] << " " << color[1] << " " << color[2] << std::endl;
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


