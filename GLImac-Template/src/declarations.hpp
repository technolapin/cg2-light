#pragma once

#include <c2ga/Mvec.hpp>
#include "c2gaTools.hpp"
#include <glimac/SDLWindowManager.hpp>
#include <glimac/Sphere.hpp>
#include <GL/glew.h>
#include <glimac/Program.hpp>
#include <glimac/FilePath.hpp>
#include <glimac/Image.hpp>


#include <iostream>
#include <unordered_map>
#include <memory> // shared ptrs
#include <assert.h>
#include <string>
#include <algorithm>
#include <string.h>
#include <math.h>
#include <random>

using namespace glimac;


const GLuint VERTEX_ATTR_POSITION = 0;
const GLuint VERTEX_ATTR_NORMAL = 1;
const GLuint VERTEX_ATTR_UV = 2;

struct Vertex
{
    glm::vec2 pos;
};

const float EPSILON = 0.000001;
using Mvec = c2ga::Mvec<float>;


#include "GLObject.hpp"
#include "instance.hpp"
#include "cga.hpp"
#include "gl_misc.hpp"

#include "optics.cpp"

#include "tree.hpp"
