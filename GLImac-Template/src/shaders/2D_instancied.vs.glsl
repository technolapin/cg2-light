#version 300 es

layout(location = 0) in vec2 aPos;
uniform mat4 trans;

void main() {
   gl_Position = trans*vec4(aPos, 0, 1);
};
