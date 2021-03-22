#version 300 es
precision mediump float;

out vec3 fFragColor;

uniform vec3 color;

void main()
{
    fFragColor = color;
}

