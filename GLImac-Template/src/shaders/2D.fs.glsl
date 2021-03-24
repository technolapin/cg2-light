#version 300 es
precision mediump float;

out vec4 fFragColor;

uniform vec3 color;

void main()
{
    float power = dot(color, color)/3.0;
    fFragColor = vec4(color, power);
}

