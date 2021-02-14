#version 300 es
precision mediump float;

in vec3 vNormale;
in vec2 vUV;

out vec3 fFragColor;

uniform sampler2D uTexture;

void main()
{
   fFragColor = texture(uTexture, vUV).rgb;
}

