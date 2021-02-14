#version 300 es
precision mediump float;

in vec3 vNormale;
in vec2 vUV;

out vec3 fFragColor;

uniform sampler2D uTexture1;
uniform sampler2D uTexture2;

void main()
{
   fFragColor = texture(uTexture1, vUV).rgb
      + texture(uTexture2, vUV).rgb;
}

