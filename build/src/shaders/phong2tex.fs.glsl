#version 300 es
precision mediump float;

in vec3 vNormale;
in vec2 vUV;
in vec3 fPos;
in vec3 fVPos;

out vec3 fFragColor;

vec3 phong(vec3 diffuse,
           vec3 ambiant,
           vec3 specular,
           float spec_exponent,
           vec3 light_intensity,
           vec3 light_direction)
{

   float diff = min(max(dot(vNormale, light_direction), 0.0), 1.);
   vec3 amb = ambiant*light_intensity;

   vec3 camera_dir = normalize(-fVPos);
   vec3 half_vec = (camera_dir + light_direction)/2.;

   float spec = max(0.0,
                    pow(dot(half_vec, vNormale),
                        spec_exponent));

   return (ambiant+ diff*diffuse+ spec*specular)*light_intensity;
}



struct Light
{
   vec3 dir;
   vec3 col;
};

uniform sampler2D uTexture1;
uniform sampler2D uTexture2;

uniform mat2x3[10] uDirLights;
uniform mat2x3[10] uPtsLights;
uniform int uDirLightsCount;
uniform int uPtsLightsCount;


void main()
{
//    vec3 diffuse_color = (texture(uTexture1, vUV).rgb + texture(uTexture2, vUV).rgb);
      vec3 diffuse_color = texture(uTexture1, vUV).rgb
          + texture(uTexture2, vUV).rgb;
    
   vec3 color = vec3(0., 0., 0);

   for (int i = 0; i < uDirLightsCount; ++i)
   {
      vec3 ldir = uDirLights[i][0];
      vec3 lcol = uDirLights[i][1];

      float ambiant_frac = 0.2;
      color += phong(diffuse_color*(1.0 - ambiant_frac),
                     diffuse_color*ambiant_frac,
                     vec3(0.1, 0.1, 0.1),
                     0.,
                     lcol,
                     normalize(-ldir));
         
      
   } 
   for (int i = 0; i < uPtsLightsCount; ++i)
   {
      vec3 lpos = uPtsLights[i][0];
      vec3 ldir = normalize(lpos-fVPos);
      float d = distance(lpos, fVPos);
      vec3 lcol = uPtsLights[i][1]/d/d;

      float ambiant_frac = 0.2;
      color += phong(diffuse_color*(1.0 - ambiant_frac),
                     diffuse_color*ambiant_frac,
                     vec3(0.1, 0.1, 0.1),
                     0.,
                     lcol,
                     ldir);
         
      
   } 
   fFragColor = color;
}

