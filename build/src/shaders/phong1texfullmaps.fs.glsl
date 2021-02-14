#version 300 es
precision mediump float;

in vec2 vUV;
in vec3 fPos;
in vec3 fVPos;
in mat3 TBN;

out vec3 fFragColor;

vec3 phong(vec3 diffuse,
           vec3 ambiant,
           vec3 specular,
           float spec_exponent,
           vec3 light_intensity,
           vec3 light_direction,
           vec3 fNormale)
{

   float diff = min(max(dot(fNormale, light_direction), 0.0), 1.);
   vec3 amb = ambiant*light_intensity;

   vec3 camera_dir = normalize(-fVPos);
   vec3 half_vec = (camera_dir + light_direction)/2.;

   float spec = max(0.0,
                    pow(dot(half_vec, fNormale),
                        spec_exponent));

   return (ambiant+ diff*diffuse+ spec*specular)*light_intensity;
}



struct Light
{
   vec3 dir;
   vec3 col;
};

uniform sampler2D uTexture;
uniform sampler2D uNormales;
uniform sampler2D uDiffuse;
uniform sampler2D uSpecular;
uniform sampler2D uSpecularLevel;

uniform mat2x3[10] uDirLights;
uniform mat2x3[10] uPtsLights;
uniform int uDirLightsCount;
uniform int uPtsLightsCount;

void main()
{
    vec3 fNormale = texture(uNormales, vUV).rgb;
    fNormale = normalize(fNormale * 2.0 - 1.0);


    
    float ambiant_frac = 0.1;

   vec3 fDiffuse = texture(uDiffuse, vUV).rgb;
   vec3 fSpecular = texture(uSpecular, vUV).rgb;
   float fSpecularLevel = texture(uSpecularLevel, vUV).r;
   
   vec3 color = vec3(0., 0., 0);

   for (int i = 0; i < uDirLightsCount; ++i)
   {
      vec3 ldir = uDirLights[i][0];
      vec3 lcol = uDirLights[i][1];


      ldir = TBN * normalize(-ldir);
      
      color += phong(fDiffuse*(1.0-ambiant_frac),
                     fDiffuse*ambiant_frac,
                     fSpecular,
                     fSpecularLevel,
                     lcol,
                     ldir,
                     fNormale);
         
      
   } 
   for (int i = 0; i < uPtsLightsCount; ++i)
   {
      vec3 lpos = uPtsLights[i][0];
      vec3 ldir = normalize(lpos-fVPos);
      float d = distance(lpos, fVPos);
      vec3 lcol = uPtsLights[i][1]/d/d;

      ldir = TBN * normalize(ldir);

      color += phong(fDiffuse*(1.0-ambiant_frac),
                     fDiffuse*ambiant_frac,
                     fSpecular,
                     fSpecularLevel,
                     lcol,
                     ldir,
                     fNormale );
         
      
   } 
   fFragColor = color;
}

