// Demonstrating my "method of semantic refactoring"
// Starting with pre-existing shader by 01000001 https://www.shadertoy.com/view/4ffXzf

#define saturation .1
// .7 is nice

const float scale = 12.;
const float pi = 3.1415926;

const int n = 4;
vec3 charges[n] = vec3[n](

    vec3(5.5, 5.5, 8),
    vec3(-5, -5, -8),
    vec3(15, -5, 8),
    vec3(-15.5, -5.5, -10)

);

vec2 field(vec2 p){

    vec2 sum = vec2(0);
    
    for (int i = n; i-->0;){
        sum += charges[i].z*normalize(p-charges[i].xy)/(length(p-charges[i].xy)*length(p-charges[i].xy));
    }
    
    return sum;
}

void mainImage( out vec4 O, vec2 U )
{
    vec2 r = iResolution.xy;
    vec2 uv = (2.*U-r)/r.y * scale;
    
    vec2 cell = floor(uv);
    vec2 cuv = (uv-cell) * 2. - 1.;

    if (iMouse.z > 0.) charges[0].xy = (2.*iMouse.xy-r)/r.y * scale;
    charges[0].z = sin(iTime)*8. + 10.;
    charges[1].xy = vec2(sin(iTime*.5), cos(iTime*1.2))*5.5;
    charges[2].xy = vec2(sin(iTime*1.3), cos(iTime*.8))*10.;

    vec2 force = field(cell);
    float x = smoothstep(.3, .0, 
        abs(dot(normalize(force), cuv*mat2(0,1,-1,0)))     // Create bar
    )                                                      // Make it look nice
        * length(force)                                    // brighter at higher force
        * smoothstep(1., .9, length(cuv))                  // limit length of bars to 1 cell width
    ;
    
    float t = atan(force.x/force.y) + (sign(force.y)>0.?-pi:0.);                            // Angle of bar

    O = vec4(sqrt(x) 
        * mix(vec4(1), (cos(t + vec4(0., 2.*pi/3., 4.*pi/3., 0.))*.5+.5), vec4(saturation)) // Colour based on angle
    );
    
    //O /= O + 1.;
}
