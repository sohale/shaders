struct Ray
{
    vec3 org;
    vec3 dir;
};

const vec3 e1=vec3(1.0, 0.0, 0.0);
const vec3 e2=vec3(0.0, 1.0, 0.0);
const vec3 e3=vec3(0.0, 0.0, 1.0);
const vec3 o0=vec3(0.0, 0.0, 0.0);

const vec4 e4=vec4(0.0, 0.0, 0.0, 1.0);


mat3 getobj() {
    // vec4 loc=vec4(0.0, 0.0, 0.0, 0.0);

    mat3 obj=mat3(e1, e2*2.0, e3);
    vec3 center= e3*6.0;  // Ellipsoid

    return obj;
}


bool solveQuadratic(in vec3 abc, out float x0, out float x1)
{
    float discr = abc.y * abc.y - 4.0 * abc.x * abc.z;
    if (discr < 0.0) return false;
    else if (discr == 0.0) x0 = x1 = - 0.5 * abc.y / abc.x;
    else {
        float q = (abc.y > 0.0) ?
            -0.5 * (abc.y + sqrt(discr)) :
            -0.5 * (abc.y - sqrt(discr));
        x0 = q / abc.x;
        x1 = abc.z / q;
    }
    if (x0 > x1) {
        float t = x0;
        x0 = x1;
        x1 = t;
    }

    return true;
}

bool intersect(in Ray ray, vec3 center, out float t)
{
    const float radius2 = 1.0;

    // analytic solution
    vec3 L = ray.org - center;
    float a = ray.dir.x * ray.dir.x + ray.dir.y * ray.dir.y + ray.dir.z * ray.dir.z ;
    float b = 2.0 * (ray.dir.x*L.x + ray.dir.y*L.y + ray.dir.z*L.z);
    float c = (L.x*L.x + L.y*L.y + L.z*L.z) - radius2;
    vec3 abc = vec3(a,b,c);

    float t0, t1; // solutions for t if the ray intersects
    if (!solveQuadratic(abc, t0, t1)) return false;

    if (t0 > t1) {
        float tt = t0;
        t0 = t1;
        t1 = tt;
    }

    if (t0 < 0.0) {
        t0 = t1; // if t0 is negative, let's use t1 instead
        if (t0 < 0.0) return false; // both t0 and t1 are negative
    }

    t = t0;

    return true;
}

float min(vec2 v) {
    return v.x > v.y ? v.x : v.y;
}

const mat3 camera_screen_mat = mat3(e1, e2, o0);
const vec3 camera_screen_center = -e3;
const vec3 camera_origin = camera_screen_center - 5.0*e3;

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float time = iGlobalTime;
    float mindim = min(iResolution.xy);
    vec2 uv = fragCoord.xy / mindim;
    vec3 uv3 = vec3(uv.x, uv.y, 0);
    vec3 s = camera_screen_mat * uv3 + camera_screen_center; //screen

    Ray r;
    r.org = camera_origin;
    r.dir = uv3 - camera_origin;

    vec3 center= e3*6.0;  // Ellipsoid


    mat3 obj = getobj();
    // mat4 invobj = inverse(obj);

    float t;
    intersect(r, center, t);
    t = t / 10.0;

    // fragColor = vec4(uv,0.5+0.5*sin(time),1.0);
    fragColor = vec4(t,t,t,1.0);
}

