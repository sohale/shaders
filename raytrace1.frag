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


const int SPHERE = 5;

struct Obj
{
    int type;
    vec3 center;
    mat3 matrix;
};

Obj getobj() {
    // vec4 loc=vec4(0.0, 0.0, 0.0, 0.0);
    Obj obj;
    obj.type = SPHERE;
    obj.matrix = mat3(e1, e2*2.0, e3);
    obj.center= e3*6.0;  // Ellipsoid
    return obj;
}


// todo: cleanup and fix careless code.
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

bool sphere_intersect(in Ray ray, vec3 center, out float t, out vec3 where)
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

    where = t * ray.dir + ray.org;
    return true;
}

vec3 sphere_normal(in Obj obj, in vec3 where) {
    vec3 d = vec3(where - obj.center);
    d = normalize(d);
    return d;
}

bool raycast(in Ray ray, in Obj obj, out float t, out vec3 where)
{
    if (obj.type == SPHERE) {
        bool did = sphere_intersect(ray, obj.center, t, where);
        return did;
    }
    // error: unrecognised object type
    return false;
}


vec3 my_reflect(in vec3 ray, in vec3 normal) {
    float cos_ = dot(ray, normal);
    vec3 p = cos_ * normal;
    return 2.0 * p - ray;
}

float my_inner(in vec3 a, in vec3 b) {
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

vec4 phong_material(in vec3 light_dir, in vec3 ray_dir, in vec3 normal) {
    float diffuse = - (light_dir.x * normal.x + light_dir.y * normal.y + light_dir.z * normal.z);
    diffuse = diffuse > 0.0 ? diffuse : 0.0;

    vec3 refl = my_reflect(-light_dir, normal);
    refl = normalize(refl);
    ray_dir = normalize(ray_dir);
    float specular0 = my_inner(refl, -ray_dir);
    specular0 = specular0 > 0.0 ? specular0 : 0.0;
    // specular0 = specular0 > 1.0 ? 1.0:specular0;
    // specular0 = specular0 > 10.0 ? specular0: 0.0;
    float specular = pow(specular0, 5.0);
    // float specular = specular0 * 1.0 - floor(specular0 * 1.0);  // nice debug tool
    return vec4(diffuse,diffuse,specular,1.0);
}

float min(vec2 v) {
    return v.x > v.y ? v.x : v.y;
}

vec2 screen_uv(vec2 fragCoord) {
    float mindim = min(iResolution.xy) / 2.0;
    vec2 center = iResolution.xy / 2.0;
    vec2 uv = (fragCoord.xy - center.xy) / mindim;
    vec2 uv2 = vec2(uv.x, uv.y);
    return uv2;
}

struct Camera {
    mat3 screen_mat;
    vec3 screen_center;
    vec3 origin;
};


Camera init_camera(vec2 mouse) {
    Camera camera;
    camera.screen_mat = mat3(e1, e2, o0);
    camera.screen_center = -e3 + mouse.x * e1 + mouse.y * e2;
    camera.origin = camera.screen_center - 5.0*e3;
    return camera;
}


Ray make_ray(Camera camera, vec2 uv2) {
    vec3 uv3 = vec3(uv2, 0.0);

    vec3 s = camera.screen_mat * uv3 + camera.screen_center; //screen

    Ray r;
    r.org = camera.origin;
    r.dir = s - camera.origin;

    return r;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{
    float time = iGlobalTime;

    vec2 mousexy = screen_uv(iMouse.zw);

    Camera camera = init_camera(mousexy);

    //vec3 uv3 = vec3(screen_uv(fragCoord), 0.0);
    // vec3 s = camera_screen_mat * uv3 + camera_screen_center; //screen

    vec2 uv2 = screen_uv(fragCoord);
    Ray r = make_ray(camera, uv2);

    Obj obj = getobj();

    // mat4 invobj = inverse(obj);

    float t;
    vec3 where;
    bool did = raycast(r, obj, t, where);

    vec3 ray_dir_normalized = normalize(r.dir);

    vec3 radial = where - obj.center;

    vec3 light_dir = vec3(-1.0, -1.0, +1.0);
    light_dir = normalize(light_dir);

    vec4 cc;
    if (did) {
        // // c = t / 5.0;
        //c = -radial.z * 1.0;

        vec3 normal = sphere_normal(obj, where);
        cc = phong_material(light_dir, ray_dir_normalized, normal);
    } else {
        //c = 0.0;  // why omitting this causes apparent noise?
        cc = vec4(0.0, 0.0, 0.0, 1.0);
    }

    // fragColor = vec4(uv,0.5+0.5*sin(time),1.0);
    // fragColor = vec4(c, c, c, 1.0);
    fragColor = cc;
}
