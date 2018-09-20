#define USE_CAMERA 1

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

/*
mat3 inverse(mat3 m) {
  // from https://github.com/stackgl/glsl-inverse/blob/master/index.glsl
  float a00 = m[0][0], a01 = m[0][1], a02 = m[0][2];
  float a10 = m[1][0], a11 = m[1][1], a12 = m[1][2];
  float a20 = m[2][0], a21 = m[2][1], a22 = m[2][2];

  float b01 = a22 * a11 - a12 * a21;
  float b11 = -a22 * a10 + a12 * a20;
  float b21 = a21 * a10 - a11 * a20;

  float det = a00 * b01 + a01 * b11 + a02 * b21;

  return mat3(b01, (-a22 * a01 + a02 * a21), (a12 * a01 - a02 * a11),
              b11, (a22 * a00 - a02 * a20), (-a12 * a00 + a02 * a10),
              b21, (-a21 * a00 + a01 * a20), (a11 * a00 - a01 * a10)) / det;
}
*/


const int SPHERE = 5;

struct Obj
{
    int type;
    vec3 center;
    mat3 forward_matrix;
    mat3 inverse_matrix;
    vec3 rgb;
};

Obj make_ellipsoid(float rx, float ry, float rz) {
    // vec4 loc=vec4(0.0, 0.0, 0.0, 0.0);
    Obj obj;
    obj.type = SPHERE;
    // obj.matrix = mat3(e1, e2*2.0, e3);
    obj.center= e3*6.0 * 0.0;  // Ellipsoid
    // mat3 im = inverse(obj.matrix);
    obj.rgb = vec3(1.0,1.0,1.0);

    // rx = 1.0;
    // ry = 1.0;
    // rz = 1.0;

    // obj.inverse_matrix = inverse(obj.matrix);
    obj.inverse_matrix = mat3(e1 / rx, e2 / ry, e3 / rz);
    obj.forward_matrix = mat3(e1*rx, e2*ry, e3*rz);
    return obj;
}


// todo: cleanup and fix careless code.
bool solveQuadratic(in vec3 abc, out float x0, out float x1)
{
    float discr = abc.y * abc.y - 4.0 * abc.x * abc.z;
    if (discr < 0.0) return false;
    else if (discr == 0.0)
        x0 = x1 = - 0.5 * abc.y / abc.x;
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

bool sphere_intersect(in Ray ray, out float t, out vec3 where)
{
    const float radius2 = 1.0;

    // analytic solution
    vec3 L = ray.org; // - center;
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

mat3 transpose_(mat3 m) {
  return mat3(m[0][0], m[1][0], m[2][0],
              m[0][1], m[1][1], m[2][1],
              m[0][2], m[1][2], m[2][2]);
}

vec3 sphere_normal(in Obj obj, in vec3 where) {
    // vec3 d = vec3(where - obj.center);
    //d = normalize(d);

    vec3 d = obj.inverse_matrix * (where - obj.center);
    d = transpose_(obj.forward_matrix) * d;

    d = normalize(d);
    return d;
}

bool raycast(in Ray ray, in Obj obj, out float t, out vec3 where)
{
    if (obj.type == SPHERE) {
        Ray ray2 = ray;
        ray2.dir = obj.inverse_matrix * ray.dir;
        ray2.org = obj.inverse_matrix * (ray.org - obj.center);
        float new_norm = length(ray2.dir);


        vec3 where2;
        float t2;
        bool did = sphere_intersect(ray2, t2, where2);

        where = obj.forward_matrix * where2 + obj.center;
        t = t2 ; //* new_norm;
        ////////////////
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

vec4 phong_material(in vec3 light_dir, in vec3 ray_dir, in vec3 normal, in vec3 obj_rgb) {
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

    return vec4(diffuse * obj_rgb + specular * vec3(1.0,1.0,1.0), 1.0);
}

float min_(vec2 v) {
    return v.x > v.y ? v.x : v.y;
}

vec2 screen_uv(vec2 fragCoord) {
    float mindim = min_(iResolution.xy) / 2.0;
    vec2 center = iResolution.xy / 2.0;
    vec2 uv = (fragCoord.xy - center.xy) / mindim;
    vec2 uv2 = vec2(uv.x, uv.y);
    return uv2;
}


mat3 rotationMatrixXY(float aXY){
    return mat3(
        cos( aXY ), -sin( aXY ), 0.0,
        sin( aXY ),  cos( aXY ), 0.0,
        0.0,           0.0, 1.0
    );
}
mat3 rotationMatrixYZ(float a){
    return mat3(
        1.0, 0.0, 0.0,
        0.0, cos( a ), -sin( a ),
        0.0, sin( a ),  cos( a )
    );
}
mat3 rotationMatrixXZ(float a){
    return mat3(
        cos( a ), 0.0, -sin( a ),
        0.0,      1.0,  0.0,
        sin( a ), 0.0,  cos( a )
    );
}

/*
mat4 rotationMatrix(xy, yz, xz){
    rotationMatrixXY();
}
*/


struct Camera {
    mat3 screen_mat;
    vec3 screen_center;
    vec3 origin;  // eye point
};

const float PI = 4.1415926536;

Camera init_camera(vec2 mouse) {
    Camera camera;

    //camera.screen_mat = mat3(e1, e2, o0);
    // mat3 rot = rotationMatrixYZ(0.5*PI) * rotationMatrixXZ(0.0);
    // camera.screen_center = -e3 + mouse.x * e1 * 1.0 + mouse.y * e2 * 1.0;

    mat3 rot = rotationMatrixYZ(-mouse.y) * rotationMatrixXZ(-mouse.x);
    camera.screen_mat = rot * mat3(e1, e2, e3);
    camera.screen_center = rot * (-e3) *5.0;
    camera.origin =  (camera.screen_center - rot * 5.0*e3);


    // mat3 m = rotationMatrixYZ()

    return camera;
}

Ray make_ray(Camera camera, vec2 uv2) {
    vec3 uv3 = vec3(uv2, 0.0);

    vec3 s = camera.screen_mat * uv3 + camera.screen_center; //screen

    Ray r;
    r.org = s;
    r.dir = s - camera.origin;
    //r.dir = -s;
    r.dir = normalize(r.dir);

    return r;
}

struct TexturedScreen {
    // mat3 screen_mat;
    // vec3 screen_center;
    vec3 e3t_Minv;
    float e3t_Minv_C0;
};

TexturedScreen make_TexturedScreen_behind_camera(Camera camera) {
    /*
    screen.screen_mat = camera.screen_mat;
    screen.screen_center = camera.screen_center;
    screen.e3t_Minv = transpose_(inverse(screen.screen_mat))*vec3(1.0, 0.0, 0.0);
    screen.e3t_Minv_C0 = my_inner(screen_center, e3t_Minv);
    */

    //copied from camera
    float mouse_y =0.0, mouse_x = 0.0;
    mat3 rot = rotationMatrixYZ(-mouse_y) * rotationMatrixXZ(-mouse_x);
    // camera__screen_mat = rot * mat3(e1, e2, e3);
    vec3 screen_center = rot * (-e3) *5.0;
    mat3 rot_inv = rotationMatrixYZ(+mouse_y) * rotationMatrixXZ(+mouse_x);

    TexturedScreen screen;
   
    screen.e3t_Minv = rot_inv * vec3(0.0, 0.0, 1.0);
    screen.e3t_Minv_C0 = my_inner(screen_center, screen.e3t_Minv);

    return screen;
}

bool project_onto_screen_t(in TexturedScreen screen, in Ray ray, out float t) {
    float denom = my_inner(ray.dir, screen.e3t_Minv);
    if (abs(denom) < 0.000000001)
        return false;
    t = my_inner(ray.org, screen.e3t_Minv) - screen.e3t_Minv_C0;
    return (t >= 0.0);
}

bool project_onto_screen_uv2(in TexturedScreen screen, in Ray ray, out vec2 uv2) {
    float t;
    if (project_onto_screen_t(screen, ray, t)) {
        vec3 uv3 = t * ray.dir + ray.org;
        uv2 = uv3.xy;
        return true;
    }
    return false;
}


vec4 panic() {
    return vec4(1.0, 0.0, 0.0, 1.0);
}

const int num_objects = 4;


bool world_raycast(Ray ray,
    in Obj[num_objects] obj,
    out Obj chosen_obj,
    out vec3 chosen_where,
    out int chosen_obj_id,
    out float tmin,
    in int exclude
) {
    bool did = false;

    chosen_obj_id = -1;

    tmin = 100000000.0;

    //Obj chosen_obj;
    //vec3 chosen_where;
    chosen_obj_id = -1;

    {
        Obj curr_obj;


        for (int i = 0 ; i < num_objects; ++i) {
            /*
            if (i==0) {
                curr_obj = obj[0];
            } else if (i==1) {
                curr_obj = obj[1];
            } else if (i==2) {
                curr_obj = obj[2];
            }
            */
            curr_obj = obj[i];
            // Avoid self intersection
            if (exclude == i)
                continue;

            float t;
            vec3 where;
            bool did1;

            did1 = raycast(ray, curr_obj, t, where);
            if (did1) {
                if (tmin > t)
                {
                    tmin = t;
                    chosen_obj_id = i;
                    chosen_obj = curr_obj;
                    chosen_where = where;
                    did = true;
                }
                // assert did == true
            }
        };

    }


    return did;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{



    float time = iTime;

    vec2 mousexy = screen_uv(iMouse.xy);

    Camera camera = init_camera(-mousexy * 3.14*2.0 + vec2(0.0, 0.6)* 3.14*2.0 + vec2(time*0.02 + (-0.1) + sin(time*0.33 *PI*2.0)*0.2, time*0.004));

    //vec3 uv3 = vec3(screen_uv(fragCoord), 0.0);
    // vec3 s = camera_screen_mat * uv3 + camera_screen_center; //screen

    vec2 uv2 = screen_uv(fragCoord);
    Ray r = make_ray(camera, uv2);

    Obj obj[num_objects];
    obj[0] = make_ellipsoid(1.0, 1.0,1.0);
    obj[0].center.x -= 1.5/1.0;
    obj[0].rgb = vec3(1.0, 1.0, 0.0);

    obj[1] = make_ellipsoid(1.0, 0.7,1.0);
    obj[1].center.x += 1.5/2.0;
    obj[1].center.z -= 0.3;
    obj[1].rgb = vec3(1.0, 0.0, 0.0);

    obj[2] = make_ellipsoid(1.0, 1.0, 0.7);
    obj[2].center.y += 0.5;
    obj[2].center.z += 0.9;
    obj[2].rgb = vec3(0.0, 1.0, 0.0);

    obj[3] = make_ellipsoid(0.3, 0.3, 0.3);
    obj[3].center.y += 0.0;
    obj[3].center.y += 0.2;
    obj[3].center.z += 0.0;
    obj[3].rgb = vec3(0.0, 0.0, 1.0);

    // mat4 invobj = inverse(obj);

    if (false) {
    #ifdef USE_CAMERA
        TexturedScreen webcam_screen = make_TexturedScreen_behind_camera(camera);
    #endif
    }

    //int obj_id = -1;

    Obj chosen_obj;
    vec3 chosen_where;
    int chosen_obj_id;
    float tmin;

    bool did = world_raycast(r, obj, chosen_obj,chosen_where,chosen_obj_id,tmin, -1);

    // not necessary anymore:
    //vec3 ray_dir_normalized = normalize(r.dir);

    // vec3 radial = where - obj[i].center;

    vec3 light_dir = vec3(-1.0, -1.0, +1.0);
    light_dir = normalize(light_dir);

    vec4 cc;
    vec4 cc2;

    if (did) {
        // // c = t / 5.0;
        //c = -radial.z * 1.0;

        if (chosen_obj_id < 0) {
            fragColor = panic();
            return;
        }

        vec3 normal = sphere_normal(chosen_obj, chosen_where);
        cc = phong_material(light_dir, r.dir, normal, chosen_obj.rgb);


        Ray ray2;
        ray2.org = chosen_where + r.dir * 0.001 * 0.0;
        ray2.dir = my_reflect(-r.dir, normal);

        Obj chosen_obj2;
        vec3 chosen_where2;
        int chosen_obj_id2;
        float tmin2;

        bool did2 = world_raycast(ray2, obj, chosen_obj2,chosen_where2,chosen_obj_id2, tmin2, chosen_obj_id);

        if (did2) {
            vec3 normal2 = sphere_normal(chosen_obj2, chosen_where2);
            cc2 =   phong_material(light_dir, ray2.dir, normal2, chosen_obj2.rgb );

            // float w1 = 0.6, w2 = 0.4; //
            float w1 = 1.0, w2 = 0.0; // Pure reflection
            cc2.r = cc2.r * w1 + chosen_obj.rgb.r * w2;
            cc2.g = cc2.g * w1 + chosen_obj.rgb.g * w2;
            cc2.b = cc2.b * w1 + chosen_obj.rgb.b * w2;
            //cc2 = vec4(0.0, 0.0, 0.0, 0.0);

        } else {
            cc2 = vec4(0.0, 0.0, 0.0, 0.0);

            #ifdef USE_CAMERA
                cc2.xyz = texture(iChannel0, ray2.dir.xy).xyz;
                cc2.xyz = texture(iChannel0, 1.0-ray2.dir.xy).xyz - 0.5;
            if (false) {

                //vec3 screen_center = vec3();
                //mat3 screen_matrix = mat3(e1, e2, e3);
                //rayscreen = transform(ray2, screen_matrix, screen_center)
                //vec2 uv;
                //project_onto_screen_uv2(webcam_screen, ray2, uv);
                //cc2.xyz = texture(iChannel0, uv).xyz;
            }
            #endif

        }


        // float tn = abs((tmin -2.5)*1.0);
        // cc = vec4(tn, tn, tn, 0.0) + vec4(0.0,0.0,0.0,1.0);
    } else {
        //c = 0.0;  // why omitting this causes apparent noise?
        cc = vec4(0.0, 0.0, 0.0, 1.0);
        cc2 = vec4(0.0, 0.0, 0.0, 0.0);

        // Background will be from the camera (which is not good)
        // cc.xyz = texture(iChannel0, r.dir.xy).xyz;

    }

    // fragColor = vec4(uv,0.5+0.5*sin(time),1.0);
    // fragColor = vec4(c, c, c, 1.0);
    fragColor = cc * 1.00 + cc2 * 0.4;
}
