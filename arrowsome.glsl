// © 2025 Sohail Siadatnejad
// Exploring some certain GLSL-native style of coding
// arrowness-sdf
// todo: cleanup after explorations
// https://www.shadertoy.com/view/WfKfDh

// A few helpers to naturally emerge ...



// program state:
vec2 permaMouse; // fakes persistence

// used for highlighting debug
bool error = false;

const float PI2 = 2.0 * 3.1415926536;
const float INFTY = float(uint(-1));


vec2 rotate90(vec2 dir) {
    vec2 n = dir.yx * vec2(-1,1);
    // rotate dir by 90° to get normal = n
   return n;
}


// polygon-specific semantics

// Polygin-SDFizer

// Beautifully breaks down into per-point updates ( in fact, per-segment ).
// poly_sdf_state , PolySdfState
struct PolySdfState {
    float minV;
    float minE;
    int num_intersections;
};

PolySdfState init_state() {
    PolySdfState st;
    st.minV = INFTY;
    st.minE = INFTY;
    st.num_intersections = 0;
    return st;
}

// Core:

// The Update per-vertex of polygon.
// Updates for vertex. (or: segment)
// Assumes a closed polygon (and: CCW)
// Avoids holsing fixed large vector-array as "state'
// ⇒ Avoids `N`.
// Can be non-convex

// Some details are partially inspired by https://www.shadertoy.com/user/01000001
void update_pervertex(inout PolySdfState state, vec2 p, vec2 V0, vec2 V1) {

    float distV0 = length(p-V0);
    state.minV = min(state.minV, distV0);

    float len01 = length(V1-V0); // Line segment's length
    vec2 dir = (V1-V0) / len01; // Normalized segment dir
    float s01_ = dot(p-V0, dir);
    bool segmCyl_between = s01_ > 0. && s01_ < len01; 

    if (segmCyl_between) {
       vec2 n = rotate90(dir);
       float cylDist = abs(dot(p-V0, n));
       state.minE = min(state.minE, cylDist);
    }

    float fy = (p.y - V0.y) / (V1.y - V0.y);
    bool fy_within01 = fy > 0. && fy < 1.;

    if (fy_within01 ) {
        // fy is valid ⇒ point is not above or below line

        if ( p.x < min(V0.x, V1.x) ) { 
            // point is left of BB
            state.num_intersections++;

        } else if ( p.x > max(V0.x, V1.x) ) { 
           // do nothing! ( right side of BB )

        } else { // x within [min, max]
            // point is in BB

            // ⇒ check for intersection:

            if ( mix(V0.x, V1.x, fy) > p.x) {
                state.num_intersections++;
            }
        }
    }
}
// Historical cadence:
// void update_point(vec2 p, vec2 V0, vec2 V1, inout PolySdfState state) {
// void update_point(inout PolySdfState state, vec2 p, vec2 V0, vec2 V1) {
// void update_vertex(inout PolySdfState state, vec2 p, vec2 V0, vec2 V1) {
// renaming update_vertex -> update_pervertex

// The sign convention dependes on CW/CCW convention.
// CW => negate it.
float conclude_sdf(in PolySdfState state) {
    bool inside = (state.num_intersections % 2 == 1 );
    float _sign = ( inside?  -1. : +1. );
    return min(state.minE, state.minV) * _sign;
}



// Kadr: a 4th layer (third "statefulness layer") / separate spemtics
/// naturally emerged as separate (laminated!)
// together (wokds in viscosity(!) with) with ph_apply_transform()
// no longer "ph" semantics. 
// Kadr-semantics

// 4th state-layer
struct Kadr {
    mat2 ph_tr;
    vec2 ph_orig;
};

// Kadr ph_tr;
Kadr ph_kadr;

// Core of Kadr-semantics
vec2 apply_kadr(vec2 v) {
   // return ph_tr * (v - ph_orig);
   // return (ph_tr * v) + ph_orig;
   return (ph_kadr.ph_tr * v) + ph_kadr.ph_orig;
}
// `ph_apply_transform`() -> `apply_kadr`

void set_kadr(Kadr kadr) {
   ph_kadr = kadr;
}
// ph_set_frame() -> set_kadr

void kill_kadr() {
   // ph_tr=mat2(0,0,0,0);
   // ph_tr=mat2(INFTY,INFTY,INFTY,INFTY);
   // vec2 ph_orig = vec2(INFTY,INFTY);
   ph_kadr.ph_tr=mat2(INFTY,INFTY,INFTY,INFTY);
   ph_kadr.ph_orig=vec2(INFTY,INFTY);
}
// kill_kadr -> wipe, reset...



// PH-semantics: Polygon-handler
// a 3rd layer (third "statefulness layer")

// PH is a driver-only.
// state of the loop (polygon-wise)
// , not the core SDF emerging-herlpers.
vec2 firstV = vec2(INFTY,INFTY);
vec2 prevV = vec2(INFTY,INFTY);
// vec2 lastV;

/*
old design:
void update_last(inout PolySdfState state, vec2 p, vec2 v, bool first=false) {
   if (first) {
      firstV = v;
   }
   update_vertex(state, p, vert, other);
}
*/
// old: void update_last(inout PolySdfState state, vec2 p, vec2 v, bool first=false) {
// update_first
// update_next
// update_last update_back_to_first

void update_first(
      inout PolySdfState state,
      vec2 p,
      // mat2 transform, vec2 orig,
      // Kadr kadr,
      vec2 V0, vec2 V1
) {
   // ph_set_frame();
   /*
   ph_tr = transform;
   ph_orig = orig;
   */
   // ph_kadr = kadr;

   // V0 = ph_apply_transform(V0);
   // V1 = ph_apply_transform(V1);
   // firstV = V0;
   // lastV = V1;
   firstV = V0;
   update_pervertex(state, p, firstV, V1);
   prevV = V1;
}
void update_next(inout PolySdfState state, vec2 p, vec2 v) {
   // v = ph_apply_transform(v);
   update_pervertex(state, p, prevV, v);
   prevV = v;
}
void update_last(inout PolySdfState state, vec2 p, vec2 v) {
   // v = ph_apply_transform(v);

   update_pervertex(state, p, prevV, v);
   update_pervertex(state, p, v, firstV);
   // reset (for debugging)
   // firstV=vec2(-1,-1);
   firstV=vec2(INFTY, INFTY);
   prevV=vec2(INFTY, INFTY);
}


// The `vec2[N]`-specific semantics!
const int N = 7;
// from my interpreatin of sdPolygon() in https://www.shadertoy.com/view/wc3fzf
float polygon_sdf( in vec2 p, in vec2[N] v)
// float sdPolygon( in vec2 p, in vec2[N] v )
{
    PolySdfState state;
    state = init_state();
    
    for( int i = 0; i < N; i++ ) {
        vec2 V0 = v[i];
        vec2 V1 = v[ (i==N-1) ? 0 : i+1 ];

        update_pervertex(state, p, V0,V1);

    }
    return conclude_sdf(state);
}




// "box"-based semantics
// Requires a vec2[2]
//   e.g., vec2 box[] = vec2[2]

bool GT2d(vec2 p, vec2 a) {
   return p.x >= a.x && p.y >= a.y;
}

bool within01box(vec2 uv) {
    // return GT2d(uv, vec2(0,0)) && GT2d(uv*(-1.0), vec2(-1.0,-1.0));
    const vec2 ZERO = vec2(0,0);
    const vec2 ONE = vec2(1.0,1.0);
    // return GT2d(uv, ZERO) && GT2d(-uv, -ONE);
    // return GT2d(uv, ZERO) && GT2d(-uv+ONE, ONE-ONE);
    return GT2d(uv, ZERO) && GT2d( ONE-uv, ZERO);
    // or use min
}


// vec2 box[] = vec2[2](vec2(40,70), vec2(180,270));
vec2 box[] = vec2[2](vec2(40,70), vec2(180,470));


vec2 withinBoxUV(vec2 pixXY) {
   // float u = (PixXY.x - boxLU.x) /
   // vec2 wh = vec2(box[1] box[0]);
   vec2 wh = box[1] - box[0];
   // float dx = (PixXY.x - box[0].x);
   // float dy = (PixXY.y - box[0].y);
   vec2 dxy = pixXY - box[0];
   // float u = dx / wh.x;
   // float v = dy / wh.x;
   vec2 uv = dxy.xy / wh.xy;
   // return uv.x >= 0.0 && uv.x <= 1.0 && uv.y >= 0.0 && uv.y <= 1.0;
   // return uv.x >= 0.0 && uv.x * (-1.0) >= -1.0 && uv.y >= 0.0 && uv.y* (-1.0) >= -1.0;
   // return uv.x >= 0.0  && uv.y >= 0.0 && uv.x * (-1.0) >= -1.0 && uv.y* (-1.0) >= -1.0;
   // return uv >= vec2(0,0) && uv.x * (-1.0) >= -1.0 && uv.y* (-1.0) >= -1.0;
   // return GT2d(uv, vec2(0,0)) && GT2d(uv*(-1.0), vec2(-1.0,-1.0));
   // return GT2d(uv, vec2(0,0)) && GT2d(uv*(-1.0), vec2(-1.0,-1.0));
   // return within01box(uv);
   // separated withinBoxUV from inBox
   return uv;
}
float withinBoxSDF2(vec2 pixXY) {
   vec2 dxy0 = -(pixXY - box[0]);
   vec2 dxy1 = -(box[1] - pixXY);
   return min(min(
      min(dxy0.x, dxy0.y),
      min(dxy1.x, dxy1.y)),
      0.0);
}

// corner-based: defining a corner-based frame for such semantics
// for right-angle SDF semantics! (LURD Corners!)

// LURD corner-types
const int _LD = 0;
const int _RD = 1;
const int _RU = 2;
const int _LU = 3;
const ivec2 CORNERS[4] =  ivec2[4](
    // CCW
    ivec2(0, 0),
    ivec2(1, 0),
    ivec2(1, 1),
    ivec2(0, 1)
  );

// diagonal matrix!
vec2 cornerToDir(ivec2 corner) {
    float signX = corner[0] == 0 ? +1.0 : -1.0; // *2.0 - 1.0;
    float signY = corner[1] == 0 ? +1.0 : -1.0;
    return vec2(signX, signY);
}

// withinBoxXY -> transformIsoCorner -> cornerTransformIso
vec2 cornerTransformIso(vec2 pixXY, ivec2 corner) {
    // vec2 dxy = pixXY - box[corner];
    // sign0
    /*
    float signX = corner[0] == 0 ? +1.0 : -1.0; // *2.0 - 1.0;
    float signY = corner[1] == 0 ? +1.0 : -1.0;
    */
    vec2 signXY = cornerToDir(corner);
    float x = pixXY.x - box[corner[0]].x;
    float y = pixXY.y - box[corner[1]].y;
    // return vec2( x * signX, y * signY );
    // return vec2( x * signXY.x, y * signXY.y );
    return vec2(x,y) *  signXY;
}

bool inBox(vec2 pixXY) {
   vec2 uv = withinBoxUV(pixXY);
   return within01box(uv);
}

void minIt(inout float sdf, in float v) {
   sdf = min( sdf, v );
}
void maxIt(inout float sdf, in float v) {
   sdf = max( sdf, v );
}

// cornerSDF2
// cornerSDF1: deprecated

// outside is negative ⇒ min = interesectio?
// inside = (+) = is.
// float cornerSDF2(vec2 p, float R) {
float cornerSDF1(vec2 p, float R) {
   // xy: keep it iso-geom
   vec2 c = vec2(R,R);
   vec2 pc = p - c;
   // float rsdf = pc.x*pc.x + pc.y*pc.y - R*R;
   // return max(rsdf, max(pc.x, pc.y)); // union
   // return max(rsdf, max(pc.x-R, pc.y-R)); // union
   // return sqrt(rsdf);
   float rsdf = pc.x*pc.x + pc.y*pc.y;
   float circleSDF = -(sqrt(rsdf) - R); // (-) ⇒ outside
    // shifted
   float rightnessX = +(p.x-R); // (+) => right (right=inside it, left=outside)
   float upnessY = (p.y-R);
   float positvX = (p.x);
   float positvY = (p.y);
 
   float sdf = -INFTY;
   // float sdf = circleSDF;
   maxIt(sdf, circleSDF);
   maxIt(sdf, min(rightnessX, positvY));
   maxIt(sdf, min(upnessY, positvX));
   // minIt(sdf, min(positvY, positvX));
   return sdf;
   
}

// (+) = inside
float cornerSDF2(vec2 p, float R) {
   // method 2 : shrink and re-expand
   // R = 0.0;
   vec2 c = vec2(R,R);
   vec2 pc = p - c;
   float rsdf = pc.x*pc.x + pc.y*pc.y;
   float circleSDF = -sqrt(rsdf);
   // half_planes
   // vec2 hp[] = vec2[2](vec2(R, 0.0), vec2(0.0,R));
   vec2 hp[] = vec2[2](vec2(0,0), vec2(0,0));
   // vec2 hp[] = vec2[2](vec2(-R,0), vec2(0,-R));

   // and their normals (inward)
   vec2 hpn[] = vec2[2](vec2(+1.0, 0.0), vec2(0.0, +1.0));

   float sdf = +INFTY;
   int ctr = 0;
   for(int i = 0 ; i < 2; i++) {
      float hp_sdf  = dot( hpn[i] , (pc - hp[i]));
      vec2 perpCCW = vec2(-hpn[i].y, hpn[i].x);
      // vec2 perpCW = vec2(hpn[i].y, -hpn[i].x);
      float where = dot( perpCCW, pc - hp[i]);
      if ( (i == 0 && where > 0.0)
         || (i == 1 && where < 0.0) )
      {
         minIt(sdf, hp_sdf);
         ctr += 1;
      }
      else {
         //minIt(sdf, circleSDF);
      }
   }
   if (ctr == 0) {
      minIt(sdf, circleSDF);
   }
   // return max(sdf, -R+circleSDF)+R;
   // float vertSDF = +(R+circleSDF);
   float vertSDF = +(circleSDF);
   // maxIt(sdf, +vertSDF);
   // minIt(sdf, +vertSDF);
   // return sdf+R;
   // minIt(sdf, vertSDF);
   return sdf + R;
   
}
/*
float boxSDF2(vec2 p, float R) {
   // xy: keep it iso-geom
   vec2 pc = p - vec2(R,R);
   return pc.x*pc.x + pc.y*pc.y - R*R;
}
*/


// Mouse and driver semantics
// To avoid semantics leaking

void sortff(inout float x1, inout float x2 ) {
   if (x1 > x2) {
      float temp = x1;
      x1 = x2;
      x2 = temp;
   }
}
void apply_on_box(vec2 permaMouse) {

    // [1] is pivot
    if (GT2d(-permaMouse, -box[1]) ) {
       // permaMouse is smaller
       box[0] = permaMouse;
    } else {
       // permaMouse is larger
       vec2 temp = box[1];
       box[1] = permaMouse;
       box[0] = temp;
    }
    sortff(box[0].x, box[1].x);
    sortff(box[0].y, box[1].y);
}

vec3 zens(vec2 pixCoord) {
    
    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = pixCoord/iResolution.xy;
    
    // Time varying pixel color
    vec3 col3 = 0.5 + 0.5*cos(iTime+uv.xyx+vec3(0,2,4));
    return col3;

}

// Visulise SDF using contours

vec4 visualise(float sdf, vec2 pixCoord, bool error) {
    float shade = (sin(sdf/ 20.0 * PI2) * 0.5 + 0.5);
    vec4 outcol;
    // if (inb && minS >= 0.0) {
    if (sdf >= 0.0) {
       vec3 col3 = zens(pixCoord);
       outcol = vec4(col3, 1.0) * shade;
    } else {
       outcol = vec4(0,0,0,0) + shade * 0.4;
    }
    return outcol;

}
mat2 rotMat2(float rad)
{
    float s = sin(rad);
    float c = cos(rad);
    return mat2( c, -s, s, c );
}
void mainImage( out vec4 pixColor, in vec2 pixCoord )
{
    error = false;
    
    // shape
    if (iMouse.z > 0.0) {
       // box[0] = iMouse.xy;
       apply_on_box(iMouse.xy);
       // does not work:
       permaMouse = iMouse.xy;
    }
    // some shape geometry:
    // vec2 midP = (box[0] + box[1]) * 0.5;
    float midx = (box[0].x + box[1].x) * 0.5;
    vec2 midP = vec2(midx, INFTY);
    // float height = 50.0;
    float height = 80.0 * (1.0+sin(iTime / 5.0 * PI2));
    midP.y = box[1].y - height;
    

    bool inb = inBox(pixCoord);

    const float R = 30.0;
    float minS = INFTY;
    // minS = min(minS, withinBoxSDF2(pixCoord));
    for(int ci=0; ci < 4; ci++) {
        // round corder's SDF^2 value
        vec2 cornerity = cornerTransformIso(pixCoord, CORNERS[ci]);
        // vec2 cornerity = cornerToDir(CORNERS[ci]) * (pixCoord-;
        float rcs2 = cornerSDF2(cornerity, R);
        // minS = min(minS, rcs2);
        minIt(minS, rcs2);
    }
 
    // float cornerity_ = cornerSDF(cornerTransformIso(pixCoord-midP, CORNERS[_LU]), R);
    float cornerity_ = cornerSDF2(cornerToDir(CORNERS[_LU]) * (pixCoord-midP), R);
    // minS =  min( minS,  -cornerity_);
    // minS =  -cornerity_;
    minIt(minS, -cornerity_);



    // vec2 p = pixCoord / 40.0 - 4.0;
    vec2 p = pixCoord; // must keep it as it is in the main program.

    // a third layer! (of statefulness!) of polygon-handler
    // ph_set_frame(mat2(0,0,0,0));

    // vec2 lastV;
    PolySdfState ss2 = init_state();

    mat2 kadrm = -30.0*rotMat2(iTime * PI2 /5.0);
    vec2 arrow_centr = vec2(-5.66/2.0, (1.39+4.45)/2.0);
    vec2 arrow_pos = vec2(400, 300);
    vec2 kadr0 = -kadrm*arrow_centr + arrow_pos;
    Kadr kadr = Kadr(kadrm, kadr0);    
    set_kadr(kadr);

    // A CW shape (hence, negation in the end)
    // Note that we avoided `vec2[N]` as global or state
    // Why "store" it, if we don't need to? => Pure SDF helper
    
    // The deriver (PH) is not good: it uses "update" of variables (hence, potentially some pure compiler-optimisaiotn cannot be done)
    // But that is not the main purpose of this program.
    // My main purpose is `update_pervertex()`
    update_first(ss2, p,
                        apply_kadr(vec2(-2.62, 4.45)),
                        apply_kadr(vec2(-0.26, 2.93)));
    update_next(ss2, p, apply_kadr(vec2(-2.68, 1.39)));
    update_next(ss2, p, apply_kadr(vec2(-2.64, 2.39)));
    update_next(ss2, p, apply_kadr(vec2(-5.66, 2.41)));
    update_next(ss2, p, apply_kadr(vec2(-5.62, 3.53)));
    update_last(ss2, p, apply_kadr(vec2(-2.66, 3.49)));
    // Last update, applies two points ^. Symmetric-join-closure with first line does one, and gets two points.
    
    kill_kadr();

    float s2 = -conclude_sdf(ss2);
    // ^ negated because it is CW ^.

    maxIt(minS, s2);
 

    // Output to screen
    pixColor = visualise(minS, pixCoord, error);
}
