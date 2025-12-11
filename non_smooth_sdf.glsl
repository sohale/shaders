// forked from: https://www.shadertoy.com/view/3sKXDR - by koktszfung in 2019-11-15
/**
 *	ray marching tools adapted from hughsk's 2D SDF Toy https://www.shadertoy.com/view/XsyGRW
 */

#define TRACE_STEPS 20
#define TRACE_RAY

// 0 = Distance Field Display
// 1 = Raymarched Edges
// 2 = Resulting Solid
// 3 = Distance Field Polarity
#define DISPLAY 0

// 0 = Angle controlled By iTime
// 1 = Angle controlled By iMouse
#define MOUSE 0
const float PI = 3.14159265359;



//==========================================



//////// Math-y part
//// Part 1: common

const int K = 2;
const int M = 6;   // (K+1)(K+2)/2 = 6
const int N = 9;

// Multi-indices α = (α1, α2) in some enumeration.
const ivec2 alphaList[M] = ivec2[M](
    ivec2(0,0),
    ivec2(1,0), ivec2(0,1),
    ivec2(2,0), ivec2(1,1), ivec2(0,2)
);

float BadNaN = -10000000.0;
// float if_val(i,i0,float val)
float cond_c(int idx, float arg1) {
    // return (idx==0?1.0:(idx==1?1.0:arg1 * (arg1-1.)));
    /*
    return
       (idx==0) ? // 0
            1.0
       :(idx==1?  // 1
            1.0
       :          // 2
            arg1 * (arg1-1.)
       );
    */
    return
       (idx==0) ? // 0
            1.0
       :(idx==1?  // 1
            1.0
       :(idx==2?  // 2
            arg1 * (arg1-1.)
       :(idx==3?  // 3
            arg1 * (arg1-1.) * (arg1-2.)
       :(idx==4?  // 4
            arg1 * (arg1-1.) * (arg1-2.) * (arg1-3.)
       :(idx==5?  // 5
            arg1 * (arg1-1.) * (arg1-2.) * (arg1-3.) * (arg1-4.)
       :          // 6+
            // 0.0f/0.0f // NaN!
            // sqrt(-0.0-1.0)
            // -10000000.0
            BadNaN
       )))));
}

// α! = α1! α2!
float alphaFactorial(ivec2 a) {
    float ax = float(a.x);
    float ay = float(a.y);
    // float ax_1 = float(a.x-1);
    // float ay_1 = float(a.x-1);
    //
    //float f1 = (a.x==0?1.0:(a.x==1?1.0:ax * ax_1)); // handles 0,1,2 only
    //float f2 = (a.y==0?1.0:(a.y==1?1.0:ay * ay_1));

    float f1 = cond_c(a.x, ax);
    float f2 = cond_c(a.y, ay);

    return f1 * f2;
}

// Part 2: Inference: approximate Taylor coeffiecnts

float fast_power_n(float x, int p) {
    /*
    float term = 1.0;
    if(p == 1) term = x;
    //else if(p == 2) term *= x*x;
    else if(p == 2) term = x*x;
    else if(p == 3) term = x*x*x;
    else if(p == 4) term = x*x;
    else if(p == 5) term = x*x;
    else if(p == 6) term = BadNaN;
    */
    if(p==0) return 1.0;
    else if(p == 1) return x;
    else if(p == 2) return x*x;
    else if(p == 3) return x*x*x;
    else if(p == 4) return x*x*x*x;
    else if(p == 5) return x*x*x*x*x;
    else if(p == 6) return BadNaN;
}
// ===============================================================
// Compute φ_α(Δx) = Δx^α / α!
// ===============================================================
float phi_alpha(vec2 dx, ivec2 a) {
    float term = 1.0;
    // dx^α = dx.x^a.x * dx.y^a.y
    /*
    if(a.x == 1) term *= dx.x;
    //else if(a.x == 2) term *= dx.x*dx.x;
    else if(a.x == 2) term *= dx.x*dx.x;
    else if(a.x == 3) term *= dx.x*dx.x*dx.x;
    else if(a.x == 4) term *= dx.x*dx.x;
    else if(a.x == 5) term *= dx.x*dx.x;
    else if(a.x == 6) term = BadNaN;
    */
    term *= fast_power_n(dx.x, a.x);

    /*
    if(a.y == 1) term *= dx.y;
    else if(a.y == 2) term *= dx.y*dx.y;
    */
    term *= fast_power_n(dx.y, a.y);

    return term / alphaFactorial(a);
}


// side inputs

// buffers for temp work
// float Awork[M][M];
float Awork_AA[M*M];
float bwork[M];

//    float ATA[M][M];
float ATA_AA[M*M];
float ATF[M];

void setAA_ATA(int i, int j, float val) {
    // void set_ATA_AA(int i, int j, float val) {
    //void set_AA(ATA_AA, i,j, 0.0) {
    ATA_AA[i+j*M] = val;
}
float getAA_ATA(int i, int j) {
    return ATA_AA[i+j*M];
}
void addAA_ATA(int i, int j, float val) {
    ATA_AA[i+j*M] += val;
}
// Awork_AA:
void setAA_Awork(int i, int j, float val) {
    Awork_AA[i+j*M] = val;
}
float getAA_Awork(int i, int j) {
    return Awork_AA[i+j*M];
}
void addAA_Awork(int i, int j, float val) {
    Awork_AA[i+j*M] += val;
}


// ===============================================================
// Compute Taylor parameters a[m] for m = 0..M-1
// using normal equations: a = inverse(A^T A) * (A^T f)
// ===============================================================
void computeTaylorCoefficients(
    // in vec2 x0,
    in int N,
    // in vec2 X[N/*N*/],
    in vec2 dx_buffer[N],
    in float F[N/*N*/],
    out float a[M]          // output coefficients
){
    // Build ATA = A^T A (M x M)
    float ATA_AA[M*M];
    float ATF[M];

    // init
    for(int i=0;i<M;i++) {
        ATF[i] = 0.0;
        for(int j=0;j<M;j++)
            // ATA[i][j] = 0.0;
            // set_AA(ATA_AA, i,j, 0.0);
            setAA_ATA(i,j, 0.0);
    }

    // Accumulate A^T A and A^T f
    for(int k=0;k<N;k++){
        // vec2 dx = X[k] - x0;
        vec2 dx = dx_buffer[k];

        // Compute φ_m(dx)
        float phiRow[M];
        for(int m=0;m<M;m++){
            phiRow[m] = phi_alpha(dx, alphaList[m]);
        }

        // Update ATA and ATF
        for(int i=0;i<M;i++){
            ATF[i] += phiRow[i] * F[k];
            for(int j=0;j<M;j++){
                // ATA[i][j] += phiRow[i] * phiRow[j];
                addAA_ATA(i,j , phiRow[i] * phiRow[j]);
            }
        }
    }

    // Solve ATA * a = ATF via naive Gaussian elimination.
    // (M is small, OK for GLSL; replace with a precomputed inverse if fixed.)
    // float Awork[M][M];
    // float bwork[M];

    for(int i=0;i<M;i++){
        bwork[i] = ATF[i];
        for(int j=0;j<M;j++)
            // Awork[i][j] = ATA[i][j];
            setAA_Awork(i,j, getAA_ATA(i,j));
    }

    // Gaussian elimination
    for(int i=0;i<M;i++){
        // pivot normalisation
        // float diag = Awork[i][i];
        float diag = getAA_Awork(i,i);
        for(int j=0;j<M;j++) { 
           // getAA_Awork(i,j) /= diag;
           setAA_Awork(i,j , getAA_Awork(i,j) / diag );
        }
        bwork[i] /= diag;

        // eliminate below and above
        for(int k=0;k<M;k++){
            if(k==i) continue;
            // float factor = Awork[k][i];
            float factor = getAA_Awork(k,i);
            for(int j=0;j<M;j++) {
                // Awork[k][j] -= factor * Awork[i][j];
                setAA_Awork(k,j, getAA_Awork(k,j) - factor * getAA_Awork(i,j));
            }
            bwork[k] -= factor * bwork[i];
        }
    }

    // bwork now holds a[]
    for(int i=0;i<M;i++)
        a[i] = bwork[i];
}




// Part 3: Evaluation

// ===============================================================
// Evaluate Taylor polynomial T(x) of order K in 2D
// ===============================================================
float evaluateTaylor(
    // in vec2 x,
    // in vec2 x0,
    in vec2 dx,
    in float a[M]      // coefficients
){
    // vec2 dx = x - x0;
    float sum = 0.0;

    // M = Taylor order?

    for(int m=0;m<M;m++){
        sum += a[m] * phi_alpha(dx, alphaList[m]);
    }

    return sum;
}

//

/* Usage:
float coeffs[M];
float SDF_buffer[N];
vec2 x_buffer[N]

vec2 dx = x - x0;
vec2 dx_buffer[N]; // = x_buffer[k] - x0;

computeTaylorCoefficients( N, dx_buffer, SDF_buffer, coeffs );
float err = evaluateTaylor(dx, coeffs);
*/
//==========================================




float smin(float a, float b, float k) {
  float h = clamp(0.5 + 0.5 * (b - a) / k, 0.0, 1.0);
  return mix(b, a, h) - k * h * (1.0 - h);
}

vec2 squareFrame(vec2 screenSize, vec2 coord) {
  vec2 position = 2.0 * (coord.xy / screenSize.xy) - 1.0;
  position.x *= screenSize.x / screenSize.y;
  return position;
}

vec2 pmod(vec2 p, float m) {
	return mod(p + m*.5, m) - m*.5;
}

float pmod(float p, float m) {
	return mod(p + m*.5, m) - m*.5;
}

vec2 rotate(vec2 p, float a){
 	return vec2(p.x*cos(a) - p.y*sin(a),
                p.x*sin(a) + p.y*cos(a));   
}



float shape_circle(vec2 p, float r) {
  return length(p) - r;
}

float shape_rect(vec2 p, vec2 b) {
  vec2 d = abs(p) - b;
  return length(max(d, 0.0)) +  min(max(d.x, d.y), 0.0);  // out + in
}

float shape_line(vec2 p, vec2 a, vec2 b) {
  vec2 dir = b - a;
  return abs(dot(normalize(vec2(dir.y, -dir.x)), a - p));
}

float shape_segment(vec2 p, vec2 a, vec2 b) {
  float d = shape_line(p, a, b);
  float d0 = dot(p - b, b - a);
  float d1 = dot(p - a, b - a);
  return d1 < 0.0 ? length(a - p) : d0 > 0.0 ? length(b - p) : d;
}

float shape_gear(vec2 p, float n, float lr, float hr, float oa) {  // point, number, low radius, high radius, teeth's offset angle
    float pr = length(p);	
    float pa = atan(p.y, p.x) + PI*2.;
    
    float ma = 2.*PI/n;  // ma: modulus angle
    pa = mod(pa, ma);  // new pa: remainder angle
    p = pr * vec2(cos(pa), sin(pa));  // take values in the first modulus angle
    
    // circle is a discont' part of the primitive circle shape
    float gr = (pa > oa && pa < ma-oa)? hr : lr;  // gear's radius at angle pa
    float circle = shape_circle(p, gr);
    
    // segment is a cont' part of the primitive segment shape
    float segment1 = shape_segment(p, vec2(lr*cos(oa), lr*sin(oa)), vec2(hr*cos(oa), hr*sin(oa)));  // start angle   
    float segment2 = shape_segment(p, vec2(lr*cos(ma-oa), lr*sin(ma-oa)), vec2(hr*cos(ma-oa), hr*sin(ma-oa)));  // end angle    
    float segment = min(segment1, segment2);
    
    // when p is inside the shape, the value is negative
    if (pr < gr) {
        return -min(-circle, segment);
    } else {
    	return min(circle, segment);
    }
}


/*
float master_sdf(vec2 uv) {
   return SAMPLER(uv,t2);
}
*/

float SAMPLER(vec2 p, float t2) {
    #if MOUSE
    	float a1 = iMouse.x / iResolution.x * PI*2.;
    	float a2 = iMouse.x / iResolution.x * -PI*2. + PI*.5;
    #endif
    #if !MOUSE
        float T = 3.;  // period
        float a1 = smoothstep(0., T*.75, mod(t2, T)) * PI*2.;  // stop for T*.25s
        float a2 = smoothstep(0., T*.75, mod(t2, T)) * -PI*2. + PI*.5;  // stop for T*.25s, +PI/2 phase
    #endif

    
    vec2 p1 = rotate(p - vec2(-.6, 0.), a2);  // left
    vec2 p2 = rotate(p - vec2(.6, 0.), a1);  // right
    
    float gear1 = shape_gear(p1, 3., .5, .7, .83);  // left
    float gear2 = shape_gear(p2, 3., .5, .7, .83);  // right

    return min(gear1, gear2);
}



vec3 draw_outline(float d, float thickness) {
  const float aa = 3.0;
  return vec3(smoothstep(0.0, aa / iResolution.y, max(0.0, abs(d) - thickness)));
}

vec3 draw_outline(float d) {
  return draw_outline(d, 0.0025);
}

float draw_solid(float d) {
  return smoothstep(0.0, 3.0 / iResolution.y, max(0.0, d));
}

vec3 draw_polarity(float d, vec2 p, float t3) {
  p += t3 * -0.1 * sign(d) * vec2(0, 1);
  p = mod(p + 0.06125, 0.125) - 0.06125;
  float s = sign(d) * 0.5 + 0.5;
  float base = draw_solid(d);
  float neg = shape_rect(p, vec2(0.045, 0.0085) * 0.5);
  float pos = shape_rect(p, vec2(0.0085, 0.045) * 0.5);
  pos = min(pos, neg);
  float pol = mix(neg, pos, s);

  float amp = abs(base - draw_solid(pol)) - 0.9 * s;

  return vec3(1.0 - amp);
}

vec3 draw_distance(float d, vec2 p) {
  float t = clamp(d * 0.85, 0.0, 1.0);
  vec3 grad = mix(vec3(1, 0.8, 0.5), vec3(0.3, 0.8, 1), t);

  float d0 = abs(1.0 - draw_outline(mod(d + 0.1, 0.2) - 0.1).x);
  float d1 = abs(1.0 - draw_outline(mod(d + 0.025, 0.05) - 0.025).x);
  float d2 = abs(1.0 - draw_outline(d).x);
  vec3 rim = vec3(max(d2 * 0.85, max(d0 * 0.25, d1 * 0.06125)));

  grad -= rim;
  grad -= mix(vec3(0.05, 0.35, 0.35), vec3(0.0), draw_solid(d));

  return grad;
}

vec3 draw_trace(float d, vec2 p, vec2 ro, vec2 rd, float t4) {
  vec3 col = vec3(0);
  vec3 line = vec3(1, 1, 1);
  vec2 _ro = ro;

  for (int i = 0; i < TRACE_STEPS; i++) {
    float t = SAMPLER(ro, t4);
    col += 0.8 * line * (1.0 - draw_outline(length(p.xy - ro) - abs(t), 0.));
    col += 0.2 * line * (1.0 - draw_solid(length(p.xy - ro) - abs(t) + 0.02));
    col += line * (1.0 - draw_solid(length(p.xy - ro) - 0.015));
    ro += rd * t;
    if (t < 0.01) break;
  }

  #ifdef TRACE_RAY
    col += 1.0 - line * draw_outline(shape_segment(p, _ro, ro), 0.);
  #endif

  return col;
}

float indotness(vec2 uv, vec2 dotuv, float r) {
   return draw_solid(length(uv - dotuv) - r);
}
float plus(vec2 uv, vec2 dotuv, float r) {
   return draw_solid(length(uv - dotuv) - r);
}

// gradient...
// this attempt led to the math part. see above.


float evaluate_smothness(float delta, vec2 x0, float t1) {
    /*
         delta: radius of neighbourhood samples. deafult: 0.01
    */

    float coeffs[M];
    float SDF_buffer[N];
    // vec2 dx === uv === x0;
    vec2 dx_buffer[N]; // = x_buffer[k] - x0;
    int GM = 3;
    {
    int ctr=0;
    for(int i=0;i<GM;i++) {
       for(int j=0;j<GM;j++) {
          dx_buffer[ctr].x = float(i-1) * delta;
          dx_buffer[ctr].y = float(j-1) * delta;
          ctr += 1;
       }
    }
    }
    // ask compiler to guarantee ctr<GM*GM
    
    //runtime assert:
    // N == GM*GM;

    for(int i=0;i<N;i++) {
       SDF_buffer[i] = SAMPLER(dx_buffer[i] + x0, t1);
    }
    
    // stores into: coeffs
    computeTaylorCoefficients( N, dx_buffer, SDF_buffer, coeffs );
    vec2 dx0 = vec2(0.); // taylor around 0 (in dx_buffer coords)
    float sdf_approx = evaluateTaylor(dx0, coeffs);
    return sdf_approx;
}

// see master_sdf()
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  float t1 = iTime * 0.0;
  float t = iTime * 0.5;
  vec2 uv = squareFrame(iResolution.xy, fragCoord);
  float d;
  vec3 col;
  vec2 ro = vec2(iMouse.xy / iResolution.xy) * 2.0 - 1.0;
  ro.x *= squareFrame(iResolution.xy, iResolution.xy).x;



  vec2 ref0 = vec2(0.2,0.3);

  vec2 rd = normalize(ref0-ro);

  d = SAMPLER(uv, t1);

  #if DISPLAY == 0
    col = vec3(draw_distance(d, uv.xy));
    #if MOUSE == 0
    	col -= (iMouse.z > 0.0 ? 1.0 : 0.0) * vec3(draw_trace(d, uv.xy, ro, rd, t1));
    #endif
  #endif
  #if DISPLAY == 1
    col = vec3(0) + 1.0 - vec3(draw_outline(d));
    #if MOUSE == 0
    	col += (iMouse.z > 0.0 ? 1.0 : 0.0) * vec3(1, 0.25, 0) * vec3(draw_trace(d, uv.xy, ro, rd, t1));
    #endif
    col = 1. - col;
  #endif
  #if DISPLAY == 2
    col = vec3(draw_solid(d));
  #endif
  #if DISPLAY == 3
    col = vec3(draw_polarity(d, uv.xy, t1));
  #endif

  float indot = indotness(uv, ro, 0.015*3.0);
  col *= 1.0-(1.0 - indot);

  col *= indotness(uv, ref0, 0.015*3.0);

  float err =  evaluate_smothness(0.01, uv.xy, t1) - SAMPLER(uv.xy, t1);
  
  // col = vec3(0.) + abs(err * 100000.00);
  // col = vec3(0.) + abs(err * 5.00);
  // col = vec3(0.) + abs(err * 500.00);
  // col = vec3(1.) - abs(err * 500.00);
  // float non_smoothness = abs(err * 500.00);
  // col = col * (vec3(1.) - non_smoothness);
  float non_smoothness = (err * 50000.00);
  // col = col * (vec3(1.) + non_smoothness*vec3(1.,0.,0.));
  vec3 red_mark = vec3(1.,0.,0.) * non_smoothness;
  // treat as alpha, using mix()
  float alpha_ = clamp(abs(non_smoothness), 0., 1.);
  col = mix(col, red_mark, alpha_);


  fragColor.rgb = col;
  fragColor.a   = 1.0;
}
