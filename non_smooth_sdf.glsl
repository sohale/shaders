
/**
 * Uses C_inf to measure non-smoothness.
 * The SDF and the Raymarching-visualiser, from:
 * A fork of https://www.shadertoy.com/view/3sKXDR - by koktszfung in 2019-11-15, which in turn, had its ray marching tools adapted from hughsk's 2D SDF Toy https://www.shadertoy.com/view/XsyGRW
 */

#define TRACE_STEPS 20
#define TRACE_RAY

// 0 = Distance Field Display
// 1 = Raymarched Edges
// 2 = Resulting Solid
// 3 = Distance Field Polarity
#define STYLE_SDF 0
#define STYLE_SOLID_EDGES 1
#define STYLE_SOLID 2
#define STYLE_POLARITY 3

#define DISPLAY_STYLE  STYLE_SDF

// 0 = Angle controlled By iTime
// 1 = Angle controlled By iMouse
#define _BY_TIME_ANIMATIONSOURCE 0
#define _BY_MOUSE_ANIMATIONSOURCE 1

#define ANIMATION_SOURCE _BY_TIME_ANIMATIONSOURCE



#define TAYLOR_BASES_MODELFAMILY 0
#define RBF_BASES_MODELFAMILY 1

#define MODELFAMILY  RBF_BASES_MODELFAMILY



//////// Math-y part

const float PI = 3.14159265359;


float BadNaN = -10000000.0;
int BadNaNInt = 100000;

//==========================================


#if MODELFAMILY == TAYLOR_BASES_MODELFAMILY


//// Part 1: common, foundations, basic vectors, factorials.

/*
Multi-dimensional Taylor series:
Has two paremters: (D, M):
  D : Number of dimenions, 
  M: to what order we want to keep trck of the smooth/analytical approximation/local model
      C^∞, but we use up to C^M.
So,
  x ∈ ℝ^D,
  f ∈ C^M.

K = number of Taylor-series terms (if matrix form).
Hene, maximum power that iappears is M (m).
power = m 
A_{k,m}, power: m
Ak,m =α(m)!(Δxk)α(m)
a=(A^T A)^{−1}A^Tf


D = dimensions: D=2

M = (number of coefficients)
   M scalars in the vector of parameters
   e.g., for D=2,
   M = (K+1)(K+2)/2 = 6

In this reperesentation of Taylor
the terms are sorted in terms of `α`s:
Each α is a vector (of dim D) of integer indices.
The sum of those integers, are the order of the term (in classical sense)
α is a "A multi–index".


Hence, it's "order", the "multi-order" or each term:
    ∣α∣ = Order:
    ∣𝛼∣ = 𝛼1 + ⋯ + 𝛼𝐷
    ∣α∣ = α1 + ⋯ + α D.
    ∣α∣ =∑^D α_i

The extended factorial of (α) is:
   𝛼! = 𝛼_1! × ⋯ × 𝛼_𝐷!
   (α)! = (α_1)! × ⋯ × (α_D)!
   = ∏^D α_i


"Monomial" power: x^α
    x^α = ∏^D pow(x, α_i)

So, we extended sum, factorial, and power to α.


Taylor-like basis function:
   𝜙_𝛼 (𝑥) = 𝑥^𝛼 / 𝛼!
   ϕ_α (x) = x^α / α!
combines extended multi-power and multi-factorial 
(extended-power and extended factrorial).

*/

// const int K = 2;
const int M = 6;   // (K+1)(K+2)/2 = 6
const int N = 9;


#endif  // TAYLOR_BASES_MODELFAMILY



#if MODELFAMILY == RBF_BASES_MODELFAMILY
// RBF Territory:

// Number of radial basis coefficients: c0, c1*r, c2*r*r
const int D = 2;  // const-snobbery
const int M = D + 1;     // = 3 when D=2
// In taylor: was: M = 6;
const int N = 9;


#endif  // RBF_BASES_MODELFAMILY



// ======== COMMON AREA ( Gauss & Jordan) =========
// Taylor's version:

// side inputs

// buffers for temp work
// Used for/in/by Gauss_Jordan_elimination() 
// AA means array of array, ie, 2d-array
float Awork_AA[M*M];
float bwork[M];

float ATA_AA[M*M];
float ATF[M];

void setAA_ATA(int i, int j, float val) {
    ATA_AA[i+j*M] = val;
}
float getAA_ATA(int i, int j) {
    return ATA_AA[i+j*M];
}
void addAA_ATA(int i, int j, float val) {
    ATA_AA[i+j*M] += val;
}
void setAA_Awork(int i, int j, float val) {
    Awork_AA[i+j*M] = val;
}
float getAA_Awork(int i, int j) {
    return Awork_AA[i+j*M];
}
void addAA_Awork(int i, int j, float val) {
    Awork_AA[i+j*M] += val;
}
//  RBF's version of ^ this part was repeated by RBF, and got merged ^.


/*
⭐️⭐️⭐️ Gauss–Jordan elimination ⭐️⭐️⭐️

Solves:
   ✨ x = (A^T A)^{−1} A^T f
Solve ATA * x = ATF via naive Gaussian elimination.
x = params (a in tylor, c in RBF)
x[M]

GLSL notes:
  (M is small, OK for GLSL; replace with a precomputed inverse if fixed.)

Inputs are implicit ("side-band"): `ATA_AA`,`ATF`.
Some temp space is used: `bwork`, `Awork`.
*/
// Got merged! (Was) good food for refactoring ( RBF's version was merged).
void Gauss_Jordan_elimination(/*in float ATA[M*M], in float ATF[M],*/  out float x[M]) {

    // float Awork[M][M]; float bwork[M];
    // Make a copy: Awork ← ATA, bwork ← ATF
    for(int i=0; i<M; i++){
        bwork[i] = ATF[i];
    }
    for(int i=0;i<M;i++){
        for(int j=0;j<M;j++)
            // Awork[i][j] ← ATA[i][j];
            setAA_Awork(i,j, getAA_ATA(i,j));
            // perturbation:
            // if (i==0 && j==0)  addAA_Awork(i, j, 0.001);
    }

    // --- Gauss–Jordan elimination ---
    for(int i=0;i<M;i++){
        // Normalize pivot row
        // diag := Awork[i][i];
        float diag = getAA_Awork(i,i);
        for(int j=0;j<M;j++) { 
           // getAA_Awork(i,j) /= diag;
           setAA_Awork(i,j , getAA_Awork(i,j) / diag );
        }
        bwork[i] /= diag;

        // Eliminate other rows below and above
        for(int k=0;k<M;k++){
            if(k==i) continue;
            float factor = getAA_Awork(k,i); // / Awork[k][i]
            // Awork[k][:] = row(k) -= f * row(i)
            for(int j=0;j<M;j++) {
                // Awork[k][j] -= factor * Awork[i][j];
                float vij = getAA_Awork(k,j) - factor * getAA_Awork(i,j);
                setAA_Awork(k,j, vij);
            }
            bwork[k] -= factor * bwork[i];
        }
    }

    // `bwork` now holds x[M] = model[M]
    // x as params: (e.g.,  a[], c[])
    // --- Copy solution  back to x ---
    for(int i=0; i<M; i++)
        x[i] = bwork[i];
}


// ===== end of COMMON AREA =====


// ======== Taylor AREA: Brook Taylor © 1715 ! =========

#if MODELFAMILY == TAYLOR_BASES_MODELFAMILY


// Multi-indices α = (α1, α2) in some enumeration.
// We have M numebr of α's.
// α are the ...
// see alphaList_with_fac[M] instead!
// deprecated!
const ivec2 alphaList[M] = ivec2[M](
    ivec2(0,0),
    ivec2(1,0), ivec2(0,1),
    ivec2(2,0), ivec2(1,1), ivec2(0,2)
);



const int[13] FACTORIALS_TABLE = int[13](
    1, // 0!,
    1,
    2,
    6,
    24,  // 4!
    120,
    720,
    5040,
    40320,
    362880,
    3628800, // 10!
    39916800,
    479001600 // 12!
    // aready overflow!!
    // 6227020800 // 13!
    // 87178291200, // 14!
);
const int[13] FAC = FACTORIALS_TABLE;

// inline
int fast_factorial(int i) {
    return FACTORIALS_TABLE[i];
}

// this is int-only!
// thm: if idx == arg1, this is essentially, exactly, the factorial.
int cond_c(int idx, int arg1) {
    // ?float if_val(i,i0,float val)

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
            1
       :(idx==1?  // 1
            1
       :(idx==2?  // 2
            arg1 * (arg1-1)
       :(idx==3?  // 3
            arg1 * (arg1-1) * (arg1-2)
       :(idx==4?  // 4
            arg1 * (arg1-1) * (arg1-2) * (arg1-3)
       :(idx==5?  // 5
            arg1 * (arg1-1) * (arg1-2) * (arg1-3) * (arg1-4)
       :          // 6+
            // 0.0f/0.0f // NaN!
            // sqrt(-0.0-1.0)
            // -10000000.0
            BadNaNInt
       )))));
}


// α! = α1! α2!
// deprecated
// this is int-only too!!
// specicic for D=2: for 2D
int alphaFactorial(ivec2 a) {
    // float ax = float(a.x);
    // float ay = float(a.y);
    // float ax_1 = float(a.x-1);
    // float ay_1 = float(a.x-1);
    //
    //float f1 = (a.x==0?1.0:(a.x==1?1.0:ax * ax_1)); // handles 0,1,2 only
    //float f2 = (a.y==0?1.0:(a.y==1?1.0:ay * ay_1));

    /*
    // it is int !
    float f1 = float(cond_c(a.x, a.x));
    float f2 = float(cond_c(a.y, a.y));

    return f1 * f2;
    */
    int f1 = cond_c(a.x, a.x);
    int f2 = cond_c(a.y, a.y);
    // cond_c(a,a) = a ! , and we only used is thi way: so:
    // hance, use alphaFactorial_fast

    return f1 * f2;
}


// see above!
// α! = α1! α2!
// specicic for D=2
// this is int-only too.
// We could make it even faster: by look-up table of two elemnts (with row-size being a power-of-2, to shift.
// or even: combine it with alphaList already!
//    its only use it:   phi_alpha(dx, alphaList[m]);
// see alphaList_with_factorial
int alphaFactorial_fast(ivec2 a) {
    // return fast_factorial(a.x) * fast_factorial(a.y);
    return FACTORIALS_TABLE[a.x] * FACTORIALS_TABLE[a.y];
}

/*
// if it was D=3, etc
int alphaFactorial_fast_d3(ivec3 a) {
    return FACTORIALS_TABLE(a.x) *
      FACTORIALS_TABLE(a.y) *
      FACTORIALS_TABLE(a.z);
}
*/

// precomputing for: phi_alpha(dx, alphaList[m]);
// shortcut
int sct(int a0,int a1) {
    return alphaFactorial_fast(ivec2(a0,a1));
    // return FACTORIALS_TABLE[a0] * FACTORIALS_TABLE[a1];
}
// alphaList_with_factorial
/*
const int alphaList_with_fac[M] = int[M](
    sct(0,0),
    sct(1,0), sct(0,1),
    sct(2,0), sct(1,1), sct(0,2)
);
*/
// ugly and fast.
// combining table of factorial with table of alphas.
const int alphaList_with_fac[M] = int[M](
    FAC[0]*FAC[0],
    FAC[1]*FAC[0], FAC[0]*FAC[1],
    FAC[2]*FAC[0], FAC[1]*FAC[1], FAC[0]*FAC[2]
);
// todo: same for fast alpha powers



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
    else if(p == 4) {
        // return x*x*x*x;
        float x2 = x*x;
        return x2*x2;
    }
    else if(p == 5) {
        //return x*x*x*x*x;
        float x2 = x*x;
        return x2*x2*x;
    }
    else if(p == 6) return BadNaN;
}

/*
float phi_alpha_denom_m(m) {

   //int denom = alphaFactorial_fast(a);
   int denom = alphaFactorial_fast(a);
   
   return denoom;

}
*/

// skips alpha
// directly the: phi_alpha(dx, alphaList[m]);
int phi_alpha_denom_directly_from_m(int mi) {

    // return term / float(alphaFactorial(a));
    // return term / float(alphaFactorial_fast(a));
    // int denom = alphaFactorial_fast(a);
    // return term / float(denom);

    //int denom = alphaFactorial_fast(a);
    //int denom = alphaFactorial_fast(a);

    // return denoom;
    
    // directly: phi_alpha(dx, alphaList[m]);
    // int denom = alphaFactorial_fast(alphaList[m]);
    // denom +=1;
    // denom = -1;
    // denom = 10; //why it is not sensitivr to this?!
    
    
    // expanding alphaFactorial_fast():
    // ivec2 a = alphaList[m];
    // int denom = FACTORIALS_TABLE[a.x] * FACTORIALS_TABLE[a.y];
    int denom = alphaList_with_fac[mi];
    return denom;
    
}

// ===============================================================
// Compute φ_α(Δx) = Δx^α / α!
// ===============================================================
float phi_alpha_numerat(vec2 dx, ivec2 a) {
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

    /*
    moved into phi_alpha_denom_m
    // return term / float(alphaFactorial(a));
    // return term / float(alphaFactorial_fast(a));
    int denom = alphaFactorial_fast(a);
    return term / float(denom);
    */
    // int denom = phi_alpha_denom_directly_from_m(a, m);
    // skip the alpha part:
    // int denom = phi_alpha_denom_directly_from_m(a, m);
    // return term / float(denom);
    return term;
}

// give me (m,) too, and I make it faster for you.
float phi_alpha_(vec2 dx, ivec2 a, int m) {
    // float phi_alpha_numerat(vec2 dx, ivec2 a) {
    float numerat = phi_alpha_numerat(dx, a);
    /*
    int denom = alphaFactorial_fast(a);
    return numerat / float(denom);
    */
    
    // now: skip the alpha part:
    // int denom = phi_alpha_denom_directly_from_m(a, m);
    // skips ther a-dependence, but adds m, which is the source of `a`:
    int denom = phi_alpha_denom_directly_from_m(m);
    return numerat / float(denom);
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
    // float ATA_AA[M*M]; float ATF[M];

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
            phiRow[m] = phi_alpha_(dx, alphaList[m], m);
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

    // Gaussian elimination
    Gauss_Jordan_elimination(/*ATA_AA, ATF,*/ a); // temp: Awork,bwork
}




// Part 3: Evaluation

// ===============================================================
// Evaluate Taylor polynomial T(x) of order K in 2D
//    with M scalars in the vector of parameters
// ===============================================================
float evaluateTaylor_(
    // in vec2 x,
    // in vec2 x0,
    in vec2 dx,
    in float a[M]      // coefficients
){
    // vec2 dx = x - x0;
    float sum = 0.0;

    // M = Taylor order?

    for(int m=0;m<M;m++){
        sum += a[m] * phi_alpha_(dx, alphaList[m], m);
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

// model() points -> params // (fit)
// demodel() params -> eval //(approx) // local-model
// evaluate model
// demodelle
float demodel( in vec2 dx, in float model[M]) {
   return evaluateTaylor_(dx, model);
}
// to model (verb)
// stores into: coeffs
// model (n.) = output coefficients
// is near zero (civinity surrounding, neighbourhood (0)
// void model(
void do_model(
    in vec2 dx_buffer[N], in float F[N],
    out float model[M]
) {
   computeTaylorCoefficients(N,  dx_buffer, F, model);
}

#endif  // MODELFAMILY == TAYLOR_BASES_MODELFAMILY


//==========================================



//==========================================
// RBF Territory:

#if MODELFAMILY == RBF_BASES_MODELFAMILY


// Number of radial basis coefficients: c0, c1*r, c2*r*r
// D = 2, M = D + 1;


// Fit f ≈ c0 + c1*r + c2*r*r  (K=2)
void fitRadialModel(
    vec2 dx[N],      // input offsets
    float F[N],      // SDF samples
    out float c[M]   // c0, c1, c2
){
    // zero
    // --- Clear ATA and ATF ---
    for (int i=0; i<M; i++) {
        ATF[i] = 0.0;
        for (int j=0; j<M; j++) {
            setAA_ATA(i, j, 0.0);
        }
    }

    // accumulate
    // --- Accumulate ATA and ATF ---
    for (int k=0; k<N; k++) {
        float r = length(dx[k]);
        // radial basis row: 1, r, r^2, ..., r^(R-1)
        /*
        float row0 = 1.0;
        float row1 = r;
        float row2 = r*r;
        float row[M] = float[M](row0, row1, row2);
        */

        float row[M] = float[M](1.0, r, r*r);
        
        for (int i=0; i<M; i++) {
            ATF[i] += row[i] * F[k];
        }
        for (int i=0; i<M; i++) {
            for (int j=0; j<M; j++) {
                addAA_ATA(i, j, row[i] * row[j]);
            }
        }
    }

    /*
    // --- Copy ATA → Awork, ATF → bwork ---
    for (int i=0; i<M; i++) {
        bwork[i] = ATF[i];
    }
    for (int i=0; i<M; i++) {
        for (int j=0; j<M; j++) {
            setAA_Awork(i, j, getAA_ATA(i,j));
            // if (i==0 && j==0)  addAA_Awork(i, j, 0.001);
        }
    }
    */

    // --- Gauss–Jordan elimination ---
    Gauss_Jordan_elimination(/*ATA_AA, ATF,*/ c); // ; Awork) temp.s: bwork

    /*
    // --- Copy solution to c[M] = model[M] ---
    for (int i=0; i<M; i++) {
        c[i] = bwork[i];
    }
    */
    /*
    // if M < M_max, (alt. (R,M) )
    // First M entries contain c0..c(M-1), rest = BadNaN
    for (int i=0; i<M_max; i++) {
        if (i < M) model[i] = bwork[i];
        else        model[i] = BadNaN;
    }
    */
}


float evalRadial(vec2 dx, float c[M]) {
    float r = length(dx);
    return c[0] + c[1]*r + c[2]*r*r;
}


float demodel( in vec2 dx, in float model[M]) {
   return evalRadial(dx, model);
}

void do_model(
    in vec2 dx_buffer[N], in float F[N],
    out float model[M]
) {
    fitRadialModel(dx_buffer, F, model);
}


#endif  // MODELFAMILY == RBF_BASES_MODELFAMILY

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

float SAMPLER(vec2 p, float anim_time) {
    // anim_time was t2
    #if ANIMATION_SOURCE == _BY_MOUSE_ANIMATIONSOURCE
    	float a1 = iMouse.x / iResolution.x * PI*2.;
    	float a2 = iMouse.x / iResolution.x * -PI*2. + PI*.5;
    #endif
    #if ANIMATION_SOURCE == _BY_TIME_ANIMATIONSOURCE
        float T = 3.;  // period
        float a1 = smoothstep(0., T*.75, mod(anim_time, T)) * PI*2.;  // stop for T*.25s
        float a2 = smoothstep(0., T*.75, mod(anim_time, T)) * -PI*2. + PI*.5;  // stop for T*.25s, +PI/2 phase
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

vec3 draw_trace(float d, vec2 p, vec2 ro, vec2 rd, float anim_time4) {
  vec3 col = vec3(0);
  vec3 line = vec3(1, 1, 1);
  vec2 _ro = ro;

  for (int i = 0; i < TRACE_STEPS; i++) {
    float t = SAMPLER(ro, anim_time4);
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


void fill_dx_buffer_using_grid(in float delta, inout vec2 dx_buffer[N], int NN) {

    bool something_wrong = false;
    //resetting
    for(int ctr=0; ctr<N; ctr++) {
       dx_buffer[ctr].xy = vec2(0.0);
       // dx_buffer[ctr].xy = vec2(0.1);
    }

    int GM = 3;
    {
    int ctr=0;
    for(int i=0;i<GM;i++) {
       for(int j=0;j<GM;j++) {
       // if (ctr<0)
       // if (!(i-1 ==0 && j-1 == 0))
       {
          dx_buffer[ctr].x = float(i-1) * delta;
          dx_buffer[ctr].y = float(j-1) * delta;

       }
       ctr += 1;
       }
    }
    if(!(N==ctr)) {
        something_wrong = true;
    }
    }
    // something_wrong = true;

     
     
    // ask compiler to guarantee ctr<GM*GM
    
    //runtime assert:
    // N == GM*GM;
    if(!(N==GM*GM) ) {
        something_wrong = true;
    }
    // something_wrong = true;
    // if(false) {
    // if (true) {
    if (something_wrong) {
        dx_buffer[0].xy = vec2(BadNaN, BadNaN);
        dx_buffer[1].xy = vec2(BadNaN, BadNaN);
    }
}
void fill_dx_buffer_using_circle(in float delta, inout vec2 dx_buffer[N], int NN) {
    for(int i = 0; i < N; i++) {
        float theta = float(i)/float(N) * 360.0 * PI / 180.0;
        // vec2 uv = vec2(float(i-1),float(j-1));
        vec2 uv = vec2(cos(theta),sin(theta));
        dx_buffer[i].x = uv.x * delta;
        dx_buffer[i].y = uv.y * delta; 
    }
}

float evaluate_smoothness(float delta, vec2 x0, float anim_time) {
    /*
         delta: radius of neighbourhood samples. default: 0.01
               (parameter for neighbourhood radius)
    */

    float coeffs[M];
    float SDF_buffer[N];
    // vec2 dx === uv === x0;
    vec2 dx_buffer[N]; // = x_buffer[k] - x0;
    fill_dx_buffer_using_grid(delta, dx_buffer, N);
    // ask compiler to guarantee ctr<GM*GM
    fill_dx_buffer_using_circle(delta, dx_buffer, N);
    
    //runtime assert:
    // N == GM*GM;

    for(int i=0;i<N;i++) {
       SDF_buffer[i] = 0.0;
    }
    for(int i=0;i<N;i++) {
       // if (i != 4)
       SDF_buffer[i] = SAMPLER(dx_buffer[i] + x0, anim_time);
       //if (i == 4)
       {
       float rnd = RAND(x0*1.0, float(i)+2.0 * iTime * iResolution.y/1000.0);
       SDF_buffer[i] += rnd * 0.000013 * 10.0 * 0.0; 
       }
    }
    
    // stores into: coeffs
    // computeTaylorCoefficients( N, dx_buffer, SDF_buffer, coeffs );
    do_model(dx_buffer, SDF_buffer, coeffs );
    vec2 dx0 = vec2(0.); // taylor around 0 (in dx_buffer coords)
    // float sdf_approx = evaluateTaylor(dx0, coeffs);
    float sdf_approx = demodel(dx0, coeffs);
    return sdf_approx;
}



vec2 project_closest_point_basedon_taylor(vec2 x0, float anim_time) {

    // 1. Gather local patch for Taylor model
    float F[N];
    vec2 DX[N];
    float r = 0.01;
    int idx = 0;

    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            vec2 d = vec2(float(i-1)*r, float(j-1)*r);
            DX[idx] = d;
            F[idx]  = SAMPLER(x0 + d, anim_time);
            idx++;
        }
    }

    // Fit Taylor series
    float a[M];
    // computeTaylorCoefficients(N, DX, F, a);
    do_model(DX, F, a);

    // Newton solve T(delta)=0
    vec2 delta = vec2(0.0);

    for(int it=0; it<3; it++){
        float f0 = demodel(delta, a);

        float eps = 1e-4;
        float fx = (demodel(delta + vec2(eps,0), a)
                  - demodel(delta - vec2(eps,0), a)) / (2.*eps);
        float fy = (demodel(delta + vec2(0,eps), a)
                  - demodel(delta - vec2(0,eps), a)) / (2.*eps);
        vec2 g = vec2(fx,fy);

        float denom = dot(g,g) + 1e-8;
        delta -= (f0/denom) * g;
    }

    // 4. Final closest point
    return x0 + delta;
}


// Based on SDF
vec2 project_closest_point(vec2 x, float t1) {
   // return x;
   return project_closest_point_basedon_taylor(x, t1);
}

vec3 visualise_discrepancy(vec3 col, float err) {
  // err = abs(err);

  // col = vec3(0.) + abs(err * 100000.00);
  // col = vec3(0.) + abs(err * 5.00);
  // col = vec3(0.) + abs(err * 500.00);
  // col = vec3(1.) - abs(err * 500.00);
  // float non_smoothness = abs(err * 500.00);
  // col = col * (vec3(1.) - non_smoothness);
  float non_smoothness = (err * 50000.00);
  // col = col * (vec3(1.) + non_smoothness*vec3(1.,0.,0.));
  // vec3 red_mark = vec3(1.,0.,0.) * non_smoothness;
  vec3 red_mark_ = clamp(vec3(1.,0.,0.) * (+non_smoothness), 0.,1.);
  vec3 blue_mark = clamp(vec3(.0,0.,1.) * (-non_smoothness), 0.,1.);
  // treat as alpha, using mix()
  float alpha_ = clamp(abs(non_smoothness), 0., 1.);
  // col = mix(col, red_mark, alpha_);
  // more pale: ice-cold
  // col = mix(1.0 - 0.4*col, red_mark, alpha_);
  col = mix(0.7 + 0.3*col, red_mark_ + blue_mark, alpha_);
  return col;
}

vec3 annotate_and_virualise(vec3 col, float d, vec2 uv,  vec2 ro, vec2 ref0, float anim_time, float t3) {

  vec2 rd = normalize(ref0-ro);

  #if ANIMATION_SOURCE == _BY_TIME_ANIMATIONSOURCE
  float mouse_button = (iMouse.z > 0.0 ? 1.0 : 0.0);
  #endif

  #if DISPLAY_STYLE == STYLE_SDF
    col = vec3(draw_distance(d, uv.xy));
    #if ANIMATION_SOURCE == _BY_TIME_ANIMATIONSOURCE
    	col -= mouse_button * vec3(draw_trace(d, uv.xy, ro, rd, anim_time));
    #endif
  #endif
  #if DISPLAY_STYLE == STYLE_SOLID_EDGES
    col = vec3(0) + 1.0 - vec3(draw_outline(d));
    #if ANIMATION_SOURCE == _BY_TIME_ANIMATIONSOURCE
    	col += mouse_button * vec3(1, 0.25, 0) * vec3(draw_trace(d, uv.xy, ro, rd, anim_time));
    #endif
    col = 1. - col;
  #endif
 
  #if DISPLAY_STYLE == STYLE_SOLID
    col = vec3(draw_solid(d));
  #endif
  #if DISPLAY_STYLE == STYLE_POLARITY
    col = vec3(draw_polarity(d, uv.xy, t3));
  #endif
  
  // Added a dot at mouse position (where the start point for annoations it)
  float indot = indotness(uv, ro, 0.015*3.0);
  col *= 1.0-(1.0 - indot);
  col *= indotness(uv, ref0, 0.015*3.0);


  return col;
}

vec2 pre_zoom(vec2 uv, vec2 screen) {
  return uv;
  // return uv * ( 1.0 + 0.04*RAND(uv*0.000000001, 2.0 * iTime * iResolution.y/1000.0));
  // uv = (uv+vec2(120.0,120.0))*0.7;
  // uv = 0.3*(uv/screen+vec2(0.3,0.3))*screen;
  float Z = 7.0;
  vec2 c = vec2(-0.8,-0.4);
  uv = (uv/screen-c)*screen / Z;
  return uv;
}

// see SAMPLER() as master_sdf
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
  fragCoord=pre_zoom(fragCoord, iResolution.xy);
  // anim_time (formerly, "t1") is anomation time ( SDF is tiem-dependent)
  // You can control animation speed here:
  float anim_time = iTime * 0.01;
  // float t3 = anim_time; // used in polarity // SoC!
  // float t3 = anim_time * 100.0; // used in polarity // SoC!
  float t3 = iTime; // used in polarity // SoC!. Now disentangled!
  // float t = iTime * 0.5; // turned out to be not used.
  vec2 uv = squareFrame(iResolution.xy, fragCoord);
  // float d;
  
  vec2 ro = vec2(iMouse.xy / iResolution.xy) * 2.0 - 1.0;
  ro.x *= squareFrame(iResolution.xy, iResolution.xy).x;
  // ro = mouse in coords/frame of "squareFrame"

  vec2 ref0 = vec2(0.2,0.3);


  float d = SAMPLER(uv, anim_time);
  // color to be plotted at pixel
  vec3 col = vec3(0.); // Detected a sin (uninitialise variable!). fixed.
  col = annotate_and_virualise( col, d, uv, ro , ref0, anim_time, t3);
  
  float err =  evaluate_smoothness(0.03, uv.xy, anim_time) - SAMPLER(uv.xy, anim_time);
  
  col = visualise_discrepancy(col, err);


  fragColor.rgb = col;
  fragColor.a   = 1.0;
}
