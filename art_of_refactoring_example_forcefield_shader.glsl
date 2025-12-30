// Demonstrating my "method of semantic refactoring"
// Starting with pre-existing shader by 01000001 https://www.shadertoy.com/view/4ffXzf


#define saturation .1
// .7 is nice


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

// uv → * scale = ChW
// pix → ( . - scr/2)/scr.y = uv

vec2 pixel_to_uv(vec2 pix_xy) {

    vec2 screen = iResolution.xy;
    vec2 uv = ( pix_xy - screen / 2.0 ) / screen.y;
    return uv;
}
vec2 uv_to_pixelxy(vec2 uv) {
    vec2 screen = iResolution.xy;
    return uv * screen.y + screen / 2.0;
}
vec2 pixeldxy_to_uv(vec2 pix_dxy) {

    vec2 screen = iResolution.xy;
    vec2 uv = ( pix_dxy ) / screen.y;
    return uv;
}
vec2 uv_to_pixeldxy(vec2 uv) {
    vec2 screen = iResolution.xy;
    return uv * screen.y;
}

vec2 uv_to_chargesworld(vec2 uv) {
    const float scale_m_ = 12.*2.0; 
    return uv * scale_m_;
}
vec2 chargesworld_to_uv(vec2 chw) {
    const float scale = 12.*2.0; 
    return chw / scale;
}

vec2 pixelxy_to_chargesworld(vec2 pix_xy) {
    // sht
    const float scale_m_ =  12.*2.0;
    // return pixel_to_uv(pix_xy) * scale_m_;
    return uv_to_chargesworld(pixel_to_uv(pix_xy));   
}
// pixeldxy_to_chw
vec2 pixeldxy_to_chargesworld(vec2 pix_dxy) {
    const float scale_m =  12.*2.0;
    // return uv_to_chargesworld(pixeldxy_to_uv(pix_dxy));
    // return pixeldxy_to_uv(pix_dxy) * scale_m;
    vec2 screen = iResolution.xy;
    return pix_dxy  / screen.y * scale_m;
}
vec2 chargesworld_to_pixeldxy(vec2 dcw) {
    const float scale_m =  12.*2.0;
    vec2 screen = iResolution.xy;
    // return pix_dxy  / screen.y * scale_m;
    return dcw / scale_m * screen.y ;
}


vec2 chargesworld_to_pixelxy(vec2 cw_coords) {
    vec2 screen = iResolution.xy;
    const float scale_m_ =  12.*2.0;
    vec2 uv = cw_coords / scale_m_;
    return uv * screen.y + screen * 0.5;
}


float atan2(in float y, in float x)
{
    // `atan2` by HuaTham - https://stackoverflow.com/a/26070411
    // Retrieved 2025-12-30, License - CC BY-SA 3.0
    bool s = (abs(x) > abs(y));
    return mix(pi/2.0 - atan(x,y), atan(y,x), s);
}
mat2 rotMat2(float rad)
{
    float s = sin(rad);
    float c = cos(rad);
    return mat2( c, -s, s, c );
}


// The MIT License. Copyright © 2017 Inigo Quilez
// Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions: The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software. THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

const uint k = 1103515245U;  // GLIB C
//const uint k = 134775813U;   // Delphi and Turbo Pascal
//const uint k = 20170906U;    // Today's date (use three days ago's dateif you want a prime)
//const uint k = 1664525U;     // Numerical Recipes

vec3 hash33( uvec3 x )
{
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;
    x = ((x>>8U)^x.yzx)*k;
    
    return vec3(x)*(1.0/float(0xffffffffU));
}
vec2 hash33_2d(ivec2 seed2d) {
   uvec3 seed3d = uvec3(seed2d, 79.0);
   return hash33(seed3d).xy;
}


// coords:
// pixel-xy , also: mousexy
// uv (scaled) -- current
// floor-able (currently the same as uv)
// another for the ROD ( 1.0-ness in rod)

//  in the input-uv
// in "charges-world"
// in rod world

// now that rod-world uses pixel-xy,
// we need to collapse one?: uv, or, charges.
// No: I just need to get the grid one in one consistent coord ( RADIUS0 ad RADIUS1 are so )


// vec2 to_griddable() { }

void locate_in_node1_deprecated(in vec2 _uv, out vec2 true_node_centre, out vec2 cuv3) {
    /*
    float alph = 1.0;
    ivec2 cellint = ivec2(floor(uv * alph));
    vec2 cell_node_center = vec2(cellint)/alph;
    // vec2 cuv_ = (uv-cell_node_center) * 2. - 1.;
    vec2 cuv_ = uv-0.5 - cell_node_center;
    vec2 true_node_centre = cell_node_center + 0.5;
    vec2 cuv3 = cuv_ * mat2(0,1,-1,0) * 2.0;
    */
    /*
    // vec2 grid_phase = 0.0*fract(vec2(0.0 + fract(iTime * 18.5), 0.0)) * 1.0;
    vec2 grid_phase = 1.0*vec2(.2, 0) * rotMat2(iTime * 2.0*pi / 2.1);
    float alph = 1.0;
    ivec2 cellint_id = ivec2(floor((uv + grid_phase) * alph));
    true_node_centre = vec2(cellint_id)/alph + 0.5 + grid_phase;
    vec2 cuv_ = uv  - true_node_centre;
    vec2 cuv3 = cuv_ * mat2(0,1,-1,0) * 2.0;
    // cellint = node identity
    //      = nodal domain id
    cuv3_ = cuv3;
    */
    
    // wow: rot90 was part of the line segment
    
    // wow: scale_ belongs INSIDE this. Because of SDF! (how?)
    const float scale_ = 12.*2.0;
    // uv = uv * scale_;

    vec2 grid_phase = 1.0*vec2(.2, 0) * rotMat2(iTime * 2.0*pi / 2.1);
    float alph = 1.0;
    
    vec2 griddable_xy = (_uv * scale_ + grid_phase) * alph;
    // keep this `ivec`: cellint_id = node identity  = nodal domain id
    ivec2 cellint_id = ivec2(floor(griddable_xy));
    vec2 cellint_id_vec2 = vec2(cellint_id);
    // true_node_centre = vec2(cellint_id)/alph + 0.5 + grid_phase;
    // true_node_centre = vec2(cellint_id)/alph + 0.5 + grid_phase;
    // anti_griddable: separated this to antigraiddable
    vec2 antigraiddable = (cellint_id_vec2/alph - grid_phase)/scale_;
    //now I split <>
    /*
    vec2 antiantigraiddable = (antigraiddable*scale_ + grid_phase)*alph;
    true_node_centre = antiantigraiddable/alph + 0.5 + grid_phase;
    */
    /*
    // and combine!
    true_node_centre = antigraiddable*scale_ + 0.5 + 2.00 * grid_phase;
    

    // vec2 cuv_ = (_uv  - true_node_centre/scale_);
    // vec2 cuv_ = (_uv * scale_  - true_node_centre)/scale_;
    vec2 cuv_ = (_uv  - true_node_centre/scale_);
    // true_node_centre = true_node_centre / scale_;
    cuv3 = cuv_ * 2.0;  // keeps it withing the radius
    // cellint_id = node identity  = nodal domain id
    */
    // and combine more
    
    true_node_centre = antigraiddable*scale_ + 0.5 + 2.00 * grid_phase;
    cuv3 = (_uv  - true_node_centre/scale_)*2.0;
}


// simplifying again
// a snap-to-grid
// void locate_in_node(in vec2 _uv, out vec2 true_node_centre_uv, out vec2 nodal_uv0, out vec2 RADIUS01_uv) {
void locate_in_node(in vec2 _uv, out vec2 true_node_centre_uv, out vec2 RADIUS01_uv) {

    const float scale_ = 12. * 2.0;
    // const float scale_ = 5.0;
    // alpha & grid_phase are about discretisation: floor.  It should have been: floor(x*alpha+phase)
    float alph = 1.0;
    // try increasing ^

    
    vec2 grid_phas_e_ = -1.0*vec2(.2, 0) * rotMat2(iTime * 2.0*pi / 2.1);
    // vec2 grid_phase2 = grid_phase - 0.5;
    // vec2 grid_phase3 = (grid_phase_ - 0.5)*alph;
    // vec2 grid_phase0_ = -grid_phase_ * alph;
    vec2 grid_phase0_ = grid_phas_e_ * alph;
    // vec2 grid_phase30 = (grid_phase0 - grid_phase3)/alph;
    // vec2 grid_phase30 = (grid_phase0 - (grid_phase_ - 0.5)*alph)/alph;
    // vec2 grid_phase30 = (grid_phase0 - grid_phase_*alph + 0.5*alph )/alph;
    // vec2 grid_phase30 = (grid_phase0 - grid_phase0 + 0.5*alph )/alph;
    // vec2 grid_phase30 = vec2(0,0)+(0.5*alph )/alph;
    // vec2 grid_phase30 = vec2(0,0) + 0.5;

    // vec2 griddable_xy = (_uv * scale_ + grid_phase) * alph;
    vec2 griddable_ixy = (_uv * scale_) * alph - grid_phase0_;

    // keep this `ivec`: cellint_id = node identity  = nodal domain id
    ivec2 cellint_id = ivec2(floor(griddable_ixy));
    vec2 cellint_id_vec2 = vec2(cellint_id);
    // ^ anchor/corner point ^
    // aha, anchor is different to the center.
    // besides phase needds alpha too.


    // back to the input-uv
    // vec2 antigraiddable_uv = (cellint_id_vec2/alph - grid_phase)/scale_;
    // vec2 h2_uv = chargesworld_to_uv(vec2(0.5, 0.5)); //  =  0.5 / scale;

    // true_node_centre_uv = antigraiddable_uv + h2_uv;
    //  = true_node_centre_uv_offs0
    
    // true_node_centre_uv = (cellint_id_vec2/alph - grid_phase)/scale_ + 0.5 / scale_;
    // true_node_centre_uv = ( cellint_id_vec2/alph - grid_phase + 0.5 )  / scale_;
    // true_node_centre_uv = ( cellint_id_vec2/alph - ( grid_phase - 0.5) )  / scale_;
    // true_node_centre_uv = ( cellint_id_vec2/alph - grid_phase2 )  / scale_;
    // true_node_centre_uv = ( cellint_id_vec2/alph - grid_phase2*alph/alph )  / scale_;
    // true_node_centre_uv = ( (cellint_id_vec2 - grid_phase2*alph)/alph )  / scale_;
    // true_node_centre_uv = ( (cellint_id_vec2 - grid_phase3)/alph )  / scale_;
    // true_node_centre_uv = ( (cellint_id_vec2 - grid_phase0 + grid_phase0- grid_phase3)/alph )  / scale_;
    // true_node_centre_uv = ( (cellint_id_vec2 - grid_phase0)/alph + (grid_phase0 - grid_phase3)/alph )  / scale_;
    // true_node_centre_uv = ( (cellint_id_vec2 - grid_phase0)/alph + grid_phase30 )  / scale_;
    // true_node_centre_uv =  (cellint_id_vec2 - grid_phase0)/alph  / scale_ + grid_phase30   / scale_;

    // vec2 ancho = (cellint_id_vec2 + grid_phase0_)/alph  / scale_;
    // true_node_centre_uv =  ancho + grid_phase30 / scale_;
    // true_node_centre_uv =  ancho + 0.5 / scale_;

    // true_node_centre_uv =  (cellint_id_vec2 + grid_phase0_ + 0.5 / scale_  * alph * scale_) /alph / scale_;
    // true_node_centre_uv =  (cellint_id_vec2 + grid_phase0_ + 0.5 * alph ) /alph / scale_;
    // true_node_centre_uv =  (cellint_id_vec2 + 0.0*0.5 * alph + grid_phase0_ ) /alph / scale_;
    // ^ bug detected (for alph):  this should not depend on alpha: It's in "floor"-space ( ivec )
    // The center (of attention:
    // Why it should not be at zero? because floor has a bias. cuts at [0,+1)
    // the anchor-center
    // in florr-able coord system:
    // vec2 node_center_ = cellint_id_vec2 + 0.5 * alph;
    // by shifting 0.5, now radius is available
    // BUG FIX:
    vec2 node_center_ = cellint_id_vec2 + 0.5;
    // this is not arbitrary. it needs tobe at center of the region, which is cuts at [0,+1).
    // the region is: [node_center_, node_center_+1)

 
    vec2 noise2 = hash33_2d(cellint_id); // (cellint_id_vec2);
    // node_center_.xy += noise2 * 0.4;
    node_center_.xy += (noise2  - 0.5) * 2.0 * 0.8;

    true_node_centre_uv =  (node_center_ + grid_phase0_ ) /alph / scale_;


    // vec2 h3_uv = h2_uv  - 0.5/24.0;
    // nodal_uv0 = true_node_centre_uv_offs0 - h3_uv;
    // nodal_uv0 = true_node_centre_uv - h3_uv;
    // nodal_uv0 = true_node_centre_uv - (h2_uv  - 0.5/24.0);
    // nodal_uv0 = (antigraiddable_uv + h2_uv) - (h2_uv  - 0.5/24.0);
    // nodal_uv0 = antigraiddable_uv + h2_uv - h2_uv  + 0.5/24.0;
    // nodal_uv0 = antigraiddable_uv + 0.5/24.0;
    // chargesworld_to_uv === div / 24.0
    // nodal_uv0 = antigraiddable_uv + chargesworld_to_uv(vec2(0.0)+0.5);
    // sht! :
    // nodal_uv0 = antigraiddable_uv + h2_uv;
    // sht!
    // nodal_uv0 = true_node_centre_uv;
    // ^ no needed anymore!

    // this turns out to be in the same coords!
    const float RADIUS1_ = 1.0 / 2.0;
    const float RADIUS0_ = 0.9 / 2.0;
    // cuv3_2 = cuv3_2 * scale_m_;
    // float RADIUS0_xy = chargesworld_to_pixeldxy(vec2(RADIUS0_, 0.0)).x;
    // float RADIUS1_xy = chargesworld_to_pixeldxy(vec2(RADIUS1_, 0.0)).x;
    vec2 RADIUS01_xy = chargesworld_to_pixeldxy(vec2(RADIUS0_, RADIUS1_));
    // what? RADIUSi were in ChW?
    // ok, now convert to uv, and back to pixel/xy:
    RADIUS01_uv = pixeldxy_to_uv(RADIUS01_xy);
    
}

const vec2 ZERO = vec2(0,0);
const vec2 OME = vec2(1,1);

//float rod_shape(vec2 local_vec_xy, vec2 force_dxy, float RADIUS0_xy, float RADIUS1_xy) {
// float rod_shape_xy(vec2 local_vec_xy, vec2 force_dxy, float RADIUS0_xy, float RADIUS1_xy) {
float rod_shape_xy(vec2 local_vec_xy, vec2 force_dxy, vec2 RADIUS01_xy, float thickness_xy) {
    // local_coords -> local_coords_cw -> local_vec_cw
    // local, relative, etc.
    // local_vec_cw -> local_vec_xy = dxy
    
    // vec2 duv = cw_to_uv(local_vec_cw)
    // vec2 dcw = uv_to_cw(duv)
    // vec2 dxy = chargesworld_to_pixelxy(local_vec_cw);
     vec2 dxy = local_vec_xy;
    // vec2 dcw = pixeldxy_to_chargesworld(dxy);
    vec2 dcw = pixeldxy_to_chargesworld(dxy);
    // vec2 dcw =  local_vec_xy  / screen.y * scale_m;
    
    vec2 radius_errorbar = pixeldxy_to_chargesworld(RADIUS01_xy);
    float RADIUS0_chw = radius_errorbar.x;
    float RADIUS1_chw = radius_errorbar.y;
    // vec2 force_ = pixeldxy_to_chargesworld(force_dxy);
    vec2 force_ch = pixeldxy_to_chargesworld(force_dxy);

    // const float THICKNESS = 0.3;
    //corrected!
    // in ChW coord syste?!
    // const float THICKNESS =  0.3/2.0;
    // const float THICKNESS =  0.3/2.0;
    // float THICKNESS0 = chargesworld_to_pixeldxy(ZERO+0.3/2.0).x;
    float THICKNESS_ = pixeldxy_to_chargesworld(ZERO+thickness_xy).x;
    const mat2 rot90 = mat2(0,1,-1,0);
    return 1.0
         //  // Create bar -> Make it look nice : side thickness 0.3
        * smoothstep(THICKNESS_, .0, 
            abs(dot(normalize(force_ch), rot90 * dcw )) 
         )

        // limit length of bars to 1 cell width radius
        // that (cell width) is called "diameter" !
        // diameter = 1 = RADIUS1_*2.0, RADIUS1_ ~ 0.5
        // limit length of bars to 1 cell width radius
        // * smoothstep(RADIUS1_chw*2.0, RADIUS0_chw*2.0, length(dcw)*2.0)
        * smoothstep(RADIUS1_chw, RADIUS0_chw, length(dcw))
    ;
}

vec3 color_from_force(vec2 force) {
    // Colour based on angle
    // Angle of bar
    // float t = atan(force.x/force.y) + (sign(force.y)>0.?-pi:0.);
    // float t = 1.0;
    float t = atan2(force.y, force.x);    
    vec4 colour = mix(vec4(1), (cos(t + vec4(0., 2.*pi/3., 4.*pi/3., 0.))*.5+.5), vec4(saturation));
    return colour.rgb;
}


float pole_dots_color(vec2 pix_xy) {

    // const float POLES_RADIUS = 15.0*3.0;
    float POLES_RADIUS = max(15.0, 15.0/600.0*iResolution.y);

    float poles_dots =  0.0;
    for(int i = 0; i < n; i++) {
        // d2_shape = max(d2_shape, smoothstep(0.2*10.0, 0.0, 1.0*distance(uv_cw, charges[1].xy)));
        // d2_shape = max(d2_shape, smoothstep(5.0, 0.0, distance(pix_xy, cw_to_pixel(charges[1].xy))));
        float d_i = distance(pix_xy, chargesworld_to_pixelxy(charges[i].xy));
        // d2_shape = max(d2_shape, smoothstep(5.0, 0.0, d_i) );
        // debugging
        //
        // d2_shape = 2.0-d_i;
        // d2_shape =  smoothstep(5.0, 0.0, d_i);
        float pole_i_sdf_shape =  smoothstep(POLES_RADIUS, POLES_RADIUS-2.0, d_i);
        poles_dots = max(poles_dots, pole_i_sdf_shape);
    }
    poles_dots = clamp(poles_dots, 0.0, 1.0);
    return poles_dots;
}

void mainImage( out vec4 O, vec2 pix_xy )
{
    vec2 uv_ =  pixel_to_uv(pix_xy);
    
    // const float scale_m_ =  12.*2.0; // 1.0 ;
    // ^ aha, now I see this is to convert from uv-world to charges-world

    // vec2 true_node_centre_uv_;
    // vec2 cuv3_;
    // locate_in_node(uv_ , true_node_centre_uv_);
    vec2 true_node_centre_uv_offs0;
    // vec2 nodal_uv0;
    vec2 RADIUS01_uv;
    // locate_in_node(uv_ , true_node_centre_uv_offs0, nodal_uv0, RADIUS01_uv);
    locate_in_node(uv_ , true_node_centre_uv_offs0, RADIUS01_uv);
    // nodal_uv0 = true_node_centre_uv_offs0;
    // ^ absorb-away / join.

        // in "charges-world" !
    // vec2 true_node_centre = true_node_centre_uv*scale_m_ + 0.5 + 2.00 * grid_phase;
    // vec2 true_node_centre = true_node_centre_uv*scale_m_ + 0.5; // + 2.00 * grid_phase;
    // lol:
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv) + 0.5; // + 2.00 * grid_phase;
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv) + vec2(0.5,0.5); // + 2.00 * grid_phase;
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv) + chargesworld_to_uv(vec2(0.5,0.5)); // + 2.00 * grid_phase;
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv + chargesworld_to_uv(vec2(0.5,0.5))); // + 2.00 * grid_phase;
    // vec2 h2_uv = chargesworld_to_uv(vec2(0.5, 0.5));
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv + h2_uv); // + 2.00 * grid_phase;
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv_ + h2_uv); // + 2.00 * grid_phase;
    // vec2 true_node_centre = uv_to_chargesworld(true_node_centre_uv_offs); // + 2.00 * grid_phase;
    
    // in rod world!
    // vec2 cuv3_ = (uv_  - true_node_centre/scale_m_)*2.0;
    // lol!
    // vec2 cuv3_ = (uv_  - true_node_centre-0.5)*2.0;
    // vec2 cuv3_ = (uv_  - (true_node_centre-0.5))*2.0;
    // LOL, lets's start over:
    // vec2 cuv3_ = (uv_  - true_node_centre/scale_m_)*2.0;
    // vec2 cuv3_ = (uv_  - (true_node_centre)/24.0)*2.0;
    // vec2 cuv3_ = (uv_  - (true_node_centre-0.5+0.5)/24.0)*2.0;
    // vec2 cuv3_ = (uv_  - (uv_to_chargesworld(true_node_centre_uv) + 0.5)/24.0)*2.0;
    // vec2 cuv3_ = (uv_  - (uv_to_chargesworld(true_node_centre_uv)/24.0 + 0.5/24.0))*2.0;
    // yea:
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre_uv + 0.5/24.0))*2.0;
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre_uv + h2_uv - h2_uv + 0.5/24.0))*2.0;
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre - h2_uv + 0.5/24.0))*2.0;
    // vec2 true_node_centre_uv1 = true_node_centre_uv + h2_uv;
    // another key step:
    // vec2 true_node_centre_uv1_ = chargesworld_to_uv(true_node_centre);
    //vec2 h2_uv = chargesworld_to_uv(vec2(0.5, 0.5));
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre_uv1 - h2_uv + 0.5/24.0))*2.0;

    // vec2 h3_uv = chargesworld_to_uv(vec2(0.5, 0.5))  - 0.5/24.0;
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre_uv1_ - h3_uv))*2.0;
    // vec2 cuv3_delta_uv = (uv_ - (true_node_centre_uv_offs - h3_uv))*2.0;
    // vec2 cuv3_delta_uv = (uv_ - true_node_centre_uv_offs + h3_uv)*2.0;
    // vec2 cuv3_delta_uv_deviant = (uv_ - true_node_centre_uv_offs + h3_uv)*2.0;



    // vec2 cuv4 = cuv3_ * scale_m_;
    
    // move below
    // vec2 cuv4 = uv_to_chargesworld(cuv3_delta_uv);


    // vec2 RADIUS01_xy2 = uv_to_pixeldxy(RADIUS01_uv);
    


    if (iMouse.z > 0.) {
       // charges[0].xy = (iMouse.xy-screen/2.0)/screen.y * scale_m_;
       charges[0].xy = pixelxy_to_chargesworld(iMouse.xy);
    }
    // modulate the charge
    // charges[0].z = sin(iTime*pi*2.0* 2.07)*8. + 10.;
    charges[0].z = sin(iTime*pi*2.0 / 5.5 )*8. + 10.;

    charges[1].xy = vec2(sin(iTime*.5), cos(iTime*1.2))*5.5;
    charges[2].xy = vec2(sin(iTime*1.3), cos(iTime*.8))*10.;

    // bug!
    // vec2 force = field(cell_node_center);
    // fixed!
    // vec2 force = field(cell_node_center+0.5);
    // vec2 force = field(true_node_centre_chw);
 

    // vec2 h3_uv = chargesworld_to_uv(vec2(0.5, 0.5))  - 0.5/24.0;
    // vec2 cuv3_delta_uv_deviant = (uv_ - true_node_centre_uv_offs + h3_uv)*2.0;
    // vec2 cuv4_delta = uv_to_chargesworld(cuv3_delta_uv);
    // vec2 cuv4_delta = uv_to_chargesworld(cuv3_delta_uv_deviant);
    // vec2 cuv4_delta_ = uv_to_chargesworld(cuv3_delta_uv_deviant);
    // vec2 cuv4_delta_xy = chargesworld_to_pixelxy(cuv4_delta_);
    // vec2 cuv4_delta_xy = chargesworld_to_pixelxy(uv_to_chargesworld(cuv3_delta_uv_deviant));
    // vec2 cuv4_delta_xy = uv_to_pixelxy(cuv3_delta_uv_deviant);
    // vec2 h3_uv = chargesworld_to_uv(vec2(0.5, 0.5))  - 0.5/24.0;
    // vec2 cuv4_delta_xy_dev_ = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs + h3_uv)*2.0);
    // vec2 nodal_dxy = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs + h3_uv)*2.0);
    // vec2 nodal_uv0 = true_node_centre_uv_offs - h3_uv;
    // vec2 nodal_dxy = uv_to_pixeldxy((uv_ - nodal_uv0)*2.0);
    // vec2 nodal_dxy = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs0)*2.0);
    // vec2 nodal_dxy2 = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs0)*2.0);
    // linear =>
    // vec2 nodal_dxy2 = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs0)*1.0)*2.0;
    // vec2 nodal_dxy1 = uv_to_pixeldxy((uv_ - true_node_centre_uv_offs0)*1.0);
    // vec2 nodal_dxy2 = nodal_dxy1*2.0;
    vec2 nodal_dxy1 = uv_to_pixeldxy(uv_ - true_node_centre_uv_offs0);


    // rod_shape is in ChW. The idea is to change it to uv, and then, pix etc
    // float d2_shape = rod_shape( cuv4_delta, force, RADIUS0, RADIUS1);
    // float d2_shape = rod_shape( cuv4_delta_xy_dev, force, RADIUS0, RADIUS1);
    

    vec2 true_node_centre_chw = uv_to_chargesworld(true_node_centre_uv_offs0); // + 2.00 * grid_phase;
    vec2 force_chw = field(true_node_centre_chw);
    
    // rod_shape_xy is now in pixel/xy coords:
    vec2 force_dxy = chargesworld_to_pixeldxy(force_chw);

    vec2 RADIUS01_xy2 = uv_to_pixeldxy(RADIUS01_uv);
    // float d2_shape = rod_shape_xy( cuv4_delta_xy_dev_, force_dxy, RADIUS01_xy2);
    // float d2_shape = rod_shape_xy( nodal_dxy2, force_dxy, RADIUS01_xy2);
    // float d2_shape = rod_shape_xy( nodal_dxy1*2.0, force_dxy, RADIUS01_xy2);
    // ^ which there reveals a subtle but satisfying bug!
    // or, may be a technique!
    // let's call it a name:
    //no, it was correct!!
    float SHRINK_FACTOR = 1.0;
    
    float thinness_pixels = chargesworld_to_pixeldxy(ZERO+0.3/2.0).x;
    // float thinness_pixels = 2.0; // pixels!
    float d2_shape = rod_shape_xy( nodal_dxy1*SHRINK_FACTOR, force_dxy, RADIUS01_xy2, thinness_pixels);


    // d2_shape = max(d2_shape, smoothstep(0.2, 0.0, distance(uv_, charges[1].xy)));
    vec2 uv_cw = uv_to_chargesworld(uv_);
    
    float poles_dots = pole_dots_color(pix_xy);
    // d2_shape = max( poles_dots, d2_shape);


    // Colour based on angle
    vec3 colour = color_from_force(force_chw);

    // brighter at higher force
    float brightness = sqrt(length(force_chw));

    vec3 colour3 = sqrt(d2_shape) * brightness * colour.rgb;
    // colour3.r = max( poles_dots, colour3.r);
    // colour3.r = max( poles_dots, colour3.r );
    // colour3 = mix(colour3, vec3(1.0,0,0)*poles_dots, poles_dots);
    vec3 poles_colour3 = vec3(1.0, 0, 0) * poles_dots;
    // colour3 = mix(colour3, poles_colour3, poles_dots);
    // colour3 = poles_colour3;
    colour3.r = max( colour3.r, poles_colour3.r );


    O = vec4(colour3, 1.0);
    // O = sqrt(d2) * brightness * colour;
    
    // O /= O + 1.;
}
