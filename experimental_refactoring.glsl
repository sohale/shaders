// experimental refactoring of an existing GLSL y someone else
// ignore tlots of notes and comments.
// TI kept the for record. This is experimental.
// code glitch-cruft!

// (private) https://www.shadertoy.com/view/3cdfz4
// refactored and tinkereed-with by Sohail

// Some refactoring experiments on fork of an existing shader: --  forked from https://www.shadertoy.com/view/ddXczM --  "golfing 567 chars  2015 shader https://shadertoy.com/view/4lXSWl "


float Pr_deeper(vec2 Rect) {
    // propability of going deeper steps: (structure of tree)
    //     depth n -> Pr = Pr_deeper^n.
    // 0.7 => stops :pr=0.3
    // float Pr_deeper = 0.7;
    // float Pr_deeper = 0.7;
    // return Pr_deeper;
    // return 0.7;
    // return 0.3; // reamins exponentially shallow!
    // return 0.7;
    // return 0.7;
    // return 0.9;
    // float thinness = min(Rect.x,Rect.y);
    // return thinness > 100. ? 0.7 : 0.2;
    return 0.7;
}

vec4 color_depth_dependent(int depth, vec3 deep_color, vec3 shallow_color) {
    /*
      `depth`: how much depth this pixel has achieved
    */
    /*
    // the more, the deeper: 1=deep black
    // (int depth_achieved)
    // float DEPTH_SCALE = 10.0; // constant
    float DEPTH_SCALE = 5.0; // constant
    // float border_alpha = 1. - exp(-1./10.0*float(depth_achieved));
    float border_alpha = 1.- exp(-1./DEPTH_SCALE*float(depth_achieved));
    border_alpha = clamp(border_alpha,0.,1.);
    // vec3 color_border_blackish = mix(vec3(0.3), vec3(0.), border_alpha);
    vec3 color_border_blackish = mix(shallow_color, deep_color, border_alpha);
    vec4 color_border = vec4(color_border_blackish.xyz,1.0);
    return color_border;
    */
    float DEPTH_SCALE = 5.0; // constant
    float alpha = 1.- exp(-1./DEPTH_SCALE*float(depth));
    alpha = clamp(alpha,0.,1.);
    vec3 c = mix(shallow_color, deep_color, alpha);
    return vec4(c, 1.0);
}

vec4 legend(vec2 uv) {
    // uv = uv / iResolution.xy;
    ivec2 numitems = ivec2(16, 2);
    ivec2 uvi = ivec2(
      int(floor(uv.x*float(numitems.x))),
      int(floor(uv.y*float(numitems.y)))
    );
   
    vec4 color_border = color_depth_dependent(uvi.x, vec3(0.), vec3(0.3));
    vec4 color_bg = color_depth_dependent(uvi.x, vec3(1.), vec3(0.7));    
    
    vec4 col = (uvi.y == 0)? color_border : color_bg;
    // return vec4(col, 1.0);
    return col;
}
bool do_legend( vec2 pixXY , out vec4 O) {
    // vec2 legend_box = vec2(0.2, 0.02);
    // vec2 legend_box = vec2(0.2, 0.02) * iResolution.xy * 2.0;
    // from top-left alignment:
    pixXY.y = iResolution.y-pixXY.y;
    vec2 legend_box = vec2(0.2, 0.07) * iResolution.xy * 2.0;
    if(pixXY.x < legend_box.x && pixXY.y < legend_box.y) {
       O = legend(pixXY/legend_box);
       bool exit = true;
       return exit;
    }
    bool exit = false;
    return exit;
}

void decide_panzoom(out vec2 U_, out vec2 Rect_, in vec2 pixU) {

    // controls animation (pan/zoom)
    // float ttimee1=iTime*0.02+2.0;
    // float FASTNESS = 0.02;
    float FASTNESS = 0.2;
    float ttimee1=iTime*FASTNESS+2.0;

    // deprecated:
    // controls LOD
    float ttimee2=iTime*0.01+7.0;
    // float mlod2 = .2+.2*sin(ttimee2);
    // size? smaller = more dense
    // float mlod2 = 0.03;
    // float mlod2 = 0.3;
    // float mlod_depth = 0.3;
    // renamed: it's a probablity!
    // moved down

    // iterative:
    // starts by whole screen
    vec2 Rect = iResolution.xy;
    // vec2 R = ^
    // vec2 anim_pos = R * (.5+.5*sin(.1*ttimee1+vec2(0,33)));
    //vec2 anim_pos;
    //anim_pos.x = R.x * (.5+.5*sin(.1*ttimee1+0.0));
    //anim_pos.y = R.y * (.5+.5*sin(.1*ttimee1+33.0));
    // vec2 m = iMouse.z>0. ?  iMouse.xy : anim_pos;
    // ^ mouse control or demo mode
    // vec2 anim_pos_dR = 
    vec2 animpos_t = (.5+.5*sin(.1*ttimee1+vec2(0,33)));
    // originally the scale was made to match its mouse-clamp/source!
    // vec2 animpos_mouse = iMouse.xy;
    vec2 animpos_mouse = iMouse.xy * 0.001;
    bool CLICKED = iMouse.z>0.;
    // vec2 animpos = CLICKED ? animpos_mouse : animpos_t;
    vec2 animpos = animpos_t;

    // fovea = centre of retina!
    vec2 screencentr = Rect/2.; 
    // vec2 m = R * anim_pos_dR;
    float zoom1 = animpos.y;
    float zoom = ( 1. - zoom1 ) * 4.;
    // x positiom: compared to center of screen
    // float x_target_fovea = 
    vec2 pan_for_target = vec2(
       8.* (Rect.x * animpos.x - screencentr.x ),
       0.0
    );
    // translate & zoom
    vec2 U = (pixU - pan_for_target) / zoom;
    // ^ coord are in Rect? (non-nrmalised?)
    // R shall die here: in that cordinate, R has no place:
    // U/R
    // vec2 Ur = U / Rect;
    // but it changes (updates) Rect!
    
    U_ = U;
    Rect_ = Rect;
    
}

void mainImageAA( out vec4 O, out int outmode, vec2 pixU_ )
{
    /*
    vec2 legend_box = vec2(0.2,0.02);
    if(pixU.x < legend_box.x && pixU.y < legend_box.y) {
       O = legend(pixU/legend_box);
    }
    */
    if (do_legend(pixU_, O)) {
        return;
    }

    vec2 U; // useful
    vec2 Rect; // by-product
    // pixU_ : consumed
    decide_panzoom(U,Rect, pixU_);
    // a randomness but fixed by (Rect,I)
    #define H(v) fract(4e4*sin(dot( v+Rect*I, vec2(13.46,41.74))+17.34)) 
    // hash(seed),I=cellId
   
    // float mlod_depth = 0.3;
    // float Pr_stop_depth = 0.3;
    // propability of going deeper steps: (structure of tree)
    //     depth n -> Pr = Pr_deeper^n.
    // float Pr_deeper = 0.7;
    
    // todo: this probability needs to depend on size, etc
    vec2 I;
    int depth_achieved = -1;
    for (int i; i++ < 16; ) {     // browse LODs
        depth_achieved = i;
        I = ceil(U/Rect)-.5;  // cell id
        // random depth (no)
        // randomly, stop
        if ( H() < (1. - Pr_deeper(Rect)) ) break;     
        int diraixs = H(.1) < .5 ? 0 : 1 ;
        // ^ random dir of subdiv
        // Rect[ diraixs] /= 3.;
        // Rect[ diraixs] /= 2.;
        // Rect[ diraixs] /= 5.; // too clump-y
        // Rect[ diraixs] /= 1.02; // almost regular!
        // Rect[ diraixs ] /= 2.5;
        // Rect[ diraixs ] /= diraixs==0?2.5:2.0;
        Rect[ diraixs ] /= 2.5;
    }
    // NOT! : ratio = Shrinkage aspect ratio (for last round? leaf?)
    vec2 rel_coords = (vec2(.5) - abs( U/Rect - I ));
    // ^ coords inside the leaf rectangle
    // true coords?: (but shifted by border?)
    vec2 U2 = Rect * rel_coords;    // coord relative to border

    // #define U EROREAORW
    // vec4 color_black = vec4(0); // border
    // the more, the deeper: 1=deep black
    vec4 color_border = color_depth_dependent(depth_achieved, vec3(0.), vec3(0.3));

    // vec4 color_white = vec4(1); // bg
    // vec4 color_bg = vec4(1); // white
    vec4 color_bg = color_depth_dependent(depth_achieved, vec3(1.), vec3(0.7));
    
    // conditional usage:
    // colored leaf
    vec4 color_random = cos(6.3*H(1.) +vec4(0,21,23,0));

    float thinness = min(U2.x,U2.y);
    float WALL_THICKNESS = 1.5/2.;
    // probablility of white? (non-color_random)
    float Pr_WHITE = .8;
    // float Pr_WHITE = .3;

    /*
    O = (thinness < WALL_THICKNESS) ?
            color_border  // border
        : H(.2) < Pr_WHITE ?
            // color_white    // bg
            color_bg
        :
            color_random;  // colored leaf
    */
   if (thinness < WALL_THICKNESS) {
        O = color_border;  // border
        outmode = 1;
   } else if ( H(.2) < Pr_WHITE ) {
        // color_white    // bg
        O = color_bg;
        outmode = 2;
    } else {
        O = color_random;  // colored leaf
        outmode = 3;
    }

            
}

// crude AA
void mainImage( out vec4 O, vec2 pixU )
{

    // remains the same across AA:
    vec2 U,Rect;
    decide_panzoom( U,Rect, pixU);
    float AA_size = min(Rect.x,Rect.y);

    vec4 OAccum=vec4(0.);
    const int Mx = 5;
    const int My = 5;
    const int AAN = Mx*My;
    vec2[AAN] shifts; // = vec2[](    vec2(0.,0.)    );
    for(int aai=0;aai<AAN;aai++) {
       float aau = float(aai)/float(AAN);
       
       // float AA_scale = 1.0;
       // float AA_scale = 1000. / AA_size;
       //
       // float AA_scale = AA_size / 10.0; // actually remains fixed?
       // actually remains fixed:
       // visible "small" circle
       // float AA_scale = AA_size / 100.0;
       float AA_scale = AA_size / 100.0 * 0.2;
       
       
       // AA-convolution modes:
       
       // rect-grid method
       vec2 aauv = vec2( float(aai%Mx) / float(Mx),float(aai/Mx)/float(My));
       // ideas:
       //    1. rotate the grid mildly?
       //    2. condende in the middle: to sample more (for thinner lines?)
       // aauv = ROTATION *...* aauv
       shifts[aai] = aauv;
       shifts[aai] += 0.01*vec2(sin(float(aai)), cos(float(aai)+iTime));
       shifts[aai] *= AA_scale * 3.0;
       // grid is alias! we need slants
       
       
       /*
       // circular method
       const float PI = 3.14159265;
       float theta = (2.0*PI)*float(aai)/float(AAN);
       float r = 1.0; // more robust
       // float r = aau; // flickers a bit
       // spiral method!
       // float r = aau*aau;
       shifts[aai]= vec2(cos(theta), sin(theta)) * r;
       // todo: this should depend on zoom.
       //  ^ done!
       shifts[aai] *= AA_scale;
       */
     }
    for(int aai=0;aai<AAN;aai++) {
       int outmode; // actually not useful
       vec4 Otemp=vec4(0.);
       mainImageAA(Otemp, outmode, pixU + shifts[aai]);
       OAccum += Otemp;
    }
    O = OAccum / float(AAN);
}
