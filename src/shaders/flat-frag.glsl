#version 300 es
precision highp float;

uniform vec3 u_Eye, u_Ref, u_Up;
uniform vec2 u_Dimensions;
uniform float u_Time;

in vec2 fs_Pos;
out vec4 out_Col;

const int MAX_RAY_STEPS = 255;
const float MIN_DIST = 0.0;
const float MAX_DIST = 100.0;
const float EPSILON = 0.01;

vec3 light1Pos = vec3(4.0, 2.0, -4.0);
vec3 light2Pos = vec3(-20.0, -15.0, -10.0);

struct Bbox
{
  vec3 min;
  vec3 max;
};

// helper functions for sdf scene

float rand(float n){return fract(sin(n) * 43758.5453123);}

float sphere (vec3 p, float r, vec3 c)
{
    return distance(p,c) - r;
}

float sphere( vec3 p, float s )
{
  return length(p) - s;
}

float torus( vec3 p, vec2 t )
{
  vec2 q = vec2(length(p.xz)-t.x,p.y);
  return length(q)-t.y;
}

float box( vec3 p, vec3 b )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0))
         + min(max(d.x,max(d.y,d.z)),0.0);
}

float roundBox( vec3 p, vec3 b, float r )
{
  vec3 d = abs(p) - b;
  return length(max(d,0.0)) - r
         + min(max(d.x,max(d.y,d.z)),0.0);
}

float plane( vec3 p, vec4 n )
{
  // n must be normalized
  return dot(p,n.xyz) + n.w;
}

float cylinder( vec3 p, vec2 h )
{
  vec2 d = abs(vec2(length(p.xz),p.y)) - h;
  return min(max(d.x,d.y),0.0) + length(max(d,0.0));
}

float roundedCylinder( vec3 p, float ra, float rb, float h )
{
    vec2 d = vec2( length(p.xz)-2.0*ra+rb, abs(p.y) - h );
    return min(max(d.x,d.y),0.0) + length(max(d,0.0)) - rb;
}

float cappedCone( vec3 p, float h, float r1, float r2 )
{
    vec2 q = vec2( length(p.xz), p.y );
    
    vec2 k1 = vec2(r2,h);
    vec2 k2 = vec2(r2-r1,2.0*h);
    vec2 ca = vec2(q.x-min(q.x,(q.y < 0.0)?r1:r2), abs(q.y)-h);
    vec2 cb = q - k1 + k2*clamp( dot(k1-q,k2)/dot(k2, k2), 0.0, 1.0 );
    float s = (cb.x < 0.0 && ca.y < 0.0) ? -1.0 : 1.0;
    return s*sqrt( min(dot(ca, ca),dot(cb,cb)) );
}

float dot2( in vec3 v ) { return dot(v,v); }
float quad( vec3 p, vec3 a, vec3 b, vec3 c, vec3 d )
{
    vec3 ba = b - a; vec3 pa = p - a;
    vec3 cb = c - b; vec3 pb = p - b;
    vec3 dc = d - c; vec3 pc = p - c;
    vec3 ad = a - d; vec3 pd = p - d;
    vec3 nor = cross( ba, ad );

    return sqrt(
    (sign(dot(cross(ba,nor),pa)) +
     sign(dot(cross(cb,nor),pb)) +
     sign(dot(cross(dc,nor),pc)) +
     sign(dot(cross(ad,nor),pd))<3.0)
     ?
     min( min( min(
     dot2(ba*clamp(dot(ba,pa)/dot2(ba),0.0,1.0)-pa),
     dot2(cb*clamp(dot(cb,pb)/dot2(cb),0.0,1.0)-pb) ),
     dot2(dc*clamp(dot(dc,pc)/dot2(dc),0.0,1.0)-pc) ),
     dot2(ad*clamp(dot(ad,pd)/dot2(ad),0.0,1.0)-pd) )
     :
     dot(nor,pa)*dot(nor,pa)/dot2(nor) );
}

float intersectSDF(float distA, float distB) {
    return max(distA, distB);
}

float unionSDF(float distA, float distB) {
    return min(distA, distB);
}

float differenceSDF(float distA, float distB) {
    return max(distA, -distB);
}

float smoothUnionSDF( float d1, float d2, float k ) {
    float h = clamp( 0.5 + 0.5*(d2-d1)/k, 0.0, 1.0 );
    return mix( d2, d1, h ) - k*h*(1.0-h); 
}

float sawtooth_wave(float x, float freq, float amplitude) {
  return (x * freq - floor(x * freq)) * amplitude - (amplitude / 2.0);
}

float triangle_wave(float x, float freq, float amplitude) {
  return abs(mod(x * freq, amplitude) - (0.5 * amplitude)) - (0.25 * amplitude);
}

float smoothblend(float d1, float d2, float a)
{
	return a * d1 + (1.0 - a) * d2;
}

mat4 rotateX(float theta) {
    float c = cos(theta);
    float s = sin(theta);

    return mat4(
        vec4(1, 0, 0, 0),
        vec4(0, c, -s, 0),
        vec4(0, s, c, 0),
        vec4(0, 0, 0, 1)
    );
}

mat4 rotateY(float theta) {
    float c = cos(theta);
    float s = sin(theta);

    return mat4(
        vec4(c, 0, s, 0),
        vec4(0, 1, 0, 0),
        vec4(-s, 0, c, 0),
        vec4(0, 0, 0, 1)
    );
}

mat4 rotateZ(float theta) {
    float c = cos(theta);
    float s = sin(theta);

    return mat4(
        vec4(c, -s, 0, 0),
        vec4(s, c, 0, 0),
        vec4(0, 0, 1, 0),
        vec4(0, 0, 0, 1)
    );
}

// noise (for wall)

float random (in vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

float noise (in vec2 st) {
    vec2 i = floor(st);
    vec2 f = fract(st);

    // Four corners in 2D of a tile
    float a = random(i);
    float b = random(i + vec2(1.0, 0.0));
    float c = random(i + vec2(0.0, 1.0));
    float d = random(i + vec2(1.0, 1.0));

    vec2 u = f * f * (3.0 - 2.0 * f);

    return mix(a, b, u.x) +
            (c - a)* u.y * (1.0 - u.x) +
            (d - b) * u.x * u.y;
}

#define OCTAVES 6
float fbm (in vec2 st) {
    // Initial values
    float value = 0.0;
    float amplitude = .5;
    float frequency = 0.;
    //
    // Loop of octaves
    for (int i = 0; i < OCTAVES; i++) {
        value += amplitude * noise(st);
        st *= 2.;
        amplitude *= .5;
    }
    return value;
}

// scene

float sceneSDF(vec3 pos, out int objHit) {
  vec3 tiltPos = (rotateY(-radians(10.0)) * vec4(pos, 1.0)).xyz;
  // tiltPos = (rotateX(radians(5.0)) * vec4(tiltpos, 1.0)).xyz;

  // surroundings
  // table
  float table = box(tiltPos + vec3(0.0, 6.77, 0.0), vec3(50.0, 1.0, 50.0));
  float tableInner = box(tiltPos + vec3(0.0, 6.75, 0.0), vec3(40.0, 1.0, 20.0));
  table = unionSDF(table, tableInner);
  // wall
  float disp = fbm(pos.xy);
  disp /= 8.0;
  float wall = box(tiltPos + vec3(0.0, 0.0, -25.0 + disp), vec3(50.0, 50.0, 1.0));
  float bg = unionSDF(table, wall);

  // switch
  // tablet
  vec3 spos = tiltPos + vec3(-2.0, 3.0, 0.0);
  // spos = pos;
  float tablet = roundBox(spos, vec3(4.5, 2.5, 0.2), 0.3);

  vec3 cpos = (rotateX(-radians(90.0)) * vec4(spos, 1.0)).xyz;
  float c1 = cylinder(cpos + vec3(4.6, 0.0, -2.7), vec2(0.2, 4.0));
  float c2 = cylinder(cpos - vec3(4.6, 0.0, 2.7), vec2(0.2, 4.0));
  float c = unionSDF(c1, c2);
  tablet = differenceSDF(tablet, c);

  // joycons
  float j1 = roundBox(spos + vec3(4.7, 0.0, 0.0), vec3(1.2, 2.5, 0.2), 0.3);
  float j2 = roundBox(spos + vec3(-4.7, 0.0, 0.0), vec3(1.2, 2.5, 0.2), 0.3);
  
  float rc = roundedCylinder(cpos + vec3(5.9, 0.0, 1.95), 0.425, 0.3, 0.2);
  j1 = unionSDF(j1, rc);
        rc = roundedCylinder(cpos + vec3(5.9, 0.0, -1.95), 0.425, 0.3, 0.2);
  j1 = unionSDF(j1, rc);
  float jj = roundBox(spos + vec3(5.96, 0.0, 0.0), vec3(0.5, 1.8, 0.2), 0.3);
  j1 = unionSDF(j1, jj);

  rc = roundedCylinder(cpos - vec3(5.9, 0.0, 1.95), 0.425, 0.3, 0.2);
  j2 = unionSDF(j2, rc);
        rc = roundedCylinder(cpos - vec3(5.9, 0.0, -1.95), 0.425, 0.3, 0.2);
  j2 = unionSDF(j2, rc);
  jj = roundBox(spos - vec3(5.96, 0.0, 0.0), vec3(0.5, 1.8, 0.2), 0.3);
  j2 = unionSDF(j2, jj);

  float jb = box(spos, vec3(4.8, 3.0, 1.0));
  j1 = differenceSDF(j1, jb);
  j2 = differenceSDF(j2, jb);

  float j = unionSDF(j1, j2);

  // directional buttons
  float rb1 = roundedCylinder(cpos + vec3(5.65, -0.3, -1.7), 0.1, 0.1, 0.3);
  float rb2 = roundedCylinder(cpos + vec3(5.65, -0.3, -0.9), 0.1, 0.1, 0.3);
  float rb3 = roundedCylinder(cpos + vec3(6.05, -0.3, -1.3), 0.1, 0.1, 0.3);
  float rb4 = roundedCylinder(cpos + vec3(5.25, -0.3, -1.3), 0.1, 0.1, 0.3);

  float rbut = unionSDF(unionSDF(unionSDF(rb1, rb2), rb3), rb4);

  float lb1 = roundedCylinder(cpos - vec3(5.65, 0.3, -0.6), 0.1, 0.1, 0.3);
  float lb2 = roundedCylinder(cpos - vec3(5.65, 0.3, 0.2), 0.1, 0.1, 0.3);
  float lb3 = roundedCylinder(cpos - vec3(6.05, 0.3, -0.2), 0.1, 0.1, 0.3);
  float lb4 = roundedCylinder(cpos - vec3(5.25, 0.3, -0.2), 0.1, 0.1, 0.3);

  float lbut = unionSDF(unionSDF(unionSDF(lb1, lb2), lb3), lb4);

  float but = unionSDF(rbut, lbut);
  j = unionSDF(j, but);

  // joysticks
  float js1 = roundedCylinder(cpos + vec3(5.65, -0.65, 0.2), 0.2, 0.1, 0.2);
  j = differenceSDF(j, js1);
  vec3 ccpos = (rotateX(radians(90.0)) * vec4(spos, 1.0)).xyz;
  js1 = cappedCone(ccpos + vec3(5.65, 0.7, -0.2), 0.2, 0.1, 0.2);
  float jt = torus(cpos + vec3(5.65, -0.6, 0.2), vec2(0.1, 0.3));
  float jtd = box(cpos + vec3(5.65, -0.55, 0.2), vec3(1.0, 0.1, 1.0));
  jt = differenceSDF(jt, jtd);
  js1 = unionSDF(js1, jt);
  j = unionSDF(j, js1);

  float js2 = roundedCylinder(cpos - vec3(5.65, 0.65, 1.4), 0.2, 0.1, 0.2);
  j = differenceSDF(j, js2);
  js2 = cappedCone(ccpos - vec3(5.65, -0.7, -1.4), 0.2, 0.1, 0.2);
  jt = torus(cpos - vec3(5.65, 0.6, 1.4), vec2(0.1, 0.3));
  jtd = box(cpos - vec3(5.65, 0.55, 1.4), vec3(1.0, 0.1, 1.0));
  jt = differenceSDF(jt, jtd);
  js2 = unionSDF(js2, jt);
  j = unionSDF(j, js2);

  // the other small buttons
  float p1 = box(spos + vec3(5.1, -2.2, 0.3), vec3(0.15, 0.05, 0.3));
  float p2 = box(spos + vec3(5.1, -2.2, 0.3), vec3(0.05, 0.15, 0.3));
  float buts = unionSDF(p1, p2);
  float m1 = box(spos + vec3(-5.1, -2.2, 0.3), vec3(0.15, 0.05, 0.3));
  buts = unionSDF(buts, m1);
  float sb1 = roundedCylinder(cpos + vec3(5.3, -0.45, 1.2), 0.1, 0.0, 0.1);
  buts = unionSDF(buts, sb1);
  float sb2 = box(spos + vec3(-5.3, 1.2, 0.45), vec3(0.15, 0.15, 0.1));
  buts = unionSDF(buts, sb2);
  j = unionSDF(j, buts);

  // shoulder triggers
  float st1 = roundedCylinder(cpos + vec3(5.9, 0.0, -2.05), 0.425, 0.1, 0.1);
  float st2 = roundedCylinder(cpos + vec3(-5.9, 0.0, -2.05), 0.425, 0.1, 0.1);
  float st = unionSDF(st1, st2);
  j = unionSDF(j, st);

  // screens
  float sc1 = roundBox(spos + vec3(0.0, 0.0, 0.205), vec3(4.3, 2.5, 0.1), 0.2);
  tablet = unionSDF(tablet, sc1);
  float sc2 = box(spos + vec3(0.0, 0.0, 0.41), vec3(3.7, 2.2, 0.1));
  tablet = unionSDF(tablet, sc2);

  // power and volume buttons
  float pw = roundedCylinder(spos + vec3(-3.5, -2.73, 0.0), 0.05, 0.0, 0.1);
  tablet = unionSDF(tablet, pw);
  float vl = roundBox(spos + vec3(-2.9, -2.45, 0.0), vec3(0.25, 0.3, 0.05), 0.1);
  float cm = cylinder(cpos + vec3(-2.9, 0.0, -2.90), vec2(0.1, 1.0));
  vl = differenceSDF(vl, cm);
  tablet = unionSDF(tablet, vl);

  // stand
  cpos = (rotateX(-radians(30.0)) * vec4(spos, 1.0)).xyz;
  float stn = box(cpos + vec3(2.5, 2.0, -0.3), vec3(0.4, 1.5, 0.05));
  tablet = unionSDF(tablet, stn);

  float swtch = unionSDF(j, tablet);

  // lamp
  float lamp = sphere(pos + vec3(-20.0, -15.0, -10.0), 3.0);
  float ln = cylinder(pos + vec3(-21.0, -5.0, -15.0), vec2(0.2, 15.0));
  float ls = cylinder(pos + vec3(-20.0, 5.0, -15.0), vec2(2.5, 1.0));
  float l = smoothUnionSDF(ls, ln, 0.8);
  lamp = unionSDF(lamp, l);

  // books
  vec3 tempPos = tiltPos + vec3(21.0, 3.0, -13.0);
  tempPos = (rotateY(radians(30.0)) * vec4(tempPos, 1.0)).xyz;
  float book1 = box(tempPos, vec3(8.0, 2.0, 6.0));
  float bm = box(tempPos, vec3(9.0, 1.7, 7.0));
  book1 = differenceSDF(book1, bm);
  float bpage = box(tempPos, vec3(7.5, 1.7, 5.8));
  bm = roundBox(tempPos - vec3(11.0, 0.0, 0.0), vec3(3.0, 0.8, 6.0), 1.0);
  bpage = differenceSDF(bpage, bm);
  book1 = unionSDF(book1, bpage);
  
  tempPos = (rotateY(radians(35.0)) * vec4(tempPos, 1.0)).xyz;
  tempPos += vec3(-2.0, -3.3, -3.0);
  float book2 = box(tempPos, vec3(7.0, 1.3, 6.0));
  bm = box(tempPos, vec3(9.0, 1.0, 7.0));
  book2 = differenceSDF(book2, bm);
  float bpage2 = box(tempPos, vec3(6.5, 1.0, 5.8));
  bm = roundBox(tempPos - vec3(11.0, 0.0, 0.0), vec3(4.0, 0.1, 6.0), 1.0);
  bpage2 = differenceSDF(bpage2, bm);
  book2 = unionSDF(book2, bpage2);

  float books = unionSDF(book1, book2);
  
  // switch shadow
  tempPos = (rotateY(radians(-30.0)) * vec4(tiltPos, 1.0)).xyz;
  float shadow = roundBox(tempPos + vec3(-1.4, 6.0, -6.5), vec3(5.0, 0.1, 4.0), 0.5);
  swtch = unionSDF(swtch, shadow);

  tempPos = tiltPos + vec3(0.0, 5.7, 3.0);
  tempPos = (rotateY(radians(-45.0)) * vec4(tempPos, 1.0)).xyz;
  float shadow2 = box(tempPos + vec3(0.0, 0.0, -0.9), vec3(9.0, 0.01, 5.0));
  swtch = unionSDF(swtch, shadow2);

  float ret = unionSDF(bg, swtch);
  ret = unionSDF(ret, lamp);
  ret = unionSDF(ret, books);

  if (ret == bg) {
    if (ret == table) {
      if (ret == tableInner) {
        objHit = 0;
      }
      else {
        objHit = 8;
      }
    }
    else if (ret == wall) {
      objHit = 7;
    }
  }
  else if (ret == swtch) {
    if (ret == j1) {
      objHit = 1;
    }
    else if (ret == j2) {
      objHit = 2;
    }
    else if (ret == sc1) {
      objHit = 3;
    }
    else if (ret == sc2) {
      objHit = 5;
    }
    else if (ret == shadow) {
      objHit = 9;
    }
    else if (ret == shadow2) {
      objHit = 14;
    }
    else {
      objHit = 4;
    }
  }
  else if (ret == lamp) {
    if (ret == ln || ret == ls) {
      objHit = 10;
    }
    else {
      objHit = 6;
    }
  }
  else if (ret == books) {
    if (ret == book1) {
      if (ret == bpage) {
        objHit = 13;
      }
      else {
        objHit = 11;
      }
    }
    else {
      if (ret == bpage2) {
        objHit = 13;
      }
      else {
        objHit = 12;
      }
    }
  }
  else {
    objHit = -1;
  }

  return ret;
}

float sceneSDF(vec3 pos) {
  int d;
  return sceneSDF(pos, d);
}

// ray marching
bool rayCubeIntersect(vec3 eye, vec3 dir, out float tnear, Bbox tcube) {
  // test for intersection of ray and scene
  tnear = -1000.0; 
  float tfar = 1000.0;

  for (int i = 0; i < 3; i++) {
    if (dir[i] == 0.0) {
      if (eye[i] < tcube.min[i] || eye[i] > tcube.max[i]) {
        return false;
      }
    }
    float t0 = (tcube.min[i] - eye[i]) / dir[i];
    float t1 = (tcube.max[i] - eye[i]) / dir[i];
    if (t0 > t1) {
      float temp = t0;
      t0 = t1;
      t1 = temp;
    }
    else if (t0 > tnear) {
      tnear = t0;
    }
    else if (t1 < tfar) {
      tfar = t1;
    }
  }

  if (tnear > tfar) return false;
  else return true;
}

bool traverseBVH(vec3 eye, vec3 dir, out float hit) {
  vec3 min = vec3(-5.45, 1.75, -0.3);
  vec3 max = vec3(1.45, 4.25, 0.3);
  vec3 rmin = (rotateY(-radians(10.0)) * vec4(min, 1.0)).xyz;
  vec3 rmax = (rotateY(-radians(10.0)) * vec4(max, 1.0)).xyz;
  Bbox switchBox = Bbox(rmin, rmax);
  if (rayCubeIntersect(eye, dir, hit, switchBox)) {
    min = vec3(-4.25, 1.75, -0.3);
    max = vec3(0.25, 4.25, 0.3);
    rmin = (rotateY(-radians(10.0)) * vec4(min, 1.0)).xyz;
    rmax = (rotateY(-radians(10.0)) * vec4(max, 1.0)).xyz;

    Bbox consoleBox = Bbox(rmin, rmax);
    if (rayCubeIntersect(eye, dir, hit, consoleBox)) {
      return true;
    }

    min = vec3(-5.45, 1.75, -0.3);
    max = vec3(-4.25, 4.25, 0.3);
    rmin = (rotateY(-radians(10.0)) * vec4(min, 1.0)).xyz;
    rmax = (rotateY(-radians(10.0)) * vec4(max, 1.0)).xyz;
    Bbox leftJCBox = Bbox(rmin, rmax);
    if (rayCubeIntersect(eye, dir, hit, leftJCBox)) {
      return true;
    }

    min = vec3(0.25, 1.75, -0.3);
    max = vec3(1.45, 4.25, 0.3);
    rmin = (rotateY(-radians(10.0)) * vec4(min, 1.0)).xyz;
    rmax = (rotateY(-radians(10.0)) * vec4(max, 1.0)).xyz;
    Bbox rightJCBox = Bbox(rmin, rmax);
    if (rayCubeIntersect(eye, dir, hit, rightJCBox)) {
      return true;
    }
    return true;
  }
  else {
    if (dir.y < 1.75) { // hit background
      return false;
    }
    else { // hit table (not sure how to check this)
      hit = 0.0;
      return true;
    }
  }
}

float raymarch(vec3 eye, vec3 dir, float start, float end, out int obj) {
  float depth = start;
  // float t;
  // bool hitScene = traverseBVH(eye, dir, t);
  // if (!hitScene) return end;
  // depth = t;

  for (int i = 0; i < MAX_RAY_STEPS; i++) {
    float dist = sceneSDF(eye + depth * dir, obj);
    if (dist < EPSILON) {
      return depth;
    }
    depth += dist;
    if (depth >= end) {
      return end;
    }
  }
  return end;
}

// noise stuff

// perlin noise

vec2 falloff(vec2 t) {
  return t*t*t*(t*(t*6.0-15.0)+10.0);
  }
vec4 permute(vec4 x){return mod(((x*34.0)+1.0)*x, 289.0);}
vec3 permute(vec3 x) {
    return mod((34.0 * x + 1.0) * x, 289.0);
  }

float pnoise(vec2 P){
  vec4 Pi = floor(P.xyxy) + vec4(0.0, 0.0, 1.0, 1.0);
  vec4 Pf = fract(P.xyxy) - vec4(0.0, 0.0, 1.0, 1.0);
  Pi = mod(Pi, 289.0); // To avoid truncation effects in permutation
  vec4 ix = Pi.xzxz;
  vec4 iy = Pi.yyww;
  vec4 fx = Pf.xzxz;
  vec4 fy = Pf.yyww;
  vec4 i = permute(permute(ix) + iy);
  vec4 gx = 2.0 * fract(i * 0.0243902439) - 1.0; // 1/41 = 0.024...
  vec4 gy = abs(gx) - 0.5;
  vec4 tx = floor(gx + 0.5);
  gx = gx - tx;
  vec2 g00 = vec2(gx.x,gy.x);
  vec2 g10 = vec2(gx.y,gy.y);
  vec2 g01 = vec2(gx.z,gy.z);
  vec2 g11 = vec2(gx.w,gy.w);
  vec4 norm = 1.79284291400159 - 0.85373472095314 * 
    vec4(dot(g00, g00), dot(g01, g01), dot(g10, g10), dot(g11, g11));
  g00 *= norm.x;
  g01 *= norm.y;
  g10 *= norm.z;
  g11 *= norm.w;
  float n00 = dot(g00, vec2(fx.x, fy.x));
  float n10 = dot(g10, vec2(fx.y, fy.y));
  float n01 = dot(g01, vec2(fx.z, fy.z));
  float n11 = dot(g11, vec2(fx.w, fy.w));
  vec2 fade_xy = falloff(Pf.xy);
  vec2 n_x = mix(vec2(n00, n01), vec2(n10, n11), fade_xy.x);
  float n_xy = mix(n_x.x, n_x.y, fade_xy.y);
  return 2.3 * n_xy;
}

void generateWoodPattern(vec3 pos, out float off) {  
  off = pnoise(vec2(pos.x, pos.y)) + pnoise(vec2(pos.x, pos.z));
  off *= 2.0;
  off = pnoise(vec2(pos.x, off)) + pnoise(vec2(pos.x, pos.y));
}

// lighting/color stuff

vec3 estimateNormal(vec3 p) {
    return normalize(vec3(
        sceneSDF(vec3(p.x + EPSILON, p.y, p.z)) - sceneSDF(vec3(p.x - EPSILON, p.y, p.z)),
        sceneSDF(vec3(p.x, p.y + EPSILON, p.z)) - sceneSDF(vec3(p.x, p.y - EPSILON, p.z)),
        sceneSDF(vec3(p.x, p.y, p.z + EPSILON)) - sceneSDF(vec3(p.x, p.y, p.z - EPSILON))
    ));
}

vec3 phongContribForLight(vec3 k_d, vec3 k_s, float alpha, vec3 p, vec3 lightPos, vec3 lightIntensity, out vec3 nor) {
  vec3 N = estimateNormal(p);
  nor = N;
  vec3 L = normalize(lightPos - p);
  vec3 V = normalize(u_Eye - p);
  vec3 R = normalize(reflect(-L, N));

  float dotLN = dot(L, N);
  float dotRV = dot(R, V);

  if (dotLN < 0.0) {
    return vec3(0.0, 0.0, 0.0);
  }
  if (dotRV < 0.0) {
    return lightIntensity * (k_d * dotLN);
  }
  return lightIntensity * (k_d * dotLN + k_s * pow(dotRV, alpha));
}

vec3 phongIllumination(vec3 k_a, vec3 k_d, vec3 k_s, float alpha, vec3 p, out vec3 nor) {
  const vec3 ambientLight = 0.5 * vec3(1.0, 1.0, 1.0);
  vec3 color = ambientLight * k_a;

  vec3 light1Intensity = vec3(0.4, 0.4, 0.4);
  color += phongContribForLight(k_d, k_s, alpha, p, light1Pos, light1Intensity, nor);

  vec3 light2Intensity = vec3(0.4, 0.4, 0.4);
  // color += phongContribForLight(k_d, k_s, alpha, p, light2Pos, light2Intensity, nor);
  
  return color;
}

vec3 calcColor(vec3 p, int obj) {
  vec3 K_a = vec3(0.1, 0.1, 0.1);
  vec3 K_d = vec3(0.0);
  if (obj == 0 || obj == 14) { // table 
    float scale;
    generateWoodPattern(p, scale);
    vec3 col1 = vec3(252.0, 205.0, 80.0) / 255.0;
    vec3 col2 = vec3(140.0, 104.0, 32.0) / 255.0;
    K_d = mix(col1, col2, scale);

    if (obj == 14) {
      return K_d / 6.0;
    }
  }
  else if (obj == 1) { // red
    K_d = vec3(229.0, 86.0, 64.0) / 255.0;
  }
  else if (obj == 2) { // blue
    K_d = vec3(64.0, 200.0, 234.0) / 255.0;
  }
  else if (obj == 3) { // black
    K_d = vec3(0.0, 0.0, 0.0) / 255.0;
  }
  else if (obj == 4) { // gray
    K_d = vec3(89.0, 89.0, 89.0) / 255.0;
  }
  else if (obj == 5) { // screen
    // return vec3(1.0 * sin(u_Time / 10.0), 1.0 * sin(u_Time / 18.0), 1.0 * cos(u_Time / 15.0));
    
    vec2 st = p.xy;

    st *= 5.0; // Scale the coordinate system
    vec2 ipos = floor(st);  // get the integer coords
    vec2 fpos = fract(st);  // get the fractional coords

    // Assign a random value based on the integer coord
    vec3 color = vec3(random( ipos ) * (sin(u_Time * fpos.x)));
    return color;
  }
  else if (obj == 6) { // lamp light
    return vec3(1.0);
  }
  else if (obj == 7) { // wall
    K_d = vec3(239.0, 224.0, 184.0) / 255.0;
  }
  else if (obj == 8) { // table border
    K_d = vec3(119.0, 87.0, 28.0) / 255.0;
  }
  else if (obj == 9) { // shadow
    return vec3(0.0);
  }
  else if (obj == 10) { // lamp body
    K_d = vec3(1.0);
  }
  else if (obj == 11) { // green
    K_d = vec3(45.0, 183.0, 0.0) / 255.0;
  }
  else if (obj == 12) { // dark blue
    K_d = vec3(0.0, 78.0, 181.0) / 255.0;
  }
  else if (obj == 13) { // pages
    // K_d = vec3(255.0, 253.0, 214.0) / 255.0;

    float y = p.y * 5.0;
    if (abs(y - floor(y)) < 0.05) {
      K_d = vec3(99.0, 99.0, 99.0) / 255.0;
    }
    else {
      K_d = vec3(255.0, 253.0, 214.0) / 255.0;
    }
  }
  vec3 K_s = vec3(1.0, 1.0, 1.0);
  float shininess = 1.0;
  if (obj == 10) {
    shininess = 10.0;
  }

  vec3 nor;
  vec3 color = phongIllumination(K_a, K_d, K_s, shininess, p, nor);

  // other lights

  vec3 lampDir = normalize(light2Pos);
  float lampDot = max(0.0, dot(nor, lampDir));
  if (obj > 5 || obj == 0) {
    color += lampDot * vec3(0.75, 0.66, 0.52);
  }
  
  return color;
}

float shadow( vec3 ro, vec3 rd, float mint, float maxt )
{
    for( float t=mint; t < maxt; )
    {
        float h = sceneSDF(ro + rd*t);
        if( h<0.001 )
            return 0.0;
        t += h;
    }
    return 1.0;
}

float softshadow( vec3 ro, vec3 rd, float mint, float maxt, float k )
{
    float res = 1.0;
    for( float t=mint; t < maxt; )
    {
        float h = sceneSDF(ro + rd*t);
        if( h<0.001 )
            return 0.0;
        res = min( res, k*h/t );
        t += h;
    }
    return res;
}

float calcShadow(vec3 pos) {

  // return shadow(pos, light1Pos, 0.1, 255.0);
  return softshadow(pos, vec3(20.0, 15.0, 10.0), 0.15, 255.0, 2.0);
}

// ray casting stuff

vec3 raycast(out vec3 screenpos) {
  float a = fs_Pos.x;
  float b = fs_Pos.y;
  float FOVY = radians(90.0);
  float aspect = u_Dimensions.x / u_Dimensions.y;
  float len = length(u_Ref - u_Eye);

  vec3 R = normalize(cross(u_Ref - u_Eye, u_Up));
  vec3 V = tan(FOVY / 2.0) * len * u_Up;
  vec3 H = aspect * len * tan(FOVY / 2.0) * R;

  screenpos = u_Ref + a * H + b * V;
  vec3 ray_dir = screenpos - u_Eye;

  return normalize(ray_dir);
}

vec3 worldToScreen(vec3 wp) {
  mat4 view = mat4(vec4(1.0, 0.0, 0.0, -1.0 * u_Eye.x),
                   vec4(0.0, 1.0, 0.0, -1.0 * u_Eye.y),
                   vec4(0.0, 0.0, 1.0, -1.0 * u_Eye.z),
                   vec4(0.0, 0.0, 0.0, 1));

  float FOVY = radians(90.0);
  float aspect = u_Dimensions.x / u_Dimensions.y;
  mat4 proj = mat4(vec4(1.0 / (aspect * tan(FOVY / 2.0)), 0.0, 0.0, 0.0),
                   vec4(0.0, 1.0 / tan(FOVY / 2.0), 0.0, 0.0),
                   vec4(0.0, 0.0, 1.0, -1.0),
                   vec4(0.0, 0.0, 1.0, 0.0));
  
  vec4 s = view * proj * vec4(wp, 1.0);

  float ndcsx = (2.0 * s.x / u_Dimensions.x) - 1.0;
  float ndcsy = (1.0 - (2.0 * s.y / u_Dimensions.y));
  return vec3(ndcsx, ndcsy, 1.0);
  // return vec3(s.x, s.y, s.z);
}

void main() {
  vec3 screenpos = vec3(0.0);
  vec3 ray_dir = raycast(screenpos);
  int obj;
  float dist = raymarch(u_Eye, ray_dir, MIN_DIST, MAX_DIST, obj);

  // if (dist > MAX_DIST - EPSILON) {
    // out_Col = vec4(vec3(221.0, 204.0, 146.0) / 255.0, 1.0);
  //   out_Col = vec4((0.5 * (ray_dir + vec3(1.0, 1.0, 1.0))), 1.0);
    // return;
  // }

  vec3 p = u_Eye + dist * ray_dir;
  
  vec3 color = calcColor(p, obj);
  // shadow - doesn't work
  // float shadow = calcShadow(p);
  // color *= shadow;

  // distance fog
  float fogT = smoothstep(35.0, 60.0, distance(p, u_Eye));
  vec3 fog_col = mix(color, vec3(56.0, 46.0, 26.0) / 255.0, fogT); 
  color = fog_col;

  // lamp glow
  vec3 screenLamp = worldToScreen(vec3(0.0, 0.0, 1.0));
  screenLamp += vec3(11.0, 7.0, 0.0);
  float d = distance(screenpos.xy, screenLamp.xy);
  if (d < 10.0) {
    float lampT = smoothstep(1.0, 10.0, d);
    vec3 new_col = mix(vec3(1.0), color, lampT);

    color = new_col;

    if (d < 2.0) {
      lampT = smoothstep(1.0, 2.0, d);
      new_col = mix(vec3(1.0), color, lampT);
      color += new_col;
    }
  }
  

  out_Col = vec4(color, 1.0);
}
