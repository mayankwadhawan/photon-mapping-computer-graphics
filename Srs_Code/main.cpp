#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <vector>
#include <algorithm>
#include <pthread.h>
#include <cstdlib>
#include <GL/glut.h>
#include <math.h>
#include <iostream>

#define WIDTH 512
#define HIEGHT 512
#define XCOOR 100
#define YCOOR 50
#define SPHERE   0
#define PLANE    1
#define NO_EFFECT 0
#define REFLECTION 1
#define REFRACTION 2
#define NO_INTERSECTION (1.0e6)

#ifndef __COORDINATE_VECTOR3_H__
#define __COORDINATE_VECTOR3_H__

namespace VectorSpace {

class CoordinateVector3 {
public:
    CoordinateVector3()
        : m_x(0.0)
        , m_y(0.0)
        , m_z(0.0)
    {
    }

    CoordinateVector3(double x, double y, double z)
        : m_x(x)
        , m_y(y)
        , m_z(z)
    {
    }

    CoordinateVector3(const float p[3])
        : m_x(p[0])
        , m_y(p[1])
        , m_z(p[2])
    {
    }

    CoordinateVector3(const double p[3])
        : m_x(p[0])
        , m_y(p[1])
        , m_z(p[2])
    {
    }

    CoordinateVector3(const CoordinateVector3 &v)
        : m_x(v.x())
        , m_y(v.y())
        , m_z(v.z())
    {
    }

    CoordinateVector3 &operator =(const CoordinateVector3 &v)
    {
        m_x = v.x();
        m_y = v.y();
        m_z = v.z();
        return (*this);
    }

    const double operator [](int caseVar) const
    {
      switch(caseVar) {
        case 0: return m_x;
        case 1: return m_y;
        case 2: return m_z;
        default: return 0.0;
      }
    }
    double operator [](int idx)
    {
      switch(idx) {
        case 0: return m_x;
        case 1: return m_y;
        case 2: return m_z;
        default: return 0.0;
      }
    }

    double abs() const
    {
        return sqrt(m_x * m_x + m_y * m_y + m_z * m_z);
    }

    bool isZero() const
    {
        return !m_x && !m_y && !m_z;
    }

    void normalize()
    {
        double absValue = abs();
        if (!absValue)
            return;

        double k = 1.0 / absValue;
        m_x *= k;
        m_y *= k;
        m_z *= k;
    }

    double x() const { return m_x; }
    double y() const { return m_y; }
    double z() const { return m_z; }

private:
    double m_x;
    double m_y;
    double m_z;
};

inline CoordinateVector3 operator+(const CoordinateVector3& v1, const CoordinateVector3& v2)
{
    return CoordinateVector3(v1.x() + v2.x(), v1.y() + v2.y(), v1.z() + v2.z());
}

inline CoordinateVector3 operator-(const CoordinateVector3& v1, const CoordinateVector3& v2)
{
    return CoordinateVector3(v1.x() - v2.x(), v1.y() - v2.y(), v1.z() - v2.z());
}

inline CoordinateVector3 operator*(double k, const CoordinateVector3& v)
{
    return CoordinateVector3(k * v.x(), k * v.y(), k * v.z());
}

inline CoordinateVector3 operator*(const CoordinateVector3& v, double k)
{
    return CoordinateVector3(k * v.x(), k * v.y(), k * v.z());
}

static inline double dot(const CoordinateVector3& v1, const CoordinateVector3& v2)
{
    return v1.x() * v2.x() + v1.y() * v2.y() + v1.z() * v2.z();
}

inline CoordinateVector3 cross(const CoordinateVector3& v1, const CoordinateVector3& v2)
{
    double x3 = v1.y() * v2.z() - v1.z() * v2.y();
    double y3 = v1.z() * v2.x() - v1.x() * v2.z();
    double z3 = v1.x() * v2.y() - v1.y() * v2.x();
    return CoordinateVector3(x3, y3, z3);
}

static inline double distance(const CoordinateVector3& v1, const CoordinateVector3& v2)
{
    return (v1 - v2).abs();
}

} 

#endif 

using std::vector;
using std::max;
using std::min;
using VectorSpace::CoordinateVector3;

static const CoordinateVector3 gOrigin;
static       CoordinateVector3 Light(0.0,1.2,3.75);   //Point Light Source Position
static const int     reflection_limit = 4;


template <typename T> inline T
constrain(T src, T lower, T upper) { return min(upper, max(src, lower)); }
inline bool checkerOdd(int x) { return x & 1; }

int sizeofImage = 512;            //rendering screen size
int surfaceType = 2;            //surface tpye = 0:SPHERE, 1:PLANE
int numOfSurfcae = 0;          //num of surface

// ----- Photon Mapping -----
int   totalPhotonEmission = 6000;     //Number of Photons Emitted
int   totalIterations = 3;        //Number of Times Each Photon Bounces
bool  photonMapping = true;  //Enable Photon Lighting?
float brightness = 100.0;     //Number of Photons Integrated at Brightest Pixel
int   photonPerSurface[64];       //Photon Count for Each Scene surface

//  Allocated Memory for Per-Object Photon Info
//  0 : location
//  1 : direction
//  2 : energy
CoordinateVector3 photons[64][5000][3];

class Surfaces {
 
  public :
    Surfaces(int type, int idx, float *cod);
    int    getSurface()   { return surfaceType; }
    int    getSpecialEffect() { return specialEffect; }
    int    getIndex()  { return index; }
    float  getRefractive()  { return refractive; }
    void   setSpecialEffect(int op)   { specialEffect = op; }
    void   setRefractive(float rf)   { this->refractive = rf; }
    void   setColor(const float *cl) {
      for(int i=0; i<3; i++) color[i] = cl[i];
    }
    float coords[9];
    float color[3];

    CoordinateVector3 sphereNormalCalculation(const CoordinateVector3 &P, const CoordinateVector3 &O);
    CoordinateVector3 planeNormalCalculation(const CoordinateVector3 &P, const CoordinateVector3 &O);

    double  calcSphereIntersection(const CoordinateVector3 &ray, const CoordinateVector3 &org);
    double  calcPlaneIntersection(const CoordinateVector3 &ray, const CoordinateVector3 &org);

 
  private :
    int surfaceType;
    int specialEffect;
    int index;
    float refractive;
};

Surfaces::Surfaces(int tp, int idx, float *cod) {

  index = idx;
  specialEffect = NO_EFFECT;
  refractive = 1.0;
  surfaceType = tp;

  int ncods = 0;
  switch(tp) {
  case SPHERE : ncods=4; break;
  case PLANE  : ncods=2; break;
  }
  for(int i=0; i<ncods; i++) coords[i] = cod[i];
  for(int i=0; i<3; i++) color[i] = 1.0;
}


CoordinateVector3 Surfaces::sphereNormalCalculation(const CoordinateVector3 &P, const CoordinateVector3 &O)
{

  CoordinateVector3 center(coords[0], coords[1], coords[2]);
  CoordinateVector3 ans = P - center;
  ans.normalize();
  return ans;
}

CoordinateVector3 Surfaces::planeNormalCalculation(const CoordinateVector3 &P, const CoordinateVector3 &O)
{
  int axis = (int) this->coords[0];
  float N[3] = {0.0,0.0,0.0};
  N[axis] = O[axis] - this->coords[1];      //Vector From Surface to Light
  CoordinateVector3 ans(N);
  ans.normalize();
  return ans;
}


double Surfaces::calcSphereIntersection(const CoordinateVector3 &r, const CoordinateVector3 &o) //Ray-Sphere Intersection: r=Ray Direction, o=Ray Origin
{
  //  s = Sphere Center Translated into Coordinate Frame of Ray Origin
  CoordinateVector3 center(coords[0], coords[1], coords[2]);
  CoordinateVector3 s = center - o;

  //radius=Sphere Radius
  float radius = coords[3];

  //Intersection of Sphere and Line
  float A = dot(r,r);
  float B = -2.0 * dot(s,r);
  float C = dot(s,s) - radius * radius;
  float D = B * B - 4 * A * C;

  if (D > 0.0) {
    float sign = (C < -0.00001) ? 1 : -1;
    return (-B + sign*sqrt(D))/(2*A);
  }
  return NO_INTERSECTION;
}

double Surfaces::calcPlaneIntersection(const CoordinateVector3 &r, const CoordinateVector3 &o)
{
  int axis = (int) coords[0];            //Determine Orientation of Axis-Aligned Plane
  if (r[axis] != 0.0){                        //Parallel Ray -> No Intersection
    return  (coords[1] - o[axis]) / r[axis]; //Solve Linear Equation (rx = p-o)
  }
  return NO_INTERSECTION;
}


/**functions**/

CoordinateVector3 reflection(
    Surfaces *ob,
    const CoordinateVector3 &point,
    const CoordinateVector3 &ray,
    const CoordinateVector3 &fromPoint);
CoordinateVector3 refraction(
    Surfaces *ob,
    const CoordinateVector3 &point,
    const CoordinateVector3 &ray,
    const CoordinateVector3 &fromPoint,
    float &ref);

CoordinateVector3 photonCollection(const CoordinateVector3 &p, Surfaces *ob);
void    emitPhotons();
void    photonArrays(Surfaces *ob,
    const CoordinateVector3 &location,
    const CoordinateVector3 &direction,
    const CoordinateVector3 &energy );
void    photonShadowCalculation(const CoordinateVector3 &ray, const CoordinateVector3 &pnt);
void    drawPhoton(const CoordinateVector3 &rgb, const CoordinateVector3 &p);


CoordinateVector3 mulColor(const CoordinateVector3 &rgbIn, Surfaces *ob);

void render();
void resetRender();

void surfaceInitialization();
void deleteSurfaces();

void display();
void resize (int w, int h);
void onKeyPress(unsigned char key, int x, int y);
void onClick(int button, int state, int x, int y);
void onDrag(int x, int y);
void onTimer(int val);

//screen status variables
//to stop drawing
bool empty = true;
//to switch Views

bool mouseDragging = false;
int  mouseX, mouseY;

//rendering pixel info
int  rowCounter, colCounter, photonCounter, resolutionMax;



std::vector<Surfaces*> surfaces;

typedef struct SurfaceIntersection {
  Surfaces    *surf;
  double  dist;
  SurfaceIntersection() {
    dist = NO_INTERSECTION;
    surf = NULL;
  }
} SurfaceIntersection;


int main(int argc, char *argv[]) {

  surfaceInitialization();

  emitPhotons();
  resetRender();

  glutInit(&argc,argv);
  glutInitWindowPosition(XCOOR, YCOOR);
  glutInitWindowSize(WIDTH, HIEGHT);
  glutInitDisplayMode(GLUT_RGBA | GLUT_DEPTH);

  glutCreateWindow("GLUT PHOTON MAPPING");

  glutDisplayFunc(display);
  glutTimerFunc(0, onTimer, 0);
  glutReshapeFunc(resize);

  glutKeyboardFunc(onKeyPress);
  glutMouseFunc(onClick);
  glutMotionFunc(onDrag);

  atexit(deleteSurfaces);

  glClear(GL_COLOR_BUFFER_BIT);
  glutMainLoop();

  return 1;
}

/*
Ray-Geometry Intersections
*/

double rayLight(Surfaces *ob, const CoordinateVector3 &r, const CoordinateVector3 &o){

  int tp = ob->getSurface();
  //switch intersection func with object type
  if      (tp == SPHERE) {
    return ob->calcSphereIntersection(r, o);
  } else if (tp == PLANE) {
    return ob->calcPlaneIntersection(r, o);
  }

  return NO_INTERSECTION;
}

float diffusionLight(const CoordinateVector3 &N, const CoordinateVector3 &P)
{
  //  Diffuse Lighting at Point P with Surface Normal N
  CoordinateVector3 L = Light - P;
  L.normalize();
  return dot(N,L);
}

CoordinateVector3 surfaceNormal(Surfaces *ob, const CoordinateVector3 &P, const CoordinateVector3 &Inside){
  if (ob->getSurface() == SPHERE)     {
    return ob->sphereNormalCalculation(P, Inside);
  } else if (ob->getSurface() == PLANE) {
    return ob->planeNormalCalculation(P, Inside);
  }
  return CoordinateVector3();
}

float light(Surfaces *ob, const CoordinateVector3 &P, float lightAmbient){
  CoordinateVector3 N = surfaceNormal(ob, P, Light);
  float   i = diffusionLight(N, P);
  //add in ambient light by constraining min value
  return min(1.0f, max(i, lightAmbient));
}

/*
Raytracing
*/

SurfaceIntersection rayTracing(const CoordinateVector3 &ray, const CoordinateVector3 &origin)
{
  //init intersection status
  SurfaceIntersection istat;

  //check intersection for each object
  for (int i=0; i<numOfSurfcae; i++) {
    double dist = rayLight(surfaces[i], ray, origin);
    if(dist < istat.dist && dist > 1.0e-5) {
      istat.dist = dist;
      istat.surf  = surfaces[i];
    }
  }
  return istat;
}

CoordinateVector3 colorPixelCalculation(float x, float y){
  CoordinateVector3 rgb(0.0,0.0,0.0);

  //generate Ray for each pixel
  //Convert Pixels to Image Plane Coordinates
  CoordinateVector3 ray(
      x / sizeofImage - 0.5 ,
    -(y / sizeofImage - 0.5),
    1.0
    //Focal Length = 1.0
  );

  float refractive = 1.0;
  CoordinateVector3 from = gOrigin;

  SurfaceIntersection istat = rayTracing(ray, from);
  if (istat.dist >= NO_INTERSECTION){ return rgb; }

  //get point of intersection
  CoordinateVector3 pnt = from + ray * istat.dist;

  int ref = 0;
  //  Mirror Surface on This Specific Object
  while (istat.surf->getSpecialEffect() != NO_EFFECT && ref < reflection_limit){
    if(istat.surf->getSpecialEffect() == REFLECTION) { ray = reflection(istat.surf, pnt, ray, from); }
    else                       /*OPT_REFRACT*/{ ray = refraction(istat.surf, pnt, ray, from, refractive); }
    ref++;

    from = pnt;
    istat = rayTracing(ray, from);             //Follow the Reflected Ray
    if (istat.dist >= NO_INTERSECTION){ return rgb; }
    else {
      pnt = from + ray * istat.dist;
    }
  }

  if (photonMapping){
    //Lighting via Photon Mapping
    rgb = photonCollection(pnt, istat.surf);
  } else {
    //Lighting via Standard Illumination Model (Diffuse + Ambient)
    //Remember Intersected Object
    SurfaceIntersection org_stat = istat;

    //If in Shadow, Use Ambient Color of Original Object
    static const float ambient = 0.1;

    //Raytrace from Light to Object
    SurfaceIntersection lht_stat = rayTracing(pnt - Light, Light);

    float intensity = ambient;
    if (lht_stat.surf == org_stat.surf) {
      //Ray from Light -> Surface Hits surface First? : not in shadow
      intensity = light(lht_stat.surf, pnt, ambient);
    }

    CoordinateVector3 energy(intensity, intensity, intensity);
    rgb = mulColor(energy, lht_stat.surf);
  }
  return rgb;
}

CoordinateVector3 reflection(
    Surfaces *ob,
    const CoordinateVector3 &point,
    const CoordinateVector3 &ray,
    const CoordinateVector3 &from)
{
  CoordinateVector3 N = surfaceNormal(ob, point, from);

  CoordinateVector3 ans = ray - N * (2 * dot(ray,N));
  ans.normalize();
  return ans;
}

CoordinateVector3 refraction(
    Surfaces *ob,
    const CoordinateVector3 &point,
    const CoordinateVector3 &ray,
    const CoordinateVector3 &from,
    float &ref)
{
  CoordinateVector3 N = surfaceNormal(ob, point, from);

  float n1 = ref;
  float n2 = ob->getRefractive();
  float s  = dot(ray, N);

  if(ob->getSurface() == SPHERE && s > 0) {
    //from inside to outside : swap n1 and n2
    float tmp = n1;
    n1 = n2;
    n2 = tmp;
  }

  float n  = n1 / n2;

  CoordinateVector3 ans = n * (ray - s * N) - N * sqrt(1 - n * n * (1 - s * s) );
  ans.normalize();

  ref = n2;
  return ans;
}

/*
Photon Mapping
*/

CoordinateVector3 photonCollection(const CoordinateVector3 &p, Surfaces *ob)
{
  //Photon Integration Area (Squared for Efficiency)
  static const float sqRadius = 0.7;

  CoordinateVector3 energy;
  int id = ob->getIndex();

  CoordinateVector3 N = surfaceNormal(ob, p, gOrigin);

  for (int i = 0; i < photonPerSurface[id]; i++) {
    //Photons Which Hit Current Object
    double cur_dist = distance(p, photons[id][i][0]);

    //Is Photon Close to Point?
    if (cur_dist < sqRadius) {
      float weight = max(0.0, -dot(N, photons[id][i][1]) );

      //Single Photon Diffuse Lighting
      //Weight by Photon-Point Distance
      weight     *= (1.0 - cur_dist) / brightness;
      CoordinateVector3 tmp = photons[id][i][2] * weight;
      energy      = energy + tmp;
    }
  }
  return energy;
}

CoordinateVector3 directionRandamize(double s)
{
  //generate vector with random derection
  double tmp[3];
  for(int i=0; i<3; i++) {
    tmp[i] = (double)rand() * 2 * s / RAND_MAX - s;
  }
  CoordinateVector3 ans(tmp);
  ans.normalize();
  return ans;
}

void emitPhotons(){

  //"randomized" photons are generated with the same properties indeed
  srand(0);

  //init photon num
  for (int t = 0; t < numOfSurfcae; t++) { photonPerSurface[t] = 0; }

  CoordinateVector3 rgb, ray, col;
  CoordinateVector3 white(1.0, 1.0, 1.0);

  //control photon num with rendering option
  const int num_photon =  totalPhotonEmission;
  for (int i = 0; i < num_photon; i++){
    int bounces = 1;

    //initialize photon properties (color, direction, location)
    rgb = white;
    ray = directionRandamize(1.0);
    CoordinateVector3 from = Light;

    //randomize photon locations
    while (from.y() >= Light.y()) {
      //+Y dir
      from = directionRandamize(1.0) * 0.75 + Light;
    }

    //photons outside of the room : invalid
    if (fabs(from.x()) > 1.5 || fabs(from.y()) > 1.2 ) {
      bounces = totalIterations + 1;
    }

    //photons inside any objects : invalid
    for(int dx = 0; dx<numOfSurfcae; dx++) {
      Surfaces *ob = surfaces[dx];

      if(ob->getSurface() != SPHERE) continue;

      CoordinateVector3 center(ob->coords);
      if(distance(from, center) < ob->coords[3]) {
        bounces = totalIterations+1;
      }
    }

    //calc intersection (1st time)
    float refractive = 1.0;
    SurfaceIntersection istat = rayTracing(ray, from);

    //calc bounced photon's intercection (2nd, 3rd, ...)
    while (istat.dist < NO_INTERSECTION && bounces <= totalIterations){
      CoordinateVector3 pnt = from + ray * istat.dist;

      //reflect or refract
      int ref = 0;
      while (istat.surf->getSpecialEffect() != NO_EFFECT && ref < reflection_limit){
        if(istat.surf->getSpecialEffect() == REFLECTION) { ray = reflection(istat.surf, pnt, ray, from); }
        else                       /*OPT_REFRACT*/{ ray = refraction(istat.surf, pnt, ray, from, refractive); }
        ref++;

        from = pnt;
        istat = rayTracing(ray, from);             //Follow the Reflected Ray
        if (istat.dist >= NO_INTERSECTION){ break; }
        else {
          pnt = from + ray * istat.dist;
        }
      }

      if(istat.dist >= NO_INTERSECTION) { continue; }

      col = mulColor(rgb, istat.surf);
      rgb = col * (1.0 / sqrt((double)bounces));

      photonArrays(istat.surf, pnt, ray, rgb);

      drawPhoton(rgb, pnt);
      photonShadowCalculation(ray, pnt);

      ray = reflection(istat.surf, pnt, ray, from);

      istat = rayTracing(ray, pnt);
      if(istat.dist >= NO_INTERSECTION){ break; }

      from = pnt;
      bounces++;
    }
  }
}

void photonArrays(Surfaces *ob, const CoordinateVector3 &location, const CoordinateVector3 &direction, const CoordinateVector3 &energy){
  int id = ob->getIndex();
  //  0 : location
  //  1 : direction
  //  2 : energy
  photons[id][photonPerSurface[id]][0] = location;
  photons[id][photonPerSurface[id]][1] = direction;
  photons[id][photonPerSurface[id]][2] = energy;
  photonPerSurface[id]++;
}

void photonShadowCalculation(const CoordinateVector3 &ray, const CoordinateVector3 &pnt){
  CoordinateVector3 shadow (-0.25,-0.25,-0.25);

  //Start Just Beyond Last Intersection
  CoordinateVector3 bumpedPoint = pnt + ray * 1.0e-5;

  //Trace to Next Intersection (In Shadow)
  SurfaceIntersection istat = rayTracing(ray, bumpedPoint);
  if(istat.dist >= NO_INTERSECTION) { return; }

  //3D Point
  CoordinateVector3 shadowPoint = bumpedPoint + ray * istat.dist;

  photonArrays(istat.surf, shadowPoint, ray, shadow);
}

CoordinateVector3 mulColor(const CoordinateVector3 &rgbIn, Surfaces *ob)
{
  //Specifies Material Color of Each Object
  return CoordinateVector3(
      ob->color[0] * rgbIn[0],
      ob->color[1] * rgbIn[1],
      ob->color[2] * rgbIn[2] );
}

void resize(int w, int h) {
  //set orthogonal view
  glViewport(0, 0, WIDTH, HIEGHT);
  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();
  glOrtho(0.0,(double)WIDTH,0.0,(double)HIEGHT,-10.0,10.0);
}

void display(){
  
    if (empty) render();
  //  else sleep(1);  //Only Draw if Image Not Fully Rendered
  
  glFlush();
}

void render(){ //Render Several Lines of Pixels at Once Before Drawing
  int x,y,iterations = 0;
  CoordinateVector3 rgb;
  resolutionMax=512;
  while (iterations < (mouseDragging ? 1024 : max(resolutionMax, 256) )){
    if (colCounter >= resolutionMax) {
      rowCounter++;
      colCounter = 0;

      if (rowCounter >= resolutionMax) {
        photonCounter++;
        rowCounter = 0;
        resolutionMax = int(pow(2.0,(double)photonCounter));
      }
    }
    float screen_ratio  = (float)sizeofImage / resolutionMax;
    bool  pNeedsDrawing = (photonCounter == 1 || checkerOdd(rowCounter) || (!checkerOdd(rowCounter) && checkerOdd(colCounter)));
    x = colCounter * screen_ratio;
    y = rowCounter * screen_ratio;
    colCounter++;

    if (pNeedsDrawing){
      iterations++;
      rgb = colorPixelCalculation(x,y);

      //render pixel by pixel
      glColor3d(rgb[0],rgb[1],rgb[2]);
      glPointSize((float)screen_ratio);
      glBegin(GL_POINTS);
      glVertex2d((double)x,(double)HIEGHT-y);
      glEnd();

    }
  }
  if (rowCounter == sizeofImage-1) {empty = false;}
}

void resetRender(){ //Reset Rendering Variables
  rowCounter=0; colCounter=0; photonCounter=1; resolutionMax=2;
  empty=true;
  if (photonMapping) emitPhotons();
}

void drawPhoton(const CoordinateVector3 &rgb, const CoordinateVector3 &p){           //Photon Visualization
  
}

/* 
Mouse and Keyboard
*/

int prevMouseX = -9999, prevMouseY = -9999, sphereIndex = -1;
float s = 130.0;

void onKeyPress(unsigned char key,int, int) {
  switch(key) {
    case 'r'  : photonMapping = false; break;
    case 'p'  : photonMapping = true; break;
    default     : return;
  }
  resetRender();
  printf("No. %d key pressed\n",key);
}

void onClick(int button,int action, int x, int y) {
  if(button != 0) return;
  if(action == 0) {
    mouseDragging = true;

    sphereIndex = numOfSurfcae;

    mouseX = x;
    mouseY = y;

    float mousecoord[] = {
       (mouseX - sizeofImage/2)/s,
      -(mouseY - sizeofImage/2)/s
    };

    for(int i=0; i<numOfSurfcae; i++) {
      Surfaces *ob = surfaces[i];
      if (ob->getSurface() != SPHERE) { continue; }
      CoordinateVector3 mouse2screen(
          mousecoord[0],
          mousecoord[1],
          ob->coords[2]
          );
      CoordinateVector3 center(ob->coords[0], ob->coords[1], ob->coords[2]);
      if (distance(mouse2screen, center) < ob->coords[3]) { sphereIndex = i; }
    }
  }
  else {
    prevMouseX = -9999;
    prevMouseY = -9999;
    mouseDragging = false;
  }
}

void onDrag(int x,int y) {
  //current mouse coords
  mouseX = x;
  mouseY = y;

  if(mouseDragging) {
    if (prevMouseX > -9999 && sphereIndex > -1){
      if (sphereIndex < numOfSurfcae){ //Drag Sphere
        surfaces[sphereIndex]->coords[0] += (mouseX - prevMouseX)/s;
        surfaces[sphereIndex]->coords[1] -= (mouseY - prevMouseY)/s;
      }else{ //Drag Light
        Light = CoordinateVector3(
            constrain(Light[0] + (mouseX - prevMouseX)/s, -1.4, 1.4),
            constrain(Light[1] - (mouseY - prevMouseY)/s, -0.4, 1.2),
            Light[2] );
      }
      resetRender();
    }
    prevMouseX = mouseX;
    prevMouseY = mouseY;
  }
}

void onTimer(int val) {
  glutPostRedisplay();
  glutTimerFunc(0, onTimer, val);
}


void surfaceInitialization() {
  static const float white[3] = {1.0,1.0,1.0};
  static const float red[3]   = {1.0,0.0,0.0};
  static const float green[3] = {0.0,1.0,0.0};
  static const float blue[3]  = {0.0,0.0,1.0};

  //{center(x,y,z), radius}
  float v_sphere[][4] = {
    { 1.0,  0.0, 4.0, 0.25},
    {-0.5,  0.0, 4.5, 0.75},
    { 0.2, -0.9, 4.0, 0.3},
  };

  //{(axis_id), (distance_from_origin)}
  //axis_id = 0:X, 1:Y, 2:Z
  float v_plane[][2]  = {
    {0,  1.5},
    {1, -1.5},
    {0, -1.5},
    {1,  1.5},
    {2,  5.0}
  };

  surfaces.resize(0);

  Surfaces *ob;

  //create spheres
  for(int i=0; i<3; i++) {
    ob = new Surfaces(SPHERE,numOfSurfcae++,v_sphere[i]);
    surfaces.push_back(ob);
  }

  //create planes
  for(int i=0; i<5; i++) {
    ob = new Surfaces(PLANE,numOfSurfcae++,v_plane[i]);
    surfaces.push_back(ob);
  }

  //set optical properties
  surfaces[0]->setSpecialEffect(REFLECTION);
  surfaces[1]->setSpecialEffect(REFLECTION);
  surfaces[2]->setSpecialEffect(REFRACTION);
  surfaces[2]->setRefractive(2.5f);
  surfaces[3]->setColor(red);
  surfaces[5]->setColor(green);
  surfaces[7]->setColor(white); 
  surfaces[4]->setColor(white);
  surfaces[6]->setColor(white);
}

void deleteSurfaces() {
  for(int i=0; i<numOfSurfcae; i++) { delete surfaces[i]; }
  surfaces.clear();
}


