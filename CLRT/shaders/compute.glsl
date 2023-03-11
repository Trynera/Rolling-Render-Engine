#version 460 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image2D screen;
uniform int accumulatedPasses;

#define MAX_NUMBER_OF_OBJECTS 128
#define MAX_LIGHT_COUNT 4
#define FLOAT_MAX 3.4028235e+38
#define FLOAT_MIN -3.4028235e+38
#define UINT_MAX 4.29497e+09
#define EPSILON 0.0001
#define PI 3.1415926535897932385
#define TWOPI 2.0 * PI
#define INFINITY 1. / 0.
#define BOUNCES 5
#define INV_UINT_MAX (1.0/4294967296.0);
#define SAMPLES_PER_PIXEL 8

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct HitPayload {
  vec3 HitPoint;
  vec3 WorldNormal;
  float HitDistance;
  float U;
  float V;
  bool FrontFace;

  int MaterialIndex;
};

struct CheckTex {
  vec3 Odd;
  vec3 Even;
};

struct Texture {
  CheckTex CTex;
  vec3 SCol;
  int Type;
};

struct Material {
  Texture Tex;
  float Fuzz;
  float IR;
  Texture Emit;
  int Type;
};

HitPayload gRec;

void SetFaceNormal(const Ray ray, const vec3 outwardNormal) {
  gRec.FrontFace = dot(ray.direction, outwardNormal) < 0;
  gRec.WorldNormal = gRec.FrontFace ? outwardNormal : -outwardNormal;
}

struct Sphere {
  vec3 Position;
	float Radius;
  bool Movable;

  int MaterialIndex;
};

struct XYRect {
  float x0;
  float x1;
  float y0;
  float y1;
  float k;
  int MaterialIndex;
};

struct XZRect {
  float x0;
  float x1;
  float z0;
  float z1;
  float k;
  int MaterialIndex;
};

struct YZRect {
  float y0;
  float y1;
  float z0;
  float z1;
  float k;
  int MaterialIndex;
};

struct BoxHitList {
  XYRect xyrects[MAX_NUMBER_OF_OBJECTS];
  XZRect xzrects[MAX_NUMBER_OF_OBJECTS];
  YZRect yzrects[MAX_NUMBER_OF_OBJECTS];
};

struct Box {
  vec3 BoxMin;
  vec3 BoxMax;
  BoxHitList sides;

  int MaterialIndex;
};

struct Aabb {
  vec3 Min;
  vec3 Max;
};

struct BVH {
  HitPayload left;
  HitPayload right;
  Aabb box;
};

struct HittableList {
  Sphere spheres[MAX_NUMBER_OF_OBJECTS];
  XYRect xyrects[MAX_NUMBER_OF_OBJECTS];
  XZRect xzrects[MAX_NUMBER_OF_OBJECTS];
  YZRect yzrects[MAX_NUMBER_OF_OBJECTS];
  Box boxes[MAX_NUMBER_OF_OBJECTS];
  Material mats[MAX_NUMBER_OF_OBJECTS];
};

struct SurfacePoint {
  vec3 Position;
  vec3 Normal;
  Material SurfaceMaterial;
};

layout(location = 1) uniform vec3 cam_origin;
layout(location = 2) uniform mat4 cam_view;
layout(location = 3) uniform mat4 cam_proj;
layout(location = 4) uniform vec3 cam_fdir;
layout(location = 5) uniform vec3 cam_up;
layout(location = 6) uniform mat4 rotationMatrix;
layout(location = 7) uniform float deltaTime;
layout(location = 8) uniform float firstTime;

int matCount;
int objectCount;

/*
  There are 4 types:
  0: Lamb
  1: Metal
  2: Diel
  3: Emissive
*/

YZRect greenWall = { 0, 555, 0, 555, 555, 2 };
YZRect redWall = { 0, 555, 0, 555, 0, 0 };
XZRect lightSrc = { 213, 343, 227, 332, 554, 3 };
XZRect whiteFloor = { 0, 555, 0, 555, 0, 1 };
XZRect whiteRoof = { 0, 555, 0, 555, 555, 1 };
XYRect whiteWall = { 0, 555, 0, 555, 555, 1 };

Box bigBoyBox;
Box smolBoyBox;

CheckTex checkVal1 = { vec3(0), vec3(1) };
CheckTex checkVal2 = { vec3(0), vec3(0) };

Texture tex1 = { checkVal2, vec3(.65, .05, .05), 0 };
Texture tex2 = { checkVal2, vec3(.73, .73, .73), 0 };
Texture tex3 = { checkVal2, vec3(.12, .45, .15), 0 };
Texture tex4 = { checkVal2, vec3(15), 1 };

Texture NoEmit = { checkVal2, vec3(0), 0 };

Material red = { tex1, 0, 0, NoEmit, 0 };
Material white = { tex2, 0, 0, NoEmit, 0 };
Material green = { tex3, 0, 0, NoEmit, 0 };
Material light = { tex4, 0, 0, tex4, 3 };

uint wang_hash(uint seed)
{
  seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
  seed *= uint(9);
  seed = seed ^ (seed >> 4);
  seed *= uint(0x27d4eb2d);
  seed = seed ^ (seed >> 15);
  return seed;
}

vec3 Emitted(const float u, const float v, const vec3 p, inout Material lightSrc) {
  
  lightSrc.Emit.SCol = vec3(u, v, p);
  return lightSrc.Emit.SCol;
}

float RandomFloat01(uint state)
{
  return float(wang_hash(state)) / 4294967296.0;
}

uint RandXor(uint state) {
  state ^= (state << 13);
  state ^= (state >> 17);
  state ^= (state << 5);
  return wang_hash(state);
}

vec3 RandomUnitVector(uint state)
{
  float z = (RandomFloat01(state) * 2.0 - 1.0);
  float a = (RandomFloat01(state) * TWOPI);
  float r = (1.0 - z * z);
  float nX = r * cos(a);
  float nY = r * sin(a);
  return vec3(nX, nY, z);
}

uint BaseHash(uvec2 p) {
  p = 1103515245U * ((p >> 1U) ^ (p.yx));
  uint h32 = 1103515245U * ((p.x) ^ (p.y>>3U));
  return h32 ^ (h32 >> 16);
}

void CreateBox(const vec3 p0, const vec3 p1, const Material mat, inout Box box) {
  box.BoxMin = p0;
  box.BoxMax = p1;

  box.sides.xyrects[0].x0 = p0.x;
  box.sides.xyrects[0].x1 = p1.x;
  box.sides.xyrects[0].y0 = p0.y;
  box.sides.xyrects[0].y1 = p1.y;
  box.sides.xyrects[0].k = p1.z;

  box.sides.xyrects[1].x0 = p0.x;
  box.sides.xyrects[1].x1 = p1.x;
  box.sides.xyrects[1].y0 = p0.y;
  box.sides.xyrects[1].x1 = p1.y;
  box.sides.xyrects[1].k = p0.z;

  box.sides.xzrects[0].x0 = p0.x;
  box.sides.xzrects[0].x1 = p1.x;
  box.sides.xzrects[0].z0 = p0.z;
  box.sides.xzrects[0].z1 = p1.z;
  box.sides.xzrects[0].k = p1.y;

  box.sides.xzrects[1].x0 = p0.x;
  box.sides.xzrects[1].x1 = p1.x;
  box.sides.xzrects[1].z0 = p0.z;
  box.sides.xzrects[1].z1 = p1.z;
  box.sides.xzrects[1].k = p0.y;

  box.sides.yzrects[0].y0 = p0.y;
  box.sides.yzrects[0].y1 = p1.y;
  box.sides.yzrects[0].z0 = p0.z;
  box.sides.yzrects[0].z1 = p1.z;
  box.sides.yzrects[0].k = p1.x;

  box.sides.yzrects[1].y0 = p0.y;
  box.sides.yzrects[1].y1 = p1.y;
  box.sides.yzrects[1].z0 = p0.z;
  box.sides.yzrects[1].z1 = p1.z;
  box.sides.yzrects[1].k = p0.x;
}

vec3 Hash(float seed) {
  uint n = BaseHash(floatBitsToUint(vec2(seed += 0.1, seed += 0.1)));
  uvec3 rz = uvec3(n, n * 16807U, n * 48271U);
  return vec3(rz & uvec3(0x7fffffffU)) / float(0x7fffffff);
}

vec3 RandomInUnitSphere(float seed) {
  vec3 h = Hash(seed) * vec3(2.0, 6.28318530718, 1.0) - vec3(1.0, 0.0, 0.0);
  float phi = h.y;
  float r = pow(h.z, 1.0/3.0);
  return r * vec3(sqrt(1.0 - h.x * h.x) * vec2(sin(phi), cos(phi)), h.x);
}

// Converts the RGB Value from 255 to a value between 0 and 1, or the other way around
vec3 EightBitIntToPercentColor(vec3 rgbColor, bool convertTo) {
  if (convertTo) return vec3(rgbColor.x / 255, rgbColor.y / 255, rgbColor.z / 255);

  return vec3(rgbColor.x * 255, rgbColor.y * 255, rgbColor.z * 255);
}

float fabsf(float x) {
  if (x < 0)
    return -x;
  return x;
}

float reflectance(const float cosine, const float refIdx) {
  float r0 = (1-refIdx) / (1+refIdx);
  r0 = r0*r0;
  return r0 + (1-r0)*pow((1-cosine),5);
}

bool NearZero(const vec3 scaDir) {
  const float s = 1e-8;
  return (fabsf(scaDir[0]) < s) && (fabsf(scaDir[1]) < s) && (fabsf(scaDir[2]) < s);
}

float fmin(float x, float y) {
  if (x > y)
    return y;
  return x;
}

float fmax(float x, float y) {
  if (x < y)
    return y;
  return x;
}

bool modified_refract(const in vec3 v, const in vec3 n, const in float ni_over_nt, 
                      out vec3 refracted) {
    float dt = dot(v, n);
    float discriminant = 1. - ni_over_nt*ni_over_nt*(1.-dt*dt);
    if (discriminant > 0.) {
        refracted = ni_over_nt*(v - n*dt) - n*sqrt(discriminant);
        return true;
    } else { 
        return false;
    }
}

vec3 CheckerValue(const float u, const float v, const vec3 p, CheckTex CTex) {
  float sines = sin(10*p.x)*sin(10*p.y)*sin(10*p.z);
  if (sines < 0)
    return CTex.Odd;
  return CTex.Even;
}

bool ScatterLamb(const Ray ray, const HitPayload rec, inout vec3 attenuation, inout Ray scattered, float rngState, const Material mat) {
  vec3 scatterDirection = gRec.WorldNormal + RandomInUnitSphere(rngState);

  if (NearZero(scatterDirection)) scatterDirection = gRec.WorldNormal;
  
  scattered.origin = gRec.HitPoint;
  scattered.direction = scatterDirection;
  attenuation = mat.Tex.SCol;
  if (mat.Tex.Type == 1) {
    attenuation = CheckerValue(gRec.U, gRec.V, gRec.HitPoint, mat.Tex.CTex);
  }

  return true;
}

bool ScatterMetal(const Ray ray, const HitPayload rec, inout vec3 attenuation, inout Ray scattered, float rngState, const Material mat) {
  vec3 reflected = reflect(ray.direction, gRec.WorldNormal);
  scattered.origin = gRec.HitPoint;
  scattered.direction = reflected + mat.Fuzz * RandomInUnitSphere(rngState);
  attenuation = mat.Tex.SCol;
  if (mat.Tex.Type == 1) {
    attenuation = CheckerValue(gRec.U, gRec.V, gRec.HitPoint, mat.Tex.CTex);
  }

  return (dot(scattered.direction, gRec.WorldNormal) > 0);
}

bool ScatterDiel(const Ray ray, const HitPayload rec, inout vec3 attenuation, inout Ray scattered, float rngState, const Material mat) {
  attenuation = vec3(1);
  float refractionRatio = gRec.FrontFace ? (1.0 / mat.IR) : mat.IR;

  float cosTheta = fmin(dot(-ray.direction, gRec.WorldNormal), 1.0);
  float sinTheta = sqrt(1.0 - cosTheta*cosTheta);

  bool cannotRefract = refractionRatio * sinTheta > 1.0;
  vec3 direction = refract(ray.direction, gRec.WorldNormal, refractionRatio);

  if (cannotRefract)
    direction = reflect(ray.direction, gRec.WorldNormal);
  
  scattered.origin = gRec.HitPoint;
  scattered.direction = direction;
  return true;
}

bool Scatter(const Ray ray, const HitPayload rec, inout vec3 attenuation, inout Ray scattered, const Material mat, float rngState, int Type) {
  switch (Type) {
  case 0: {
    vec3 tempAtten;
    return ScatterLamb(ray, rec, tempAtten, scattered, rngState, mat);
  }
  case 1: return ScatterMetal(ray, rec, attenuation, scattered, rngState, mat);
  case 2: return ScatterDiel(ray, rec, attenuation, scattered, rngState, mat);
  default: return false;
  }
}

vec3 at(float t, const Ray ray) { return ray.origin + t * ray.direction; }

void GetSphereUV(const vec3 p, inout float U, inout float V, const Sphere sphere) {
  float theta = acos(-p.y);
  float phi = atan(p.x, -p.z) + PI;

  U = phi / (2*PI);
  V = theta / PI;
}

bool Hit(const Ray ray, float t_min, float t_max, HitPayload rec, const Sphere sphere) {
  vec3 oc = ray.origin - sphere.Position;
  float a = dot(length(ray.direction), length(ray.direction));
  float half_b = dot(oc, ray.direction);
  float c = dot(length(oc), length(oc)) - sphere.Radius * sphere.Radius;

  float discriminant = half_b * half_b - a * c;
  if (discriminant < 0) return false;
  float sqrtd = sqrt(discriminant);

  float root = (-half_b - sqrtd) / a;
  if (root < t_min || t_max < root) {
    root = (-half_b + sqrtd) / a;
    if (root < t_min || t_max < root)
      return false;
  }

  gRec.HitDistance = root;
  gRec.HitPoint = at(gRec.HitDistance, ray);
  vec3 outwardNormal = (gRec.HitPoint - sphere.Position) / sphere.Radius;
  SetFaceNormal(ray, outwardNormal);
  GetSphereUV(outwardNormal, gRec.U, gRec.V, sphere);

  return true;
}

bool Hit(const Ray ray, float t_min, float t_max, HitPayload rec, const XYRect rect) {
  float t = (rect.k-ray.origin.z / ray.direction.z);
  
  if (t < t_min || t > t_max)
    return false;

  float x = ray.origin.x + t*ray.direction.x;
  float y = ray.origin.y + t*ray.direction.y;

  if (x < rect.x0 || x > rect.x1 || y < rect.y0 || y > rect.y1)
    return false;
  
  gRec.U = (x-rect.x0)/(rect.x1-rect.x0);
  gRec.V = (y-rect.y0)/(rect.y1-rect.y0);
  gRec.HitDistance = t;
  vec3 outwardNormal = vec3(0, 0, 1);
  SetFaceNormal(ray, outwardNormal);
  gRec.MaterialIndex = rect.MaterialIndex;
  gRec.HitPoint = at(t, ray);
  return true;
}

bool Hit(const Ray ray, float t_min, float t_max, HitPayload rec, const XZRect rect) {
  float t = (rect.k-ray.origin.y) / ray.direction.y;
  if (t < t_min || t > t_max)
    return false;
  
  float x = ray.origin.x + t*ray.direction.x;
  float z = ray.origin.z + t*ray.direction.z;

  if (x < rect.x0 || x > rect.x1 || z < rect.z0 || z > rect.z1)
    return false;
  
  gRec.U = (x-rect.x0)/(rect.x1-rect.x0);
  gRec.V = (z-rect.z0)/(rect.z1-rect.z0);
  gRec.HitDistance = t;
  vec3 outwardNormal = vec3(0, 1, 0);
  SetFaceNormal(ray, outwardNormal);
  gRec.MaterialIndex = rect.MaterialIndex;
  gRec.HitPoint = at(t, ray);
  return true;
}

bool Hit(const Ray ray, float t_min, float t_max, HitPayload rec, const YZRect rect) {
  float t = (rect.k-ray.origin.x) / ray.direction.x;
  if (t < t_min || t > t_max)
    return false;
  
  float y = ray.origin.y + t*ray.direction.y;
  float z = ray.origin.z + t*ray.direction.z;
  
  if (y < rect.y0 || y > rect.y1 || z < rect.z0 || z > rect.z1)
    return false;
  
  gRec.U = (y-rect.y0)/(rect.y1-rect.y0);
  gRec.V = (z-rect.z0)/(rect.z1-rect.z0);
  gRec.HitDistance = t;
  vec3 outwardNormal = vec3(1, 0, 0);
  SetFaceNormal(ray, outwardNormal);
  gRec.MaterialIndex = rect.MaterialIndex;
  gRec.HitPoint = at(t, ray);
  return true;
}

void Swap(inout float x, inout float y) {
  const float z = x; x=y; y=z;
}

bool Hit(const Ray ray, float t_min, float t_max, HitPayload rec, HittableList hitList, const Aabb aabb) {
  for (int i = 0; i < 3; i++) {
    float invD = 1.0 / ray.direction[i];
    float t0 = (aabb.Min[i] - ray.origin[i]) * invD;
    float t1 = (aabb.Max[i] - ray.origin[i]) * invD;
    if (invD < 0.0)
      Swap(t0, t1);
    t_min = t0 > t_min ? t0 : t_min;
    t_max = t1 < t_max ? t1 : t_max;
    if (t_max <= t_min)
      return false;
  }
  return true;
}

bool BVHHit(const Ray ray, float t_min, float t_max, const BVH mainBVH, HittableList hitList, Sphere sphere) {
  if (!Hit(ray, t_min, t_max, gRec, hitList, mainBVH.box)) return false;

  HitPayload tempGRec = gRec;
  gRec = mainBVH.left;
  bool hitLeft = Hit(ray, t_min, t_max, gRec, sphere);
  gRec = mainBVH.right;
  bool hitRight = Hit(ray, t_min, hitLeft ? gRec.HitDistance : t_max, gRec, sphere);
  gRec = tempGRec;

  return hitLeft || hitRight;
}

Aabb SurroundingBox(Aabb box0, Aabb box1) {
  vec3 Small = vec3(fmin(box0.Min.x, box1.Min.x),
             fmin(box0.Min.y, box1.Min.y),
             fmin(box0.Min.z, box1.Min.z));

  vec3 Big = vec3(fmax(box0.Max.x, box1.Max.x),
           fmax(box0.Max.y, box1.Max.y),
           fmax(box0.Max.z, box1.Max.z));
  
  return Aabb(Small, Big);
}

bool BoundingBox(const float time0, const float time1, inout Aabb outBox, const XYRect mRect) {
  outBox.Min = vec3(mRect.x0, mRect.y0, mRect.k-0.0001);
  outBox.Max = vec3(mRect.x1, mRect.y1, mRect.k-0.0001);
  return true;
}

bool BoundingBox(const float time0, const float time1, inout Aabb outBox, const XZRect mRect) {
  outBox.Min = vec3(mRect.x0, mRect.z0, mRect.k-0.0001);
  outBox.Max = vec3(mRect.x1, mRect.z1, mRect.k-0.0001);
  return true;
}

bool BoundingBox(const float time0, const float time1, inout Aabb outBox, const YZRect mRect) {
  outBox.Min = vec3(mRect.y0, mRect.z0, mRect.k-0.0001);
  outBox.Max = vec3(mRect.y1, mRect.z1, mRect.k-0.0001);
  return true;
}

bool SphereBoundingBox(const float time0, const float time1, inout Aabb outBox, const Sphere sphere) {
  if (!sphere.Movable) {
    outBox = Aabb(
      sphere.Position - vec3(sphere.Radius),
      sphere.Position + vec3(sphere.Radius));
    return true;
  }

  Aabb box0 = Aabb(
    vec3(time0) - vec3(sphere.Radius),
    vec3(time0) + vec3(sphere.Radius));
  Aabb box1 = Aabb(
    vec3(time1) - vec3(sphere.Radius),
    vec3(time1) + vec3(sphere.Radius));
  outBox = SurroundingBox(box0, box1);
  return true;
}

bool BoxCompare(const Sphere a, const Sphere b, int axis) {
  Aabb boxA;
  Aabb boxB;

  if (SphereBoundingBox(0, 0, boxA, a) || SphereBoundingBox(0, 0, boxB, b))
    return false;

  return boxA.Min[axis] < boxB.Min[axis];
}

bool BoxXCompare(const Sphere a, const Sphere b) { return BoxCompare(a, b, 0); }
bool BoxYCompare(const Sphere a, const Sphere b) { return BoxCompare(a, b, 1); }
bool BoxZCompare(const Sphere a, const Sphere b) { return BoxCompare(a, b, 2); }

bool ListBoundingBox(const float time0, const float time1, inout Aabb outBox, const HittableList world) {
  if (world.spheres[0].Radius == 0) return false;

  Aabb tempBox;
  bool firstBox = true;

  for (int i = 0; i < objectCount; i++) {
    if (!SphereBoundingBox(time0, time1, tempBox, world.spheres[i])) return false;
    outBox = firstBox ? tempBox : SurroundingBox(outBox, tempBox);
    firstBox = false;
  }

  return true;
}

bool BVHBoundingBox(const float time0, const float time1, inout Aabb outBox, const BVH mainBVH) {
  outBox = mainBVH.box;
  return true;
}

bool ListHit(const Ray ray, float t_min, float t_max, HitPayload rec, HittableList hitList) {
  HitPayload tempRec;
  bool hitAnything = false;
  float closestSoFar = t_max;

  for (int i = 0; i < objectCount; i++) {
    if (Hit(ray, t_min, closestSoFar, tempRec, hitList.spheres[i])) {
      hitAnything = true;
      closestSoFar = gRec.HitDistance;
      gRec = tempRec;
      gRec.MaterialIndex = hitList.spheres[i].MaterialIndex;
      break;
    }
    if (Hit(ray, t_min, closestSoFar, tempRec, hitList.xyrects[i])) {
      hitAnything = true;
      closestSoFar = gRec.HitDistance;
      gRec = tempRec;
      gRec.MaterialIndex = hitList.xyrects[i].MaterialIndex;
      break;
    }
    if (Hit(ray, t_min, closestSoFar, tempRec, hitList.xzrects[i])) {
      hitAnything = true;
      closestSoFar = gRec.HitDistance;
      gRec = tempRec;
      gRec.MaterialIndex = hitList.xzrects[i].MaterialIndex;
      break;
    }
    if (Hit(ray, t_min, closestSoFar, tempRec, hitList.yzrects[i])) {
      hitAnything = true;
      closestSoFar = gRec.HitDistance;
      gRec = tempRec;
      gRec.MaterialIndex = hitList.yzrects[i].MaterialIndex;
      break;
    }
  }

  return hitAnything;
}

float rand(vec2 co) {
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

vec3 randVec(vec2 co) {
  return vec3(rand(co), rand(co), rand(co));
}

vec3 RayColor(float rngState, Ray ray, const vec3 background, HittableList world) {
  HitPayload rec;
  Ray scattered;
  vec3 attenuation = vec3(1);
  vec3 emitted = Emitted(gRec.U, gRec.V, gRec.HitPoint, world.mats[gRec.MaterialIndex]);
  vec3 color = vec3(1);

  float multiplier = 1;
  int maxBounces = 5;
  for (int i = 0; i < maxBounces; i++) {
    emitted = Emitted(gRec.U, gRec.V, gRec.HitPoint, world.mats[gRec.MaterialIndex]);
    if (!ListHit(ray, 0.001, INFINITY, rec, world))
        return background;

    if (!Scatter(ray, gRec, attenuation, scattered, world.mats[gRec.MaterialIndex], rngState, world.mats[gRec.MaterialIndex].Type))
      return emitted + attenuation * color;

    ray = scattered;
    attenuation *= multiplier;
    color = emitted + attenuation * color;
    emitted = emitted + attenuation * emitted;
    continue;
  }

  return vec3(0);
}

void main() {
  vec4 pixel = vec4(0.0);
  ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

  ivec2 dims = imageSize(screen);
  
  float x = (float((pixel_coords.x + 1.0) * 2 - dims.x) / dims.x);
  float y = (float((pixel_coords.y + 0.1) * 2 - dims.y) / dims.x);

  HittableList world;

  world.mats[0] = red;
  world.mats[1] = white;
  world.mats[2] = green;
  world.mats[3] = light;

  CreateBox(vec3(130, 0, 65), vec3(295, 165, 230), white, smolBoyBox);
  CreateBox(vec3(265, 0, 295), vec3(430, 330, 460), white, bigBoyBox);

  world.yzrects[0] = greenWall;
  world.yzrects[1] = redWall;
  world.yzrects[2] = smolBoyBox.sides.yzrects[0];
  world.yzrects[3] = smolBoyBox.sides.yzrects[1];
  world.yzrects[5] = bigBoyBox.sides.yzrects[0];
  world.yzrects[6] = bigBoyBox.sides.yzrects[1];

  world.xzrects[0] = lightSrc;
  world.xzrects[1] = whiteFloor;
  world.xzrects[2] = whiteRoof;
  world.xzrects[3] = smolBoyBox.sides.xzrects[0];
  world.xzrects[4] = smolBoyBox.sides.xzrects[1];
  world.xzrects[5] = bigBoyBox.sides.xzrects[0];
  world.xzrects[6] = bigBoyBox.sides.xzrects[1];

  world.xyrects[0] = whiteWall;
  world.xyrects[1] = smolBoyBox.sides.xyrects[0];
  world.xyrects[2] = smolBoyBox.sides.xyrects[1];
  world.xyrects[3] = bigBoyBox.sides.xyrects[0];
  world.xyrects[4] = bigBoyBox.sides.xyrects[1];

  vec3 background = vec3(0);
  float fov = 40.0;

  // Defines the required Variables to shoot the Rays into random Directions (I made this extra thing to make it more "Random" than it was before)
  // If I didn't add the last 2 things, I would have had the same values everytime!
  uint rngState = uint(uint(pixel_coords.x) * uint(1973) + uint(pixel_coords.y) * uint(9277) + uint(pixel) * uint(26699)) | uint(1);
  float fRngState = rngState;
  fRngState += firstTime * fRngState;

  Ray ray = { cam_origin, (normalize(vec4(vec2(x, y), -1.0, 0.0)) * rotationMatrix).xyz };

  objectCount = 18;

  vec3 rayColor = RayColor(fRngState, ray, background, world);
  rayColor = (sqrt(rayColor));
  
  pixel = vec4(rayColor, 1.0);

  imageStore(screen, pixel_coords, pixel);
}