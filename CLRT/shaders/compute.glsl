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

// Albedo is Color guys
struct Material {
  vec3 Albedo;
  float Roughness;
  float Metallic;
};

struct DirectLight {
  vec3 Position;
  float Radius;
  vec3 Albedo;
  float Power;
  float Reach;
};

struct SurfacePoint {
  vec3 Position;
  vec3 Normal;
  Material SurfaceMaterial;
};

struct Ray {
  vec3 origin;
  vec3 direction;
};

struct HitPayload {
  vec3 HitPoint;
  vec3 WorldNormal;
  float HitDistance;
  bool FrontFace;
};

HitPayload gRec;

void SetFaceNormal(const Ray ray, const vec3 outwardNormal) {
  gRec.FrontFace = dot(ray.direction, outwardNormal) < 0;
  gRec.WorldNormal = gRec.FrontFace ? outwardNormal : -outwardNormal;
}

struct Sphere {
  vec3 Position;
	float Radius;

  int MaterialIndex;
};

struct HittableList {
  Sphere spheres[MAX_NUMBER_OF_OBJECTS];
};

layout(location = 1) uniform vec3 cam_origin;
layout(location = 2) uniform mat4 cam_view;
layout(location = 3) uniform mat4 cam_proj;
layout(location = 4) uniform vec3 cam_fdir;
layout(location = 5) uniform vec3 cam_up;
layout(location = 6) uniform mat4 rotationMatrix;
layout(location = 7) uniform float deltaTime;
layout(location = 8) uniform float firstTime;

int objectCount;

Sphere firstSphere = { vec3(0.0, 0.0, -1.0), 1.0, 0 };
Sphere secondSphere = { vec3(0.0, -101.0, -1.0), 100.0, 1 };

Material firstMaterial = { vec3(1.0, 0.0, 1.0), 0.0, 0.0 };
Material secondMaterial = { vec3(0.2, 0.3, 1.0), 1, 0.0 };

DirectLight firstLight = { vec3(-1.0, 1.0, 1.0), 1.0, vec3(0, 0, 0), 1.0, 100.0 };

uint wang_hash(uint seed)
{
  seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
  seed *= uint(9);
  seed = seed ^ (seed >> 4);
  seed *= uint(0x27d4eb2d);
  seed = seed ^ (seed >> 15);
  return seed;
}

float RandomFloat01(uint state)
{
  return float(wang_hash(state)) / 4294967296.0;
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

// Converts the RGB Value from 255 to a value between 0 and 1, or the other way around
vec3 EightBitIntToPercentColor(vec3 rgbColor, bool convertTo) {
  if (convertTo) return vec3(rgbColor.x / 255, rgbColor.y / 255, rgbColor.z / 255);

  return vec3(rgbColor.x * 255, rgbColor.y * 255, rgbColor.z * 255);
}

vec3 at(float t, const Ray ray) { return ray.origin + t * ray.direction; }

bool SphereHit(const Ray ray, float t_min, float t_max, HitPayload rec, const Sphere sphere) {
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

  return true;
}

bool ListHit(const Ray ray, float t_min, float t_max, HitPayload rec, HittableList hitList) {
  HitPayload tempRec;
  bool hitAnything = false;
  float closestSoFar = t_max;

  for (int i = 0; i < objectCount; i++) {
    if (SphereHit(ray, t_min, closestSoFar, tempRec, hitList.spheres[i])) {
      hitAnything = true;
      closestSoFar = gRec.HitDistance;
      gRec = tempRec;
    }
  }

  return hitAnything;
}

vec3 RayColor(float rngState, Ray ray, const HittableList world, int depth) {
  HitPayload rec;
  if (ListHit(ray, 0, INFINITY, rec, world))
    return 0.5 * (gRec.WorldNormal + vec3(1));

  float t = 0.5 * (ray.direction.y + 1.0);
  return (1.0 - t) * vec3(1) + t * vec3(0.5, 0.7, 1.0);
}

void main() {
  vec4 pixel = vec4(0.0);
  ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

  ivec2 dims = imageSize(screen);
  float x = (float(pixel_coords.x * 2 - dims.x) / dims.x);
  float y = (float(pixel_coords.y * 2 - dims.y) / dims.x);
  
  // This is for Antialiasing!!!
  const int samples_per_pixel = 50;

  HittableList world;
  world.spheres[0] = firstSphere;
  world.spheres[1] = secondSphere;

  float fov = 90.0;
  int max_depth = 50;

  for (int i = 0; i < world.spheres.length(); i++)
    if (dot(world.spheres[i].Position, world.spheres[i].Position) + world.spheres[i].Radius > 0.0)
      objectCount = i+1;

  // Defines the required Variables to shoot the Rays into random Directions (I made this extra thing to make it more "Random" than it was before)
  // If I didn't add the last 2 things, I would have had the same values everytime!
  uint rngState = uint(uint(pixel_coords.x) * uint(1973) + uint(pixel_coords.y) * uint(9277) + uint(pixel) * uint(26699)) | uint(1);
  float fRngState = rngState;
  fRngState += firstTime * fRngState;

  Ray ray = { cam_origin, (normalize(vec4(vec2(x, y), -1.0, 0.0)) * rotationMatrix).xyz };

  world.spheres[0] = firstSphere;
  world.spheres[1] = secondSphere;

  vec3 rayColor = vec3(0);
  rayColor += RayColor(uint(fRngState), ray, world, max_depth);

  pixel = vec4(rayColor, 1.0);

  imageStore(screen, pixel_coords, pixel);
}