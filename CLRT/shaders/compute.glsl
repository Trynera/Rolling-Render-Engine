#version 460 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image2D screen;

#define MAX_NUMBER_OF_OBJECTS 128
#define MAX_LIGHT_COUNT 4
#define FLOAT_MAX 3.4028235e+38
#define FLOAT_MIN -3.4028235e+38
#define UINT_MAX 4.29497e+09
#define EPSILON 0.0001
#define PI 3.14159265
#define TWOPI 2.0*PI

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
  float HitDistance;
  vec3 WorldPosition;
  vec3 WorldNormal;

  int ObjectIndex;
};

struct Sphere {
  vec3 Position;
	float Radius;

  int MaterialIndex;
};

struct Scene {
  Sphere spheres[MAX_NUMBER_OF_OBJECTS];
  Material materials[MAX_NUMBER_OF_OBJECTS];
  DirectLight lightSources[MAX_LIGHT_COUNT];
};

layout(location = 1) uniform vec3 cam_origin;
layout(location = 2) uniform mat4 cam_view;
layout(location = 3) uniform mat4 cam_proj;
layout(location = 4) uniform vec3 cam_fdir;
layout(location = 5) uniform vec3 cam_up;
layout(location = 6) uniform mat4 rotationMatrix;
layout(location = 7) uniform float deltaTime;

uint frameIndex = 1;

Scene currentScene;
Sphere firstSphere = { vec3(0.0, 0.0, -1.0), 1.0, 0 };
Sphere secondSphere = { vec3(0.0, -101.0, -1.0), 100.0, 1 };

Material firstMaterial = { vec3(1.0, 0.0, 1.0), 0.0, 0.0 };
Material secondMaterial = { vec3(0.2, 0.3, 1.0), 0.04, 0.0 };

uint wang_hash(inout uint seed)
{
    seed = uint(seed ^ uint(61)) ^ uint(seed >> uint(16));
    seed *= uint(9);
    seed = seed ^ (seed >> 4);
    seed *= uint(0x27d4eb2d);
    seed = seed ^ (seed >> 15);
    return seed;
}

float RandomFloat01(inout uint state)
{
    return float(wang_hash(state)) / 4294967296.0;
}

vec3 RandomUnitVector(inout uint state)
{
    float z = RandomFloat01(state) * 2.0 - 1.0;
    float a = RandomFloat01(state) * TWOPI;
    float r = sqrt(1.0 - z * z);
    float nX = r * cos(a);
    float nY = r * sin(a);
    return vec3(nX, nY, z);
}

HitPayload Miss(const Ray ray) {
  HitPayload payload;
  payload.HitDistance = -1.0;
  return payload;
}

HitPayload ClosestHit(const Ray ray, float hitDistance, int objectIndex) {
  HitPayload payload;
  payload.HitDistance = hitDistance;
  payload.ObjectIndex = objectIndex;

  const Sphere closestSphere = currentScene.spheres[objectIndex];

  vec3 origin = ray.origin - closestSphere.Position;
  payload.WorldPosition = origin + ray.direction * hitDistance;
  payload.WorldNormal = normalize(payload.WorldPosition);

  payload.WorldPosition += closestSphere.Position;

  return payload;
}

HitPayload TraceRay(const Ray ray) {
  int closestSphere = -1;
  float hitDistance = FLOAT_MAX;

  for (int i = 0; i <= currentScene.spheres.length(); i++) {
    const Sphere sphere = currentScene.spheres[i];
    vec3 origin = ray.origin - sphere.Position;

    float a = dot(ray.direction, ray.direction);
    float b = 2.0 * dot(origin, ray.direction);
    float c = dot(origin, origin) - sphere.Radius * sphere.Radius;
    
    float discriminant = b * b - 4.0 * a * c;

    if (discriminant < 0.0) continue;

    float closestT = (-b - sqrt(discriminant)) / (2.0 * a);
    if (closestT > 0.0 && closestT < hitDistance) {
      hitDistance = closestT;
      closestSphere = i;
    }
  }

  if (closestSphere < 0) return Miss(ray);

  return ClosestHit(ray, hitDistance, closestSphere);
}

vec4 PerPixel(float x, float y, uint rngState) {
  Ray ray = { cam_origin, (normalize(vec4(vec2(x, y), -1.0, 0.0)) * rotationMatrix).xyz };

	vec3 color = vec3(0.0);
	float multiplier = 1.0;

	int bounces = 5;
	for (int i = 0; i < bounces; i++) {
		HitPayload payload = TraceRay(ray);
		if (payload.HitDistance < 0.0) {
			vec3 skyColor = vec3(0.6, 0.7, 0.9);
			color += skyColor * multiplier;
			break;
		}

		vec3 lightDir = normalize(vec3(-1.0, -1.0, -1.0));
		float lightIntensity = max(dot(payload.WorldNormal, -lightDir), 0.0);

		const Sphere sphere = currentScene.spheres[payload.ObjectIndex];
		const Material material = currentScene.materials[sphere.MaterialIndex];
		vec3 sphereColor = material.Albedo;
		sphereColor *= lightIntensity;
		color += sphereColor * multiplier;

		multiplier *= 0.5;

		ray.origin = payload.WorldPosition + payload.WorldNormal * 0.0001;
		ray.direction = reflect(ray.direction,
			payload.WorldNormal + material.Roughness * RandomUnitVector(rngState));
	}

	return vec4(color, 1.0);
}

void main() {
    vec4 pixel = vec4(0.0);
    ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

    ivec2 dims = imageSize(screen);
    float x = (float(pixel_coords.x * 2 - dims.x) / dims.x);
    float y = (float(pixel_coords.y * 2 - dims.y) / dims.x);

    float fov = 90.0;

    uint rngState = uint(uint(pixel_coords.x) * uint(1973) + uint(pixel_coords.y) * uint(9277) + uint(pixel) * uint(26699)) | uint(1);

    currentScene.spheres[0] = firstSphere;
    currentScene.spheres[1] = secondSphere;

    currentScene.materials[0] = firstMaterial;
    currentScene.materials[1] = secondMaterial;

    vec4 color = PerPixel(x, y, rngState);

    vec4 accumlationData[1280 * 720];

    accumlationData[int(x) + int(y) * 1280] += color;

    vec4 accumlatedColor = accumlationData[int(x) + int(y) * 1280];
    accumlatedColor /= 1.0;

    accumlatedColor = clamp(accumlatedColor, vec4(0.0), vec4(1.0));

    pixel = accumlatedColor;
    imageStore(screen, pixel_coords, pixel);
}