#version 460 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image2D screen;

#define MAX_NUMBER_OF_OBJECTS 2
#define FLOAT_MAX 340282346638528859811704183484516925440.0

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
  vec3 SphereColor;
  float Roughness;
};

struct Scene {
  Sphere spheres[MAX_NUMBER_OF_OBJECTS];
};

layout(location = 1) uniform vec3 cam_origin;
layout(location = 2) uniform mat4 cam_view;
layout(location = 3) uniform mat4 cam_proj;
layout(location = 4) uniform vec3 cam_fdir;
layout(location = 5) uniform vec3 cam_up;
layout(location = 6) uniform mat4 rotationMatrix;
layout(location = 7) uniform float deltaTime;

uint frameIndex = 1;

float hitDistance = FLOAT_MAX;

Scene currentScene;
Sphere firstSphere = { vec3(0.0, 0.0, 0.0), 1.0, vec3(1.0, 0.0, 1.0), 0.0 };
Sphere secondSphere = { vec3(0.0, -101.0, 0.0), 100.0, vec3(0.2, 0.3, 1.0), 0.1 };

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
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

    for (int i = 0; i <= 1; i++) {
      const Sphere sphere = currentScene.spheres[i];
      vec3 origin = ray.origin - sphere.Position;

      float a = dot(ray.direction, ray.direction);
      float b = 2.0 * dot(origin, ray.direction);
      float c = dot(origin, origin) - sphere.Radius * sphere.Radius;
    
      float discriminant = b * b - 4.0 * a * c;
      if (discriminant >= 0.0) {
        // float t0 = (-b + sqrt(discriminant)) / (2.0 * a);
        float closestT = (-b - sqrt(discriminant)) / (2.0 * a);
        if (closestT > 0.0 && closestT < hitDistance) {
          hitDistance = closestT;
          closestSphere = i;
        }

        vec3 hitPoint = origin + ray.direction * closestT;
        vec3 normal = normalize(hitPoint);

        vec3 lightDir = normalize(vec3(-1.0, -1.0, -1.0));

        float d = max(dot(normal, -lightDir), 0.0);

        vec3 sphereColor = sphere.SphereColor;
        sphereColor *= d;
      }
    }

    if (closestSphere < 0) return Miss(ray);

    return ClosestHit(ray, hitDistance, closestSphere);
}

void main() {
    vec4 pixel = vec4(0.0);
    ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

    ivec2 dims = imageSize(screen);
    float x = (float(pixel_coords.x * 2 - dims.x) / dims.x);
    float y = (float(pixel_coords.y * 2 - dims.y) / dims.x);

    float fov = 90.0;
    Ray ray = { cam_origin, (normalize(vec4(vec2(x, y), -1.0, 0.0)) * rotationMatrix).xyz };

    currentScene.spheres[0] = firstSphere;
    currentScene.spheres[1] = secondSphere;

    float multiplier = 1.0;

    vec3 color = vec3(0.0);

    int bounces = 5;
    for (int i = 0; i < bounces; i++) {
      HitPayload payload = TraceRay(ray);
      if (payload.HitDistance < 0.0) {
        vec3 skyColor = vec3(0.6, 0.7, 0.9);
        color += skyColor * multiplier;
        break;
      }

      vec3 lightDir = normalize(vec3(-1.0));
      float lightIntensity = max(dot(payload.WorldNormal, -lightDir), 0.0);

      const Sphere sphere = currentScene.spheres[payload.ObjectIndex];
      vec3 sphereColor = sphere.SphereColor;
      sphereColor *= lightIntensity;
      color += sphereColor * multiplier;

      multiplier *= 0.5;

      ray.origin = payload.WorldPosition + payload.WorldNormal * 0.0001;
      ray.direction = reflect(ray.direction, payload.WorldNormal + sphere.Roughness * vec3(rand(vec2(-0.5, 0.5))));
    }

    vec4 accumlationData[1280 * 720];
    accumlationData[int(x) + int(y)] += vec4(color, 1.0);

    vec4 accumlatedColor = accumlationData[int(x) + int(y)];
    accumlatedColor /= float(frameIndex);

    accumlatedColor = clamp(accumlatedColor, vec4(0.0), vec4(1.0));

    pixel = accumlatedColor;
    imageStore(screen, pixel_coords, pixel);

    frameIndex++;
}