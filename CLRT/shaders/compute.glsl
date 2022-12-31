#version 460 core
layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;
layout(rgba32f, binding = 0) uniform image2D screen;

#define MAX_NUMBER_OF_OBJECTS 2
#define FLOAT_MAX 340282346638528859811704183484516925440.0

struct HitPayload {
  float hitDistance;
  vec3 WorldPosition;
  vec3 WorldNormal;

  int ObjectIndex;
};

struct Sphere {
  vec3 Position;
	float Radius;
  vec3 SphereColor;
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

float hitDistance = FLOAT_MAX;

void main() {
    vec4 pixel = vec4(0.075, 0.133, 0.173, 1.0);
    ivec2 pixel_coords = ivec2(gl_GlobalInvocationID.xy);

    ivec2 dims = imageSize(screen);
    float x = (float(pixel_coords.x * 2 - dims.x) / dims.x);
    float y = (float(pixel_coords.y * 2 - dims.y) / dims.x);

    float fov = 90.0;
    vec3 ray_o = cam_origin;
    vec3 ray_d = (normalize(vec4(vec2(x, y), -1.0, 0.0)) * rotationMatrix).xyz;

    Scene currentScene;
    Sphere firstSphere = { vec3(0.0, 0.0, 0.0), 1.0, vec3(1.0, 0.0, 1.0) };
    Sphere secondSphere = { vec3(0.0, -101.0, 0.0), 100.0, vec3(0.2, 0.3, 1.0) };
    currentScene.spheres[0] = firstSphere;
    currentScene.spheres[1] = secondSphere;

    for (int i = 0; i <= 1; i++) {
      const Sphere sphere = currentScene.spheres[i];
      vec3 origin = ray_o - sphere.Position;

      float a = dot(ray_d, ray_d);
      float b = 2.0 * dot(origin, ray_d);
      float c = dot(origin, origin) - sphere.Radius * sphere.Radius;
    
      float discriminant = b * b - 4.0 * a * c;
      if (discriminant >= 0.0) {
        // float t0 = (-b + sqrt(discriminant)) / (2.0 * a);
        float closestT = (-b - sqrt(discriminant)) / (2.0 * a);

        vec3 hitPoint = origin + ray_d * closestT;
        vec3 normal = normalize(hitPoint);

        vec3 lightDir = normalize(vec3(-1.0, -1.0, -1.0));

        float d = max(dot(normal, -lightDir), 0.0);

        vec3 sphereColor = sphere.SphereColor;
        sphereColor *= d;

        pixel = vec4(sphereColor, 1.0);
      }

    imageStore(screen, pixel_coords, pixel);
    }
}