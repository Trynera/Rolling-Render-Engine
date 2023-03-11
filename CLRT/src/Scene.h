#pragma once

#include <glm.hpp>

#include <vector>

struct Material {
	glm::vec3 Albedo{ 1.0f };
	float Fuzz = 0.0f;
	float IR = 0.0f;
	int Type = 0;
};

struct Sphere {
	glm::vec3 Position{ 0.0f };
	float Radius = 0.5f;

	int MaterialIndex = 0;
};