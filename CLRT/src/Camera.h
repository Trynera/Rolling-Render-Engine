#pragma once

struct Camera {
    glm::vec3 position;
    glm::vec3 lookAt;
    glm::vec3 up;
    glm::vec3 forwardDir;
    float cameraYaw;
    float cameraPitch;
};