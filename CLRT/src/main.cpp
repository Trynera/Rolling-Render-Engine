// compute shaders tutorial
// Dr Anton Gerdelan <gerdela@scss.tcd.ie>
// Trinity College Dublin, Ireland
// 26 Feb 2016. latest v 2 Mar 2016

#include "gl_utils.h"
#include "Scene.h"
#include "Camera.h"
#include "Ray.h"

#include <gtc/type_ptr.hpp>
#include <gtc/quaternion.hpp>
#include <gtx/quaternion.hpp>

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

Camera camera;
glm::mat4 projection;

bool mouseAbsorbed = false;
bool refreshRequired = false;

glm::mat4 rotationMatrix(1);

void keyCallback(GLFWwindow* window, int button, int action, int mods) {
    if (action == GLFW_PRESS) {
        if (button == GLFW_MOUSE_BUTTON_RIGHT) {
            glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
            glfwSetCursorPos(window, 1280.0f / 2.0f, 720.0f / 2.0f);
            mouseAbsorbed = true;
        }
    }
    else {
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
        mouseAbsorbed = false;
    }
}

// Thanks carl-vbn lmao
bool handleInput(GLFWwindow* window, double deltaTime, glm::vec3& cameraPosition, float& cameraYaw, float& cameraPitch, glm::mat4* rotationMatrix) {
    bool moved = false;

    double mouseX;
    double mouseY;
    glfwGetCursorPos(window, &mouseX, &mouseY);
    glfwSetCursorPos(window, 1280.0 / 2.0, 720.0 / 2.0);

    float xOffset = (float)(mouseX - 1280.0 / 2.0);
    float yOffset = (float)(mouseY - 720.0 / 2.0);

    if (xOffset != 0.0f || yOffset != 0.0f) moved = true;

    cameraYaw += xOffset * 0.002f;
    cameraPitch += yOffset * 0.002f;

    if (cameraPitch > 1.5707f)
        cameraPitch = 1.5707f;
    if (cameraPitch < -1.5707f)
        cameraPitch = -1.5707f;

    *rotationMatrix = glm::rotate(glm::rotate(glm::mat4(1), cameraPitch, glm::vec3(1.0f, 0.0f, 0.0f)), cameraYaw, glm::vec3(0.0f, 1.0f, 0.0f));

    glm::vec3 forward = glm::vec3(glm::vec4(0.0f, 0.0f, -1.0f, 0.0f) * (*rotationMatrix));

    glm::vec3 up(0.0f, 1.0f, 0.0f);
    glm::vec3 right = glm::cross(forward, up);

    glm::vec3 movementDirection(0);
    float multiplier = 1;

    if (glfwGetKey(window, GLFW_KEY_W)) {
        movementDirection += forward;
    }
    if (glfwGetKey(window, GLFW_KEY_S)) {
        movementDirection -= forward;
    }
    if (glfwGetKey(window, GLFW_KEY_D)) {
        movementDirection += right;
    }
    if (glfwGetKey(window, GLFW_KEY_A)) {
        movementDirection -= right;
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE)) {
        movementDirection += up;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT)) {
        movementDirection -= up;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_CONTROL)) {
        multiplier = 2.5f;
    }

    if (glm::length(movementDirection) > 0.0f) {
        cameraPosition += glm::normalize(movementDirection) * (float)deltaTime * (float)multiplier;
        moved = true;
    }

    return moved;
}

// Code stolen from https://antongerdelan.net/opengl/glcontext2.html, thanks bud.
void _update_fps_counter(GLFWwindow* window) {
    static double previous_seconds = glfwGetTime();
    static int frame_count;
    double current_seconds = glfwGetTime();
    double elapsed_seconds = current_seconds - previous_seconds;
    if (elapsed_seconds > 0.25) {
        previous_seconds = current_seconds;
        double fps = (double)frame_count / elapsed_seconds;
        char tmp[128];
        sprintf_s(tmp, "OpenGL Ray Tracer @ FPS: %.2f", fps);
        glfwSetWindowTitle(window, tmp);
        frame_count = 0;
    }
    frame_count++;
}

int main() {
    (start_gl()); // just starts a 4.3 GL context+window

    projection = glm::perspective(90.0f, 1280.0f / 720.0f, 0.1f, 100.0f);

    std::string computeCode;
    std::ifstream cShaderFile;
    cShaderFile.exceptions(std::ifstream::failbit | std::ifstream::badbit);

    try {
        cShaderFile.open("shaders\\compute.glsl");
        std::stringstream cShaderStream;
        cShaderStream << cShaderFile.rdbuf();
        cShaderFile.close();
        computeCode = cShaderStream.str();
    }
    catch (std::ifstream::failure& e) {
        std::cout << "ERROR::SHADER::FILE_NOT_SUCCESFULLY_READ: " << e.what() << std::endl;
    }
    const char* compute_shader_str = computeCode.c_str();

    // set up shaders and geometry for full-screen quad
    // moved code to gl_utils.cpp
    GLuint quad_vao = create_quad_vao();
    GLuint quad_program = create_quad_program();

    GLuint ray_program = 0;
    { // create the compute shader
        GLuint ray_shader = glCreateShader(GL_COMPUTE_SHADER);
        glShaderSource(ray_shader, 1, &compute_shader_str, NULL);
        glCompileShader(ray_shader);
        (check_shader_errors(ray_shader)); // code moved to gl_utils.cpp
        ray_program = glCreateProgram();
        glAttachShader(ray_program, ray_shader);
        glLinkProgram(ray_program);
        (check_program_errors(ray_program)); // code moved to gl_utils.cpp
    }

    // texture handle and dimensions
    GLuint tex_output = 0;
    int tex_w = 1280, tex_h = 720;
    { // create the texture
        glGenTextures(1, &tex_output);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex_output);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
        // linear allows us to scale the window up retaining reasonable quality
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        // same internal format as compute shader input
        glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA32F, tex_w, tex_h, 0, GL_RGBA, GL_FLOAT, NULL);
        // bind to image unit so can write to specific pixels from the shader
        glBindImageTexture(0, tex_output, 0, GL_FALSE, 0, GL_WRITE_ONLY, GL_RGBA32F);
    }

    { // query up the workgroups
        int work_grp_size[3], work_grp_inv;
        // maximum global work group (total work in a dispatch)
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &work_grp_size[0]);
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &work_grp_size[1]);
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &work_grp_size[2]);
        printf("max global (total) work group size x:%i y:%i z:%i\n", work_grp_size[0], work_grp_size[1], work_grp_size[2]);
        // maximum local work group (one shader's slice)
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 0, &work_grp_size[0]);
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 1, &work_grp_size[1]);
        glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_SIZE, 2, &work_grp_size[2]);
        printf("max local (in one shader) work group sizes x:%i y:%i z:%i\n", work_grp_size[0], work_grp_size[1], work_grp_size[2]);
        // maximum compute shader invocations (x * y * z)
        glGetIntegerv(GL_MAX_COMPUTE_WORK_GROUP_INVOCATIONS, &work_grp_inv);
        printf("max computer shader invocations %i\n", work_grp_inv);
    }

    Camera camera = { glm::vec3(0.0f, 0.0f, 2.0f), glm::vec3(0.0f), glm::vec3(0.0f), glm::vec3(0.0f, 0.0f, -1.0f), 0.0f, 0.0f };

    glm::mat4 view = glm::lookAt(camera.position, camera.lookAt, camera.up);
    
    float deltaTime = 0.0f;
    float lastFrame = 0.0f;

    // Loops until the Window is closed
    while (!glfwWindowShouldClose(window)) {
        _update_fps_counter(window);

        view = glm::lookAt(camera.position, camera.lookAt, camera.up);

        {                                          // launch compute shaders!
            glUseProgram(ray_program);
            glDispatchCompute(ceil((GLuint)tex_w / 8), ceil((GLuint)tex_h / 8), 1);
        }
        
        {
            glUniform3f(1, camera.position.x, camera.position.y, camera.position.z);
            glUniformMatrix4fv(2, 1, GL_FALSE, glm::value_ptr(view));
            glUniformMatrix4fv(3, 1, GL_FALSE, glm::value_ptr(projection));
            glUniform3f(4, camera.forwardDir.x, camera.forwardDir.y, camera.forwardDir.z);
            glUniform3f(5, camera.up.x, camera.up.y, camera.up.z);
            glUniformMatrix4fv(6, 1, GL_FALSE, glm::value_ptr(rotationMatrix));
        }

        glMemoryBarrier(GL_ALL_BARRIER_BITS);

        // This is responsible for drawing the Texture we rendered with the Compute Shader
        glClear(GL_COLOR_BUFFER_BIT);
        glUseProgram(quad_program);
        glBindVertexArray(quad_vao);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, tex_output);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

        // Checks for PollEvents
        glfwPollEvents();

        // The Window just closes when Escape is pressed
        if (GLFW_PRESS == glfwGetKey(window, GLFW_KEY_ESCAPE)) { glfwSetWindowShouldClose(window, 1); }

        float currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;

        glfwSetMouseButtonCallback(window, keyCallback);

        if (mouseAbsorbed) {
            if (handleInput(window, deltaTime, camera.position, camera.cameraYaw, camera.cameraPitch, &rotationMatrix)) {
                refreshRequired = true;
            }
        }
        else {
            rotationMatrix = glm::rotate(glm::rotate(glm::mat4(1), camera.cameraPitch, glm::vec3(1.0f, 0.0f, 0.0f)), camera.cameraYaw, glm::vec3(0.0f, 1.0f, 0.0f));
        }

        // Swaps the Buffers for the Window we just created
        glfwSwapBuffers(window);
    }

    stop_gl();
    return 0;
}