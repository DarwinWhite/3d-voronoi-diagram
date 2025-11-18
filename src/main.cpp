#include <iostream>
#include <vector>
#include <memory>
#include <random>

#include "visualization/Window.h"
#include "visualization/Camera.h"
#include "geometry/Point3D.h"
#include "algorithms/VoronoiStructures.h"

#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

namespace Voronoi3D {

/**
 * @brief Main application class
 */
class VoronoiApp {
private:
    std::unique_ptr<Window> window_;
    std::unique_ptr<Camera> camera_;
    std::unique_ptr<VoronoiDiagram> voronoi_diagram_;
    
    // Input state
    bool mouse_pressed_;
    double last_mouse_x_;
    double last_mouse_y_;
    
    // Rendering state
    std::vector<Point3D> test_points_;
    GLuint shader_program_;
    GLuint vao_;
    GLuint vbo_;
    
public:
    VoronoiApp() : mouse_pressed_(false), last_mouse_x_(0.0), last_mouse_y_(0.0),
                   shader_program_(0), vao_(0), vbo_(0) {
    }
    
    ~VoronoiApp() {
        cleanup();
    }
    
    bool initialize() {
        // Create window
        window_ = std::make_unique<Window>(1200, 800, "3D Voronoi Diagram Visualization");
        if (!window_->initialize()) {
            std::cerr << "Failed to initialize window" << std::endl;
            return false;
        }
        
        // Setup callbacks
        window_->setUserPointer(this);
        window_->setKeyCallback(keyCallback);
        window_->setMouseButtonCallback(mouseButtonCallback);
        window_->setCursorPosCallback(cursorPosCallback);
        window_->setScrollCallback(scrollCallback);
        window_->setFramebufferSizeCallback(framebufferSizeCallback);
        
        // Create camera
        camera_ = std::make_unique<Camera>();
        camera_->setPerspective(45.0f, window_->getAspectRatio(), 0.1f, 100.0f);
        camera_->setOrbit(10.0f, 0.0f, 0.3f); // Distance, azimuth, elevation
        
        // Create Voronoi diagram
        voronoi_diagram_ = std::make_unique<VoronoiDiagram>();
        
        // Setup OpenGL
        if (!setupOpenGL()) {
            std::cerr << "Failed to setup OpenGL" << std::endl;
            return false;
        }
        
        // Generate test data
        generateTestData();
        
        std::cout << "Application initialized successfully!" << std::endl;
        std::cout << "Controls:" << std::endl;
        std::cout << "  Left mouse + drag: Rotate camera" << std::endl;
        std::cout << "  Right mouse + drag: Pan camera" << std::endl;
        std::cout << "  Mouse wheel: Zoom" << std::endl;
        std::cout << "  Space: Generate new random points" << std::endl;
        std::cout << "  R: Reset camera" << std::endl;
        std::cout << "  ESC: Exit" << std::endl;
        
        return true;
    }
    
    void run() {
        while (!window_->shouldClose()) {
            window_->pollEvents();
            
            update();
            render();
            
            window_->swapBuffers();
        }
    }
    
private:
    bool setupOpenGL() {
        // Create basic shader program
        const char* vertexShaderSource = R"(
            #version 330 core
            layout (location = 0) in vec3 aPos;
            uniform mat4 uMVP;
            void main() {
                gl_Position = uMVP * vec4(aPos, 1.0);
                gl_PointSize = 8.0;
            }
        )";
        
        const char* fragmentShaderSource = R"(
            #version 330 core
            out vec4 FragColor;
            void main() {
                FragColor = vec4(1.0, 0.5, 0.2, 1.0);
            }
        )";
        
        // Compile shaders
        GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vertexShaderSource);
        GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fragmentShaderSource);
        
        if (vertexShader == 0 || fragmentShader == 0) {
            return false;
        }
        
        // Create program
        shader_program_ = glCreateProgram();
        glAttachShader(shader_program_, vertexShader);
        glAttachShader(shader_program_, fragmentShader);
        glLinkProgram(shader_program_);
        
        // Check linking
        GLint success;
        glGetProgramiv(shader_program_, GL_LINK_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetProgramInfoLog(shader_program_, 512, nullptr, infoLog);
            std::cerr << "Shader program linking failed: " << infoLog << std::endl;
            return false;
        }
        
        glDeleteShader(vertexShader);
        glDeleteShader(fragmentShader);
        
        // Setup VAO and VBO
        glGenVertexArrays(1, &vao_);
        glGenBuffers(1, &vbo_);
        
        glBindVertexArray(vao_);
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        
        // Position attribute
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
        glEnableVertexAttribArray(0);
        
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
        
        // Enable point rendering
        glEnable(GL_PROGRAM_POINT_SIZE);
        
        return true;
    }
    
    GLuint compileShader(GLenum type, const char* source) {
        GLuint shader = glCreateShader(type);
        glShaderSource(shader, 1, &source, nullptr);
        glCompileShader(shader);
        
        GLint success;
        glGetShaderiv(shader, GL_COMPILE_STATUS, &success);
        if (!success) {
            char infoLog[512];
            glGetShaderInfoLog(shader, 512, nullptr, infoLog);
            std::cerr << "Shader compilation failed: " << infoLog << std::endl;
            return 0;
        }
        
        return shader;
    }
    
    void generateTestData() {
        test_points_.clear();
        
        // Generate random points in a cube
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<double> dist(-3.0, 3.0);
        
        const size_t num_points = 8; // Smaller number for better visualization
        for (size_t i = 0; i < num_points; ++i) {
            test_points_.emplace_back(dist(gen), dist(gen), dist(gen));
        }
        
        std::cout << "Generated " << test_points_.size() << " test points" << std::endl;
        
        // Compute Voronoi diagram
        try {
            voronoi_diagram_->compute(test_points_);
            std::cout << "Voronoi diagram computed successfully!" << std::endl;
            std::cout << voronoi_diagram_->getStatistics() << std::endl;
        } catch (const std::exception& e) {
            std::cout << "Failed to compute Voronoi diagram: " << e.what() << std::endl;
            std::cout << "Continuing with point visualization only." << std::endl;
        }
        
        // Update VBO with new data
        updatePointBuffer();
    }
    
    void updatePointBuffer() {
        std::vector<float> vertices;
        for (const auto& point : test_points_) {
            vertices.push_back(static_cast<float>(point.x));
            vertices.push_back(static_cast<float>(point.y));
            vertices.push_back(static_cast<float>(point.z));
        }
        
        glBindBuffer(GL_ARRAY_BUFFER, vbo_);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }
    
    void update() {
        // Application logic updates would go here
    }
    
    void render() {
        // Clear screen
        glClearColor(0.1f, 0.1f, 0.15f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        
        // Use shader program
        glUseProgram(shader_program_);
        
        // Set MVP matrix
        glm::mat4 model = glm::mat4(1.0f);
        glm::mat4 mvp = camera_->getViewProjectionMatrix() * model;
        
        GLint mvpLocation = glGetUniformLocation(shader_program_, "uMVP");
        glUniformMatrix4fv(mvpLocation, 1, GL_FALSE, glm::value_ptr(mvp));
        
        // Draw points
        glBindVertexArray(vao_);
        glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(test_points_.size()));
        glBindVertexArray(0);
        
        glUseProgram(0);
    }
    
    void cleanup() {
        if (vao_ != 0) {
            glDeleteVertexArrays(1, &vao_);
            vao_ = 0;
        }
        if (vbo_ != 0) {
            glDeleteBuffers(1, &vbo_);
            vbo_ = 0;
        }
        if (shader_program_ != 0) {
            glDeleteProgram(shader_program_);
            shader_program_ = 0;
        }
    }
    
    // Static callback functions
    static void keyCallback(GLFWwindow* window, int key, int scancode, int action, int mods) {
        VoronoiApp* app = static_cast<VoronoiApp*>(glfwGetWindowUserPointer(window));
        app->handleKey(key, scancode, action, mods);
    }
    
    static void mouseButtonCallback(GLFWwindow* window, int button, int action, int mods) {
        VoronoiApp* app = static_cast<VoronoiApp*>(glfwGetWindowUserPointer(window));
        app->handleMouseButton(button, action, mods);
    }
    
    static void cursorPosCallback(GLFWwindow* window, double xpos, double ypos) {
        VoronoiApp* app = static_cast<VoronoiApp*>(glfwGetWindowUserPointer(window));
        app->handleCursorPos(xpos, ypos);
    }
    
    static void scrollCallback(GLFWwindow* window, double xoffset, double yoffset) {
        VoronoiApp* app = static_cast<VoronoiApp*>(glfwGetWindowUserPointer(window));
        app->handleScroll(xoffset, yoffset);
    }
    
    static void framebufferSizeCallback(GLFWwindow* window, int width, int height) {
        glViewport(0, 0, width, height);
        VoronoiApp* app = static_cast<VoronoiApp*>(glfwGetWindowUserPointer(window));
        app->handleFramebufferSize(width, height);
    }
    
    // Event handlers
    void handleKey(int key, int /*scancode*/, int action, int /*mods*/) {
        if (action == GLFW_PRESS) {
            switch (key) {
                case GLFW_KEY_ESCAPE:
                    glfwSetWindowShouldClose(window_->getGLFWWindow(), GLFW_TRUE);
                    break;
                case GLFW_KEY_SPACE:
                    generateTestData();
                    break;
                case GLFW_KEY_R:
                    camera_->reset();
                    break;
            }
        }
    }
    
    void handleMouseButton(int button, int action, int /*mods*/) {
        if (button == GLFW_MOUSE_BUTTON_LEFT || button == GLFW_MOUSE_BUTTON_RIGHT) {
            if (action == GLFW_PRESS) {
                mouse_pressed_ = true;
                glfwGetCursorPos(window_->getGLFWWindow(), &last_mouse_x_, &last_mouse_y_);
            } else if (action == GLFW_RELEASE) {
                mouse_pressed_ = false;
            }
        }
    }
    
    void handleCursorPos(double xpos, double ypos) {
        if (mouse_pressed_) {
            double dx = xpos - last_mouse_x_;
            double dy = ypos - last_mouse_y_;
            
            if (glfwGetMouseButton(window_->getGLFWWindow(), GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS) {
                // Orbit camera
                camera_->orbit(static_cast<float>(dx) * 0.01f, static_cast<float>(-dy) * 0.01f);
            } else if (glfwGetMouseButton(window_->getGLFWWindow(), GLFW_MOUSE_BUTTON_RIGHT) == GLFW_PRESS) {
                // Pan camera
                camera_->pan(static_cast<float>(-dx) * 0.01f, static_cast<float>(dy) * 0.01f);
            }
            
            last_mouse_x_ = xpos;
            last_mouse_y_ = ypos;
        }
    }
    
    void handleScroll(double /*xoffset*/, double yoffset) {
        camera_->zoom(static_cast<float>(yoffset));
    }
    
    void handleFramebufferSize(int width, int height) {
        camera_->setAspectRatio(static_cast<float>(width) / static_cast<float>(height));
    }
};

} // namespace Voronoi3D

int main() {
    try {
        Voronoi3D::VoronoiApp app;
        
        if (!app.initialize()) {
            std::cerr << "Failed to initialize application" << std::endl;
            return -1;
        }
        
        app.run();
        
    } catch (const std::exception& e) {
        std::cerr << "Exception: " << e.what() << std::endl;
        return -1;
    }
    
    return 0;
}