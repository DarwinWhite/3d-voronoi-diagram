#include "visualization/Window.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>
#include <iostream>
#include <stdexcept>

namespace Voronoi3D {

// Static member to track GLFW initialization
static bool glfw_initialized = false;
static int glfw_window_count = 0;

Window::Window(int width, int height, const std::string& title)
    : window_(nullptr), width_(width), height_(height), title_(title), vsync_enabled_(true) {
}

Window::~Window() {
    cleanup();
}

Window::Window(Window&& other) noexcept
    : window_(other.window_), width_(other.width_), height_(other.height_), 
      title_(std::move(other.title_)), vsync_enabled_(other.vsync_enabled_) {
    other.window_ = nullptr;
}

Window& Window::operator=(Window&& other) noexcept {
    if (this != &other) {
        cleanup();
        window_ = other.window_;
        width_ = other.width_;
        height_ = other.height_;
        title_ = std::move(other.title_);
        vsync_enabled_ = other.vsync_enabled_;
        other.window_ = nullptr;
    }
    return *this;
}

bool Window::initialize() {
    // Initialize GLFW if not already done
    if (!glfw_initialized) {
        glfwSetErrorCallback(errorCallback);
        
        if (!glfwInit()) {
            std::cerr << "Failed to initialize GLFW" << std::endl;
            return false;
        }
        glfw_initialized = true;
    }

    // Set OpenGL version (3.3 Core)
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // For macOS compatibility
    
    // Additional hints
    glfwWindowHint(GLFW_SAMPLES, 4); // 4x MSAA
    glfwWindowHint(GLFW_RESIZABLE, GLFW_TRUE);

    // Create window
    window_ = glfwCreateWindow(width_, height_, title_.c_str(), nullptr, nullptr);
    if (!window_) {
        std::cerr << "Failed to create GLFW window" << std::endl;
        return false;
    }

    // Make context current
    glfwMakeContextCurrent(window_);

    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadfunc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD" << std::endl;
        glfwDestroyWindow(window_);
        window_ = nullptr;
        return false;
    }

    // Set VSync
    setVSync(vsync_enabled_);

    // Enable depth testing
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    // Enable multisampling
    glEnable(GL_MULTISAMPLE);

    // Set viewport
    int fb_width, fb_height;
    getFramebufferSize(fb_width, fb_height);
    glViewport(0, 0, fb_width, fb_height);

    ++glfw_window_count;

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << std::endl;
    std::cout << "GLSL Version: " << glGetString(GL_SHADING_LANGUAGE_VERSION) << std::endl;
    
    return true;
}

bool Window::shouldClose() const {
    return window_ && glfwWindowShouldClose(window_);
}

void Window::swapBuffers() {
    if (window_) {
        glfwSwapBuffers(window_);
    }
}

void Window::pollEvents() {
    glfwPollEvents();
}

void Window::getSize(int& width, int& height) const {
    if (window_) {
        glfwGetWindowSize(window_, &width, &height);
    } else {
        width = width_;
        height = height_;
    }
}

void Window::getFramebufferSize(int& width, int& height) const {
    if (window_) {
        glfwGetFramebufferSize(window_, &width, &height);
    } else {
        width = width_;
        height = height_;
    }
}

void Window::setSize(int width, int height) {
    width_ = width;
    height_ = height;
    if (window_) {
        glfwSetWindowSize(window_, width, height);
    }
}

void Window::setVSync(bool enable) {
    vsync_enabled_ = enable;
    if (window_) {
        glfwSwapInterval(enable ? 1 : 0);
    }
}

float Window::getAspectRatio() const {
    int width, height;
    getFramebufferSize(width, height);
    return height > 0 ? static_cast<float>(width) / static_cast<float>(height) : 1.0f;
}

void Window::setUserPointer(void* ptr) {
    if (window_) {
        glfwSetWindowUserPointer(window_, ptr);
    }
}

void* Window::getUserPointer() const {
    return window_ ? glfwGetWindowUserPointer(window_) : nullptr;
}

void Window::setKeyCallback(void (*callback)(GLFWwindow*, int, int, int, int)) {
    if (window_) {
        glfwSetKeyCallback(window_, callback);
    }
}

void Window::setMouseButtonCallback(void (*callback)(GLFWwindow*, int, int, int)) {
    if (window_) {
        glfwSetMouseButtonCallback(window_, callback);
    }
}

void Window::setCursorPosCallback(void (*callback)(GLFWwindow*, double, double)) {
    if (window_) {
        glfwSetCursorPosCallback(window_, callback);
    }
}

void Window::setScrollCallback(void (*callback)(GLFWwindow*, double, double)) {
    if (window_) {
        glfwSetScrollCallback(window_, callback);
    }
}

void Window::setFramebufferSizeCallback(void (*callback)(GLFWwindow*, int, int)) {
    if (window_) {
        glfwSetFramebufferSizeCallback(window_, callback);
    }
}

void Window::cleanup() {
    if (window_) {
        glfwDestroyWindow(window_);
        window_ = nullptr;
        --glfw_window_count;
    }
    
    if (glfw_initialized && glfw_window_count == 0) {
        glfwTerminate();
        glfw_initialized = false;
    }
}

void Window::errorCallback(int error, const char* description) {
    std::cerr << "GLFW Error " << error << ": " << description << std::endl;
}

} // namespace Voronoi3D