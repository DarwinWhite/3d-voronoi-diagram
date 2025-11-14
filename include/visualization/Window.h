#pragma once

#include <string>
#include <memory>

struct GLFWwindow;

namespace Voronoi3D {

/**
 * @brief Window management class using GLFW
 */
class Window {
private:
    GLFWwindow* window_;
    int width_;
    int height_;
    std::string title_;
    bool vsync_enabled_;

public:
    /**
     * @brief Constructor
     * @param width Window width
     * @param height Window height
     * @param title Window title
     */
    Window(int width = 1200, int height = 800, const std::string& title = "3D Voronoi Visualization");

    /**
     * @brief Destructor
     */
    ~Window();

    // Non-copyable
    Window(const Window&) = delete;
    Window& operator=(const Window&) = delete;

    // Moveable
    Window(Window&& other) noexcept;
    Window& operator=(Window&& other) noexcept;

    /**
     * @brief Initialize the window and OpenGL context
     * @return True if successful
     */
    bool initialize();

    /**
     * @brief Check if window should close
     * @return True if window should close
     */
    bool shouldClose() const;

    /**
     * @brief Swap front and back buffers
     */
    void swapBuffers();

    /**
     * @brief Poll for events
     */
    void pollEvents();

    /**
     * @brief Get window size
     * @param width Output width
     * @param height Output height
     */
    void getSize(int& width, int& height) const;

    /**
     * @brief Get framebuffer size (may differ from window size on high DPI displays)
     * @param width Output width
     * @param height Output height
     */
    void getFramebufferSize(int& width, int& height) const;

    /**
     * @brief Set window size
     * @param width New width
     * @param height New height
     */
    void setSize(int width, int height);

    /**
     * @brief Enable or disable VSync
     * @param enable True to enable VSync
     */
    void setVSync(bool enable);

    /**
     * @brief Get aspect ratio
     * @return Width/height ratio
     */
    float getAspectRatio() const;

    /**
     * @brief Get GLFW window handle
     * @return GLFW window pointer
     */
    GLFWwindow* getGLFWWindow() const { return window_; }

    /**
     * @brief Set window user pointer for callbacks
     * @param ptr User pointer
     */
    void setUserPointer(void* ptr);

    /**
     * @brief Get window user pointer
     * @return User pointer
     */
    void* getUserPointer() const;

    // Callback setters
    void setKeyCallback(void (*callback)(GLFWwindow*, int, int, int, int));
    void setMouseButtonCallback(void (*callback)(GLFWwindow*, int, int, int));
    void setCursorPosCallback(void (*callback)(GLFWwindow*, double, double));
    void setScrollCallback(void (*callback)(GLFWwindow*, double, double));
    void setFramebufferSizeCallback(void (*callback)(GLFWwindow*, int, int));

private:
    void cleanup();
    static void errorCallback(int error, const char* description);
};

} // namespace Voronoi3D