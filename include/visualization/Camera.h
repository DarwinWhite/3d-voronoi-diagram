#pragma once

#include "geometry/Point3D.h"
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

namespace Voronoi3D {

/**
 * @brief 3D camera class for OpenGL rendering
 */
class Camera {
public:
    enum class ProjectionType {
        PERSPECTIVE,
        ORTHOGRAPHIC
    };

private:
    // Camera position and orientation
    glm::vec3 position_;
    glm::vec3 target_;
    glm::vec3 up_;
    
    // View parameters
    float distance_;
    float azimuth_;    // Horizontal rotation (around Y axis)
    float elevation_;  // Vertical rotation (around X axis)
    
    // Projection parameters
    ProjectionType projection_type_;
    float fov_;
    float near_plane_;
    float far_plane_;
    float aspect_ratio_;
    
    // Orthographic parameters
    float ortho_size_;
    
    // Movement parameters
    float movement_speed_;
    float rotation_speed_;
    float zoom_speed_;
    
    // Cached matrices
    mutable glm::mat4 view_matrix_;
    mutable glm::mat4 projection_matrix_;
    mutable bool view_matrix_dirty_;
    mutable bool projection_matrix_dirty_;

public:
    /**
     * @brief Constructor
     */
    Camera();

    /**
     * @brief Set camera position
     * @param position Camera position
     */
    void setPosition(const glm::vec3& position);

    /**
     * @brief Set camera target
     * @param target Target point to look at
     */
    void setTarget(const glm::vec3& target);

    /**
     * @brief Set up vector
     * @param up Up direction
     */
    void setUp(const glm::vec3& up);

    /**
     * @brief Look at a specific point
     * @param position Camera position
     * @param target Target to look at
     * @param up Up direction
     */
    void lookAt(const glm::vec3& position, const glm::vec3& target, const glm::vec3& up);

    /**
     * @brief Set orbit parameters (spherical coordinates around target)
     * @param distance Distance from target
     * @param azimuth Horizontal angle in radians
     * @param elevation Vertical angle in radians
     */
    void setOrbit(float distance, float azimuth, float elevation);

    /**
     * @brief Set projection parameters
     * @param fov Field of view in degrees (for perspective)
     * @param aspect_ratio Aspect ratio (width/height)
     * @param near_plane Near clipping plane
     * @param far_plane Far clipping plane
     */
    void setPerspective(float fov, float aspect_ratio, float near_plane, float far_plane);

    /**
     * @brief Set orthographic projection
     * @param size Orthographic view size
     * @param aspect_ratio Aspect ratio
     * @param near_plane Near clipping plane
     * @param far_plane Far clipping plane
     */
    void setOrthographic(float size, float aspect_ratio, float near_plane, float far_plane);

    /**
     * @brief Update aspect ratio (e.g., when window is resized)
     * @param aspect_ratio New aspect ratio
     */
    void setAspectRatio(float aspect_ratio);

    // Movement methods
    /**
     * @brief Orbit around target
     * @param delta_azimuth Change in azimuth (radians)
     * @param delta_elevation Change in elevation (radians)
     */
    void orbit(float delta_azimuth, float delta_elevation);

    /**
     * @brief Pan camera (move target and position together)
     * @param delta_x Horizontal pan
     * @param delta_y Vertical pan
     */
    void pan(float delta_x, float delta_y);

    /**
     * @brief Zoom in/out (change distance from target)
     * @param delta_zoom Zoom delta (positive = zoom in)
     */
    void zoom(float delta_zoom);

    /**
     * @brief Reset camera to default position
     */
    void reset();

    // Matrix accessors
    /**
     * @brief Get view matrix
     * @return View matrix
     */
    const glm::mat4& getViewMatrix() const;

    /**
     * @brief Get projection matrix
     * @return Projection matrix
     */
    const glm::mat4& getProjectionMatrix() const;

    /**
     * @brief Get view-projection matrix
     * @return Combined view-projection matrix
     */
    glm::mat4 getViewProjectionMatrix() const;

    // Property accessors
    glm::vec3 getPosition() const { return position_; }
    glm::vec3 getTarget() const { return target_; }
    glm::vec3 getUp() const { return up_; }
    float getDistance() const { return distance_; }
    float getAzimuth() const { return azimuth_; }
    float getElevation() const { return elevation_; }
    float getFOV() const { return fov_; }
    float getNearPlane() const { return near_plane_; }
    float getFarPlane() const { return far_plane_; }
    ProjectionType getProjectionType() const { return projection_type_; }

    // Speed settings
    void setMovementSpeed(float speed) { movement_speed_ = speed; }
    void setRotationSpeed(float speed) { rotation_speed_ = speed; }
    void setZoomSpeed(float speed) { zoom_speed_ = speed; }

    float getMovementSpeed() const { return movement_speed_; }
    float getRotationSpeed() const { return rotation_speed_; }
    float getZoomSpeed() const { return zoom_speed_; }

    // Utility methods
    /**
     * @brief Convert screen coordinates to world ray
     * @param screen_x Screen X coordinate (0 to width)
     * @param screen_y Screen Y coordinate (0 to height)
     * @param screen_width Screen width
     * @param screen_height Screen height
     * @return Ray direction in world space
     */
    glm::vec3 screenToWorldRay(float screen_x, float screen_y, 
                               float screen_width, float screen_height) const;

    /**
     * @brief Get camera forward direction
     * @return Forward vector
     */
    glm::vec3 getForward() const;

    /**
     * @brief Get camera right direction
     * @return Right vector
     */
    glm::vec3 getRight() const;

private:
    void updateViewMatrix() const;
    void updateProjectionMatrix() const;
    void updateFromOrbit();
    void clampElevation();
};

} // namespace Voronoi3D