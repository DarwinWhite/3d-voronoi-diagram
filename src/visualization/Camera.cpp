#include "visualization/Camera.h"
#include <algorithm>
#include <cmath>

namespace Voronoi3D {

Camera::Camera()
    : position_(0.0f, 0.0f, 5.0f)
    , target_(0.0f, 0.0f, 0.0f)
    , up_(0.0f, 1.0f, 0.0f)
    , distance_(5.0f)
    , azimuth_(0.0f)
    , elevation_(0.0f)
    , projection_type_(ProjectionType::PERSPECTIVE)
    , fov_(45.0f)
    , near_plane_(0.1f)
    , far_plane_(100.0f)
    , aspect_ratio_(16.0f / 9.0f)
    , ortho_size_(5.0f)
    , movement_speed_(1.0f)
    , rotation_speed_(1.0f)
    , zoom_speed_(1.0f)
    , view_matrix_dirty_(true)
    , projection_matrix_dirty_(true) {
}

void Camera::setPosition(const glm::vec3& position) {
    position_ = position;
    view_matrix_dirty_ = true;
    
    // Update orbit parameters
    glm::vec3 to_target = target_ - position_;
    distance_ = glm::length(to_target);
    if (distance_ > 0.001f) {
        to_target /= distance_;
        azimuth_ = std::atan2(to_target.x, to_target.z);
        elevation_ = std::asin(to_target.y);
    }
}

void Camera::setTarget(const glm::vec3& target) {
    target_ = target;
    view_matrix_dirty_ = true;
    
    // Update orbit parameters
    glm::vec3 to_target = target_ - position_;
    distance_ = glm::length(to_target);
    if (distance_ > 0.001f) {
        to_target /= distance_;
        azimuth_ = std::atan2(to_target.x, to_target.z);
        elevation_ = std::asin(to_target.y);
    }
}

void Camera::setUp(const glm::vec3& up) {
    up_ = glm::normalize(up);
    view_matrix_dirty_ = true;
}

void Camera::lookAt(const glm::vec3& position, const glm::vec3& target, const glm::vec3& up) {
    position_ = position;
    target_ = target;
    up_ = glm::normalize(up);
    view_matrix_dirty_ = true;
    
    // Update orbit parameters
    glm::vec3 to_target = target_ - position_;
    distance_ = glm::length(to_target);
    if (distance_ > 0.001f) {
        to_target /= distance_;
        azimuth_ = std::atan2(to_target.x, to_target.z);
        elevation_ = std::asin(to_target.y);
    }
}

void Camera::setOrbit(float distance, float azimuth, float elevation) {
    distance_ = std::max(0.1f, distance);
    azimuth_ = azimuth;
    elevation_ = elevation;
    clampElevation();
    updateFromOrbit();
    view_matrix_dirty_ = true;
}

void Camera::setPerspective(float fov, float aspect_ratio, float near_plane, float far_plane) {
    projection_type_ = ProjectionType::PERSPECTIVE;
    fov_ = fov;
    aspect_ratio_ = aspect_ratio;
    near_plane_ = near_plane;
    far_plane_ = far_plane;
    projection_matrix_dirty_ = true;
}

void Camera::setOrthographic(float size, float aspect_ratio, float near_plane, float far_plane) {
    projection_type_ = ProjectionType::ORTHOGRAPHIC;
    ortho_size_ = size;
    aspect_ratio_ = aspect_ratio;
    near_plane_ = near_plane;
    far_plane_ = far_plane;
    projection_matrix_dirty_ = true;
}

void Camera::setAspectRatio(float aspect_ratio) {
    aspect_ratio_ = aspect_ratio;
    projection_matrix_dirty_ = true;
}

void Camera::orbit(float delta_azimuth, float delta_elevation) {
    azimuth_ += delta_azimuth * rotation_speed_;
    elevation_ += delta_elevation * rotation_speed_;
    clampElevation();
    updateFromOrbit();
    view_matrix_dirty_ = true;
}

void Camera::pan(float delta_x, float delta_y) {
    glm::vec3 right = getRight();
    glm::vec3 up = getUp();
    
    glm::vec3 offset = (right * delta_x + up * delta_y) * movement_speed_;
    position_ += offset;
    target_ += offset;
    view_matrix_dirty_ = true;
}

void Camera::zoom(float delta_zoom) {
    distance_ -= delta_zoom * zoom_speed_;
    distance_ = std::max(0.1f, distance_);
    updateFromOrbit();
    view_matrix_dirty_ = true;
}

void Camera::reset() {
    distance_ = 5.0f;
    azimuth_ = 0.0f;
    elevation_ = 0.0f;
    target_ = glm::vec3(0.0f, 0.0f, 0.0f);
    up_ = glm::vec3(0.0f, 1.0f, 0.0f);
    updateFromOrbit();
    view_matrix_dirty_ = true;
}

const glm::mat4& Camera::getViewMatrix() const {
    if (view_matrix_dirty_) {
        updateViewMatrix();
    }
    return view_matrix_;
}

const glm::mat4& Camera::getProjectionMatrix() const {
    if (projection_matrix_dirty_) {
        updateProjectionMatrix();
    }
    return projection_matrix_;
}

glm::mat4 Camera::getViewProjectionMatrix() const {
    return getProjectionMatrix() * getViewMatrix();
}

glm::vec3 Camera::screenToWorldRay(float screen_x, float screen_y, 
                                   float screen_width, float screen_height) const {
    // Convert to normalized device coordinates
    float ndc_x = (2.0f * screen_x) / screen_width - 1.0f;
    float ndc_y = 1.0f - (2.0f * screen_y) / screen_height;
    
    // Get inverse view-projection matrix
    glm::mat4 inv_vp = glm::inverse(getViewProjectionMatrix());
    
    // Ray in world space
    glm::vec4 ray_world = inv_vp * glm::vec4(ndc_x, ndc_y, -1.0f, 1.0f);
    ray_world /= ray_world.w;
    
    glm::vec3 ray_direction = glm::normalize(glm::vec3(ray_world) - position_);
    return ray_direction;
}

glm::vec3 Camera::getForward() const {
    return glm::normalize(target_ - position_);
}

glm::vec3 Camera::getRight() const {
    return glm::normalize(glm::cross(getForward(), up_));
}

void Camera::updateViewMatrix() const {
    view_matrix_ = glm::lookAt(position_, target_, up_);
    view_matrix_dirty_ = false;
}

void Camera::updateProjectionMatrix() const {
    if (projection_type_ == ProjectionType::PERSPECTIVE) {
        projection_matrix_ = glm::perspective(glm::radians(fov_), aspect_ratio_, near_plane_, far_plane_);
    } else {
        float half_width = ortho_size_ * aspect_ratio_ * 0.5f;
        float half_height = ortho_size_ * 0.5f;
        projection_matrix_ = glm::ortho(-half_width, half_width, -half_height, half_height, near_plane_, far_plane_);
    }
    projection_matrix_dirty_ = false;
}

void Camera::updateFromOrbit() {
    // Convert spherical coordinates to Cartesian
    float x = distance_ * std::sin(azimuth_) * std::cos(elevation_);
    float y = distance_ * std::sin(elevation_);
    float z = distance_ * std::cos(azimuth_) * std::cos(elevation_);
    
    position_ = target_ + glm::vec3(x, y, z);
}

void Camera::clampElevation() {
    const float max_elevation = glm::pi<float>() / 2.0f - 0.01f; // Avoid gimbal lock
    elevation_ = std::clamp(elevation_, -max_elevation, max_elevation);
}

} // namespace Voronoi3D