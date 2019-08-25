use crate::color::Color;
use crate::geometry::{Point, reflect, Vector};
use crate::light::PointLight;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Material {
    pub color: Color,
    pub ambient: f64,
    pub diffuse: f64,
    pub specular: f64,
    pub shininess: f64,
}

impl Material {
    pub fn new() -> Self {
        Material {
            color: Color::new(1., 1., 1.),
            ambient: 0.1,
            diffuse: 0.9,
            specular: 0.9,
            shininess: 200.,
        }
    }

    /// Compute light intensity using Phong shader.
    pub fn lighting(&self, light: &PointLight, position: &Point,
                    eyev: &Vector, normalv: &Vector,
                    in_shadow: bool) -> Color {
        let effective_color = self.color * light.intensity;
        let lightv = (light.position - position).normalize();
        let ambient = effective_color * self.ambient;

        let diffuse;
        let specular;
        let light_dot_normal = lightv.dot(normalv);
        if in_shadow || light_dot_normal < 0. {
            diffuse = Color::black();
            specular = Color::black();
        } else {
            diffuse = effective_color * self.diffuse * light_dot_normal;

            let reflectv = reflect(&-lightv, &normalv);
            let reflect_dot_eye = reflectv.dot(eyev);
            if reflect_dot_eye <= 0. {
                specular = Color::black();
            } else {
                let factor = reflect_dot_eye.powf(self.shininess);
                specular = light.intensity * self.specular * factor;
            }
        }

        ambient + diffuse + specular
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_material() {
        let m = Material::new();
        assert_eq!(m.color, Color::new(1., 1., 1.));
        assert_eq!(m.ambient, 0.1);
        assert_eq!(m.diffuse, 0.9);
        assert_eq!(m.specular, 0.9);
        assert_eq!(m.shininess, 200.);
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
        let in_shadow = false;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 1.9;
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity));
    }

    #[test]
    fn lighting_with_eye_between_light_and_surface_eye_offset_45_degree() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let x = 1./f64::sqrt(2.);
        let eyev = Vector::new(0., x, -x);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
        let in_shadow = false;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 1.0;
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity));
    }

    #[test]
    fn lighting_with_eye_opposite_surface_light_offset_45_degree() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::new(1., 1., 1.));
        let in_shadow = false;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 0.1 + 0.9 * 1./f64::sqrt(2.) + 0.; // approximately 0.7364
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity));
    }

    #[test]
    fn lighting_with_eye_in_path_of_reflection_vector() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let x = 1./f64::sqrt(2.);
        let eyev = Vector::new(0., -x, -x);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 10., -10.), Color::new(1., 1., 1.));
        let in_shadow = false;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 0.1 + 0.9 * x + 0.9; // approximately 1.6364
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity), max_relative = 1e-12);
    }

    #[test]
    fn lighting_with_light_behind_surface() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., 10.), Color::new(1., 1., 1.));
        let in_shadow = false;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 0.1;
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity));
    }

    #[test]
    fn lighting_with_surface_in_shadow() {
        let m = Material::new();
        let position = Point::new(0., 0., 0.);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::new(1., 1., 1.));
        let in_shadow = true;
        let result = m.lighting(&light, &position, &eyev, &normalv, in_shadow);
        let intensity = 0.1;
        assert_relative_eq!(result, Color::new(intensity, intensity, intensity));
    }
}
