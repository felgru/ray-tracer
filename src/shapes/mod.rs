mod cube;
mod sphere;
mod plane;

use crate::color::Color;
use crate::geometry::{identity, Point, Transform, Vector};
use crate::light::PointLight;
use crate::material::Material;
use crate::rays::Ray;

pub struct Shape {
    pub data: ShapeData,
    pub shape_type: ShapeTypeMarker,
}

impl Shape {
    pub fn new(shape_type: ShapeTypeMarker) -> Self {
        Shape{data: ShapeData::new(), shape_type}
    }

    pub fn sphere() -> Self {
        Self::new(ShapeTypeMarker::Sphere)
    }

    pub fn plane() -> Self {
        Self::new(ShapeTypeMarker::Plane)
    }

    pub fn cube() -> Self {
        Self::new(ShapeTypeMarker::Cube)
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let shape_type = self.get_shape_type();
        self.data.intersect(ray, shape_type)
    }

    pub fn normal_at(&self, point: &Point) -> Vector {
        let shape_type = self.get_shape_type();
        self.data.normal_at(point, shape_type)
    }

    pub fn color_at(&self, point: &Point) -> Color {
        let pattern = self.get_material().pattern;
        let object_point  = self.get_transform().inverse() * point;
        let pattern_point = pattern.transform.inverse() * object_point;
        pattern.local_color_at(&pattern_point)
    }

    pub fn lighting(&self, light: &PointLight,
                    position: &Point, eyev: &Vector, normalv: &Vector,
                    in_shadow: bool) -> Color {
        let color = self.color_at(position);
        let shader = self.get_material().shader;
        shader.lighting(color, light, position, eyev, normalv, in_shadow)
    }

    fn get_shape_type(&self) -> &'static dyn ShapeType {
        match self.shape_type {
            ShapeTypeMarker::Sphere => &sphere::SPHERE,
            ShapeTypeMarker::Plane => &plane::PLANE,
            ShapeTypeMarker::Cube => &cube::CUBE,
        }
    }

    pub fn set_material(&mut self, material: Material) {
        self.data.material = material;
    }

    pub fn get_material(&self) -> &Material {
        &self.data.material
    }

    pub fn set_transform(&mut self, transform: Transform) {
        self.data.transform = transform;
    }

    pub fn get_transform(&self) -> &Transform {
        &self.data.transform
    }
}

pub struct ShapeData {
    pub transform: Transform,
    pub material: Material,
}

#[derive(Copy, Clone, PartialEq, Eq)]
pub enum ShapeTypeMarker {
    Sphere,
    Plane,
    Cube,
}

pub trait ShapeType {
    fn local_intersect(&self, ray: &Ray) -> Vec<f64>;
    fn local_normal_at(&self, point: &Point) -> Vector;
}

impl ShapeData {
    pub fn new() -> Self {
        ShapeData{transform: identity(), material: Material::new()}
    }

    pub fn intersect(&self, ray: &Ray, shape_type: &dyn ShapeType) -> Vec<f64> {
        let local_ray = ray.transform(&self.transform.inverse());
        shape_type.local_intersect(&local_ray)
    }

    pub fn normal_at(&self, point: &Point, shape_type: &dyn ShapeType) -> Vector {
        let inverse_transform = self.transform.inverse();
        let local_point = inverse_transform * point;
        let local_normal = shape_type.local_normal_at(&local_point);
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        (inv_transp * local_normal).normalize()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::patterns::Pattern;

    #[test]
    fn default_transformation_of_a_sphere() {
        let s = Shape::sphere();
        assert_relative_eq!(s.get_transform(), &identity());
    }

    #[test]
    fn changing_transformation_of_a_sphere() {
        let mut s = Shape::sphere();
        use crate::geometry::translation;
        let t = translation(&Vector::new(2., 3., 4.));
        s.set_transform(t);
        assert_relative_eq!(s.get_transform(), &t);
    }

    #[test]
    fn intersecting_scaled_sphere_with_a_ray() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Shape::sphere();
        use crate::geometry::scaling;
        s.set_transform(scaling(2., 2., 2.));
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], 3.0);
        assert_relative_eq!(xs[1], 7.0);
    }

    #[test]
    fn intersecting_translated_sphere_with_a_ray() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Shape::sphere();
        use crate::geometry::translation;
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn normal_of_translated_sphere() {
        let mut s = Shape::sphere();
        use crate::geometry::translation;
        s.set_transform(translation(&Vector::new(0., 1., 0.)));
        let n = s.normal_at(&Point::new(0., 1.70711, -0.70711));
        assert_relative_eq!(n, Vector::new(0., 0.70711, -0.70711), max_relative = 1e-5);
    }

    #[test]
    fn normal_of_transformed_sphere() {
        let mut s = Shape::sphere();
        use crate::geometry::{scaling, rotation_z};
        s.set_transform(scaling(1., 0.5, 1.) * rotation_z(std::f64::consts::PI/5.));
        let x = 1./f64::sqrt(2.);
        let n = s.normal_at(&Point::new(0., x, -x));
        assert_relative_eq!(n, Vector::new(0., 0.97014, -0.24254), max_relative = 1e-4);
    }

    #[test]
    fn default_sphere_material() {
        let s = Shape::sphere();
        assert_eq!(s.get_material(), &Material::new());
    }

    #[test]
    fn assign_material_to_sphere() {
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.shader.ambient = 1.;
        s.set_material(m);
        assert_eq!(s.get_material(), &m);
    }

    #[test]
    fn stripes_with_transformed_sphere() {
        let mut object = Shape::sphere();
        object.set_transform(scaling(2., 2., 2.));
        let pattern = Pattern::stripes(Color::white(), Color::black());
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern() {
        let mut object = Shape::sphere();
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = scaling(2., 2., 2.);
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern_and_object() {
        let mut object = Shape::sphere();
        object.set_transform(scaling(2., 2., 2.));
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = translation(&Vector::new(0.5, 0., 0.));
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(2.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    fn set_pattern_of_object(object: &mut Shape, pattern: Pattern) {
        let mut material = object.get_material().clone();
        material.pattern = pattern;
        object.set_material(material);
    }

    #[test]
    fn lighting_with_stripe_pattern() {
        let mut m = Material::new();
        m.pattern = Pattern::stripes(Color::white(), Color::black());
        m.shader.ambient = 1.;
        m.shader.diffuse = 0.;
        m.shader.specular = 0.;
        let mut obj = Shape::sphere();
        obj.set_material(m);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let c1 = obj.lighting(&light, &Point::new(0.9, 0., 0.),
                              &eyev, &normalv, false);
        let c2 = obj.lighting(&light, &Point::new(1.1, 0., 0.),
                              &eyev, &normalv, false);
        assert_relative_eq!(c1, Color::white());
        assert_relative_eq!(c2, Color::black());
    }
}
