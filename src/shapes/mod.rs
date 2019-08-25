mod sphere;
mod plane;

use crate::geometry::{identity, Point, Transform, Vector};
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

    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let shape_type = self.get_shape_type();
        self.data.intersect(ray, shape_type)
    }

    pub fn normal_at(&self, point: &Point) -> Vector {
        let shape_type = self.get_shape_type();
        self.data.normal_at(point, shape_type)
    }

    fn get_shape_type(&self) -> &'static dyn ShapeType {
        match self.shape_type {
            ShapeTypeMarker::Sphere => &sphere::SPHERE,
            ShapeTypeMarker::Plane => &plane::PLANE,
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
        m.ambient = 1.;
        s.set_material(m);
        assert_eq!(s.get_material(), &m);
    }
}
