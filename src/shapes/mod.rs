pub mod bounds;
mod cube;
mod sphere;
mod plane;

use bounds::Bounds;
use crate::geometry::{Point, Vector};
use crate::material::Material;
use crate::rays::Ray;

#[derive(Debug, PartialEq, Eq, PartialOrd, Ord, Copy, Clone)]
pub struct ShapeIndex {
    pub index: usize,
}

impl ShapeIndex {
    pub fn new(index: usize) -> Self {
        ShapeIndex{index}
    }
}

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

    pub fn local_intersect(&self, ray: &Ray) -> Vec<f64> {
        let shape_type = self.get_shape_type();
        shape_type.local_intersect(ray)
    }

    pub fn local_normal_at(&self, local_point: &Point) -> Vector {
        let shape_type = self.get_shape_type();
        shape_type.local_normal_at(&local_point)
    }

    pub fn local_bounds(&self) -> Bounds {
        let shape_type = self.get_shape_type();
        shape_type.local_bounds()
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
}

pub struct ShapeData {
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
    fn local_bounds(&self) -> Bounds;
}

impl ShapeData {
    pub fn new() -> Self {
        ShapeData{material: Material::new()}
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::object_store::ObjectStore;

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
    fn creating_a_transformed_sphere() {
        use crate::geometry::translation;
        let t = translation(&Vector::new(2., 3., 4.));
        let mut objects = ObjectStore::new();
        let s = objects.add_shape(Shape::sphere(), t);
        assert_relative_eq!(objects.get_transform_of_object(s.into()), &t);
    }

    #[test]
    fn intersecting_scaled_sphere_with_a_ray() {
        use crate::geometry::scaling;
        let mut objects = ObjectStore::new();
        let s = objects.add_shape(Shape::sphere(), scaling(2., 2., 2.));
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = objects.intersect(s.into(), &r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0].t, 3.0);
        assert_relative_eq!(xs[1].t, 7.0);
    }

    #[test]
    fn intersecting_translated_sphere_with_a_ray() {
        use crate::geometry::translation;
        let mut objects = ObjectStore::new();
        let s = objects.add_shape(Shape::sphere(),
                                  translation(&Vector::new(5., 0., 0.)));
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = objects.intersect(s.into(), &r);
        assert!(xs.is_empty());
    }
}
