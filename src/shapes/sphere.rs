use crate::geometry::{Point, Vector};
use crate::rays::Ray;

use super::bounds::Bounds;
use super::ShapeType;

pub struct Sphere {}

pub const SPHERE: Sphere = Sphere{};

impl ShapeType for Sphere {
    fn local_intersect(&self, ray: &Ray) -> Vec<f64> {
        let sphere_to_ray: Vector = ray.origin - Point::new(0., 0., 0.);

        let a: f64 = ray.direction.dot(&ray.direction);
        let b: f64 = 2. * ray.direction.dot(&sphere_to_ray);
        let c: f64 = sphere_to_ray.dot(&sphere_to_ray) - 1.;

        let discriminant = b*b - 4. * a * c;

        if discriminant < 0. {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);
            vec![t1, t2]
        }
    }

    fn local_normal_at(&self, point: &Point) -> Vector {
        point - Point::new(0., 0., 0.)
    }

    fn local_bounds(&self) -> Bounds {
        Bounds::new(Point::new(-1., -1., -1.), Point::new(1., 1., 1.))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_intersecting_sphere_at_two_points() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = SPHERE.local_intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], 4.0);
        assert_relative_eq!(xs[1], 6.0);
    }

    #[test]
    fn ray_intersecting_sphere_at_a_tangent() {
        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let xs = SPHERE.local_intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], 5.0);
        assert_relative_eq!(xs[1], 5.0);
    }

    #[test]
    fn ray_missing_sphere() {
        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let xs = SPHERE.local_intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_originating_inside_sphere() {
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = SPHERE.local_intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], -1.0);
        assert_relative_eq!(xs[1],  1.0);
    }

    #[test]
    fn ray_originating_behind_sphere() {
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let xs = SPHERE.local_intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], -6.0);
        assert_relative_eq!(xs[1], -4.0);
    }

    #[test]
    fn normal_of_sphere_on_x_axis() {
        let n = SPHERE.local_normal_at(&Point::new(1., 0., 0.));
        assert_relative_eq!(n, Vector::new(1., 0., 0.));
    }

    #[test]
    fn normal_of_sphere_on_y_axis() {
        let n = SPHERE.local_normal_at(&Point::new(0., 1., 0.));
        assert_relative_eq!(n, Vector::new(0., 1., 0.));
    }

    #[test]
    fn normal_of_sphere_on_z_axis() {
        let n = SPHERE.local_normal_at(&Point::new(0., 0., 1.));
        assert_relative_eq!(n, Vector::new(0., 0., 1.));
    }

    #[test]
    fn normal_of_sphere_at_nonaxial_point() {
        let x = 1./f64::sqrt(3.);
        let n = SPHERE.local_normal_at(&Point::new(x, x, x));
        assert_relative_eq!(n, Vector::new(x, x, x));
    }
}
