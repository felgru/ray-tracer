use std::mem;

use crate::geometry::{Point, Vector};
use crate::rays::Ray;

use super::ShapeType;

pub struct Cube {}

pub const CUBE: Cube = Cube{};

impl ShapeType for Cube {
    fn local_intersect(&self, ray: &Ray) -> Vec<f64> {
        let (xtmin, xtmax) = check_axis(ray.origin[0], ray.direction[0]);
        let (ytmin, ytmax) = check_axis(ray.origin[1], ray.direction[1]);
        let (ztmin, ztmax) = check_axis(ray.origin[2], ray.direction[2]);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        if tmin <= tmax {
            vec![tmin, tmax]
        } else {
            Vec::new()
        }
    }

    fn local_normal_at(&self, point: &Point) -> Vector {
        let maxc = point[0].abs().max(point[1].abs()).max(point[2].abs());

        if maxc == point[0].abs() {
            Vector::new(point[0], 0., 0.)
        } else if maxc == point[1].abs() {
            Vector::new(0., point[1], 0.)
        } else {
            Vector::new(0., 0., point[2])
        }
    }
}

fn check_axis(origin: f64, direction: f64) -> (f64, f64) {
    let tmin_numerator = -1. - origin;
    let tmax_numerator = 1. - origin;

    let (mut tmin, mut tmax);
    if direction.abs() >= 1e-12 {
        tmin = tmin_numerator / direction;
        tmax = tmax_numerator / direction;
    } else {
        tmin = tmin_numerator * std::f64::INFINITY;
        tmax = tmax_numerator * std::f64::INFINITY;
    }

    if tmin > tmax {
        mem::swap(&mut tmin, &mut tmax);
    }
    (tmin, tmax)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn ray_intersects_a_cube() {
        // +x
        ray_intersection_test(Point::new(5., 0.5, 0.),
                              Vector::new(-1., 0., 0.), 4., 6.);
        // -x
        ray_intersection_test(Point::new(-5., 0.5, 0.),
                              Vector::new(1., 0., 0.), 4., 6.);
        // +y
        ray_intersection_test(Point::new(0.5, 5., 0.),
                              Vector::new(0., -1., 0.), 4., 6.);
        // -y
        ray_intersection_test(Point::new(0.5, -5., 0.),
                              Vector::new(0., 1., 0.), 4., 6.);
        // +z
        ray_intersection_test(Point::new(0.5, 0., 5.),
                              Vector::new(0., 0., -1.), 4., 6.);
        // -z
        ray_intersection_test(Point::new(0.5, 0., -5.),
                              Vector::new(0., 0., 1.), 4., 6.);
        // inside
        ray_intersection_test(Point::new(0., 0.5, 0.),
                              Vector::new(0., 0., 1.), -1., 1.);
    }

    fn ray_intersection_test(origin: Point, direction: Vector,
                             t1: f64, t2: f64) {
        let ray = Ray::new(origin, direction);
        let xs = CUBE.local_intersect(&ray);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0], t1);
        assert_relative_eq!(xs[1], t2);
    }

    #[test]
    fn ray_misses_a_cube() {
        ray_miss_test(Point::new(-2., 0., 0.), Vector::new(0.2673, 0.5345, 0.8018));
        ray_miss_test(Point::new(0., -2., 0.), Vector::new(0.8018, 0.2673, 0.5345));
        ray_miss_test(Point::new(0., 0., -2.), Vector::new(0.5345, 0.8018, 0.2673));
        ray_miss_test(Point::new(2., 0., 2.), Vector::new(0., 0., -1.));
        ray_miss_test(Point::new(0., 2., 2.), Vector::new(0., -1., 0.));
        ray_miss_test(Point::new(2., 2., 0.), Vector::new(-1., 0., 0.));
    }

    fn ray_miss_test(origin: Point, direction: Vector) {
        let ray = Ray::new(origin, direction);
        let xs = CUBE.local_intersect(&ray);
        assert!(xs.is_empty());
    }

    #[test]
    fn normal_on_the_surface_of_a_cube() {
        normal_test(Point::new(1., 0.5, -0.8), Vector::new(1., 0., 0.));
        normal_test(Point::new(-1., -0.2, 0.9), Vector::new(-1., 0., 0.));
        normal_test(Point::new(-0.4, 1., -0.1), Vector::new(0., 1., 0.));
        normal_test(Point::new(0.3, -1., -0.7), Vector::new(0., -1., 0.));
        normal_test(Point::new(-0.6, 0.3, 1.), Vector::new(0., 0., 1.));
        normal_test(Point::new(0.4, 0.4, -1.), Vector::new(0., 0., -1.));
        normal_test(Point::new(1., 1., 1.), Vector::new(1., 0., 0.));
        normal_test(Point::new(-1., -1., -1.), Vector::new(-1., 0., 0.));
    }

    fn normal_test(point: Point, normal: Vector) {
        assert_relative_eq!(CUBE.local_normal_at(&point), normal);
    }
}
