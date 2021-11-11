// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::geometry::{Point, Vector};
use crate::rays::Ray;

use super::bounds::Bounds;
use super::ShapeType;

pub struct Plane {}

pub const PLANE: Plane = Plane{};

impl ShapeType for Plane {
    fn local_intersect(&self, ray: &Ray) -> Vec<f64> {
        if ray.direction[1].abs() < 1e-12 {
            return Vec::new();
        }

        let t = -ray.origin[1] / ray.direction[1];
        vec![t]
    }

    fn local_normal_at(&self, _point: &Point) -> Vector {
        Vector::new(0., 1., 0.)
    }

    fn local_bounds(&self) -> Bounds {
        let infty = std::f64::INFINITY;
        Bounds::new(Point::new(-infty, 0., -infty),
                    Point::new(infty, 0., infty))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn normal_of_a_plane_is_invariant_of_point() {
        let n1 = PLANE.local_normal_at(&Point::new(0., 0., 0.));
        let n2 = PLANE.local_normal_at(&Point::new(10., 0., -10.));
        let n3 = PLANE.local_normal_at(&Point::new(-5., 0., 150.));
        assert_eq!(n1, Vector::new(0., 1., 0.));
        assert_eq!(n2, Vector::new(0., 1., 0.));
        assert_eq!(n3, Vector::new(0., 1., 0.));
    }

    #[test]
    fn intersect_with_parrallel_ray() {
        let r = Ray::new(Point::new(0., 10., 0.), Vector::new(0., 0., 1.));
        let xs = PLANE.local_intersect(&r);
        assert!(xs.is_empty());
    }

    #[test]
    fn intersect_with_coplanar_ray() {
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = PLANE.local_intersect(&r);
        assert!(xs.is_empty());
    }

    #[test]
    fn ray_intersecting_from_above() {
        let r = Ray::new(Point::new(0., 1., 0.), Vector::new(0., -1., 0.));
        let xs = PLANE.local_intersect(&r);
        assert_eq!(xs.len(), 1);
        assert_relative_eq!(xs[0], 1.);
    }

    #[test]
    fn ray_intersecting_from_below() {
        let r = Ray::new(Point::new(0., -1., 0.), Vector::new(0., 1., 0.));
        let xs = PLANE.local_intersect(&r);
        assert_eq!(xs.len(), 1);
        assert_relative_eq!(xs[0], 1.);
    }
}
