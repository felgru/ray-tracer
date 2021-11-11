// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::geometry::{Point, Transform, Vector};

pub struct Ray {
    pub origin: Point,
    pub direction: Vector,
}

impl Ray {
    pub fn new(origin: Point, direction: Vector) -> Ray {
        Ray{origin, direction}
    }

    pub fn position(&self, t: f64) -> Point {
        self.origin + t * self.direction
    }

    pub fn transform(&self, m: &Transform) -> Ray {
        Ray::new(m * self.origin, m * self.direction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn creating_and_querying_a_ray() {
        let origin = Point::new(1., 2., 3.);
        let direction = Vector::new(4., 5., 6.);
        let r = Ray::new(origin, direction);
        assert_eq!(r.origin, origin);
        assert_eq!(r.direction, direction);
    }

    #[test]
    fn position_after_given_time_of_flight() {
        let r = Ray::new(Point::new(2., 3., 4.), Vector::new(1., 0., 0.));
        assert_relative_eq!(r.position(0.), Point::new(2., 3., 4.));
        assert_relative_eq!(r.position(1.), Point::new(3., 3., 4.));
        assert_relative_eq!(r.position(-1.), Point::new(1., 3., 4.));
        assert_relative_eq!(r.position(2.5), Point::new(4.5, 3., 4.));
    }

    #[test]
    fn translating_a_ray() {
        let r = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        use crate::geometry::translation;
        let m = translation(&Vector::new(3., 4., 5.));
        let r2 = r.transform(&m);
        assert_relative_eq!(r2.origin, Point::new(4., 6., 8.));
        assert_relative_eq!(r2.direction, Vector::new(0., 1., 0.));
    }

    #[test]
    fn scaling_a_ray() {
        let r = Ray::new(Point::new(1., 2., 3.), Vector::new(0., 1., 0.));
        use crate::geometry::scaling;
        let m = scaling(2., 3., 4.);
        let r2 = r.transform(&m);
        assert_relative_eq!(r2.origin, Point::new(2., 6., 12.));
        assert_relative_eq!(r2.direction, Vector::new(0., 3., 0.));
    }
}
