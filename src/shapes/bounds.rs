// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use std::mem;

use approx::{AbsDiffEq, RelativeEq};

use crate::geometry::{Point, Transform};
use crate::rays::Ray;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Bounds {
    min: Point,
    max: Point,
}

impl Bounds {
    pub fn new(min: Point, max: Point) -> Self {
        Bounds{min, max}
    }

    pub fn empty() -> Self {
        let nan = std::f64::NAN;
        Bounds::new(Point::new(nan, nan, nan),
                    Point::new(nan, nan, nan))
    }

    pub fn from_points(points: &[Point]) -> Self {
        let mut min = points[0];
        let mut max = points[1];

        for p in points[1..].iter() {
            for i in 0..3 {
                min[i] = min[i].min(p[i]);
                max[i] = max[i].max(p[i]);
            }
        }
        Bounds::new(min, max)
    }

    pub fn intersects(&self, ray: &Ray) -> bool {
        let (xtmin, xtmax) = check_axis(ray.origin[0], ray.direction[0],
                                        self.min[0], self.max[0]);
        let (ytmin, ytmax) = check_axis(ray.origin[1], ray.direction[1],
                                        self.min[1], self.max[1]);
        let (ztmin, ztmax) = check_axis(ray.origin[2], ray.direction[2],
                                        self.min[2], self.max[2]);

        let tmin = xtmin.max(ytmin).max(ztmin);
        let tmax = xtmax.min(ytmax).min(ztmax);

        tmin <= tmax
    }

    pub fn merge_bounds(&mut self, other: &Bounds) {
        for i in 0..3 {
            self.min[i] = self.min[i].min(other.min[i]);
            self.max[i] = self.max[i].max(other.max[i]);
        }
    }

    pub fn to_global(self, t: &Transform) -> Self {
        let mut corners = Vec::with_capacity(8);
        for &x in &[self.min[0], self.max[0]] {
            for &y in &[self.min[1], self.max[1]] {
                for &z in &[self.min[2], self.max[2]] {
                    let p = t * Point::new(x, y, z);
                    corners.push(p);
                }
            }
        }
        Self::from_points(&corners)
    }
}

fn check_axis(origin: f64, direction: f64, min: f64, max: f64) -> (f64, f64) {
    let tmin_numerator = min - origin;
    let tmax_numerator = max - origin;

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

impl AbsDiffEq for Bounds {
    type Epsilon = <f64 as AbsDiffEq>::Epsilon;

    fn default_epsilon() -> Self::Epsilon {
        f64::default_epsilon()
    }

    fn abs_diff_eq(&self, other: &Self, epsilon: Self::Epsilon) -> bool {
        Point::abs_diff_eq(&self.min, &other.min, epsilon) &&
        Point::abs_diff_eq(&self.max, &other.max, epsilon)
    }
}

impl RelativeEq for Bounds {
    fn default_max_relative() -> Self::Epsilon {
        f64::default_max_relative()
    }

    fn relative_eq(&self, other: &Self, epsilon: Self::Epsilon,
                   max_relative: Self::Epsilon) -> bool {
        Point::relative_eq(&self.min, &other.min, epsilon, max_relative) &&
        Point::relative_eq(&self.max, &other.max, epsilon, max_relative)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;

    #[test]
    fn empty_bounds_never_intersect() {
        let bounds = Bounds::empty();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        assert!(!bounds.intersects(&r));
    }

    #[test]
    fn adding_bounds_to_empty_bounds() {
        let mut bounds1 = Bounds::empty();
        let bounds2 = Bounds::new(Point::new(-1., -1., -1.),
                                  Point::new(1., 1., 1.));
        bounds1.merge_bounds(&bounds2);
        assert_relative_eq!(bounds1, bounds2);
    }
}
