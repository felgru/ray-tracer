use std::ops::Add;
use std::ops::Sub;
use std::ops::Mul;

use crate::float::Approx;

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Color {
    pub red: f64,
    pub green: f64,
    pub blue: f64,
}

impl Color {
    pub fn new(red: f64, green: f64, blue: f64) -> Self {
        Color{red, green, blue}
    }

    pub fn black() -> Self {
        Color::new(0., 0., 0.)
    }

    pub fn clamp(self) -> Self {
        Color {
            red: Self::clamp_component(self.red),
            green: Self::clamp_component(self.green),
            blue: Self::clamp_component(self.blue),
        }
    }

    fn clamp_component(comp: f64) -> f64 {
        let mut x = comp;
        if x < 0. { x = 0.; }
        if x > 1. { x = 1.; }
        x
    }
}

impl Default for Color {
    fn default() -> Self {
        Color::new(0., 0., 0.)
    }
}

impl Add for Color {
    type Output = Self;

    fn add(self, rhs: Self) -> Self {
        Color {
            red: self.red + rhs.red,
            green: self.green + rhs.green,
            blue: self.blue + rhs.blue,
        }
    }
}

impl Sub for Color {
    type Output = Self;

    fn sub(self, rhs: Self) -> Self {
        Color {
            red: self.red - rhs.red,
            green: self.green - rhs.green,
            blue: self.blue - rhs.blue,
        }
    }
}

impl Mul for Color {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self {
        Color {
            red: self.red * rhs.red,
            green: self.green * rhs.green,
            blue: self.blue * rhs.blue,
        }
    }
}

impl Mul<f64> for Color {
    type Output = Self;

    fn mul(self, rhs: f64) -> Self {
        Color {
            red: self.red * rhs,
            green: self.green * rhs,
            blue: self.blue * rhs,
        }
    }
}

impl Approx for Color {
    fn approx(self, rhs: Self, eps: f64) -> bool {
        self.red.approx(rhs.red, eps)
            && self.green.approx(rhs.green, eps)
            && self.blue.approx(rhs.blue, eps)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_approx;

    const EPS: f64 = 1e-8;

    #[test]
    fn color_values() {
        let c = Color::new(-0.5, 0.4, 1.7);
        assert_eq!(c.red, -0.5);
        assert_eq!(c.green, 0.4);
        assert_eq!(c.blue, 1.7);
    }

    #[test]
    fn default_color() {
        let c = Color::default();
        assert_eq!(c.red, 0.);
        assert_eq!(c.green, 0.);
        assert_eq!(c.blue, 0.);
    }

    #[test]
    fn add_colors() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_approx!(c1 + c2, Color::new(1.6, 0.7, 1.0), EPS);
    }

    #[test]
    fn subtract_colors() {
        let c1 = Color::new(0.9, 0.6, 0.75);
        let c2 = Color::new(0.7, 0.1, 0.25);
        assert_approx!(c1 - c2, Color::new(0.2, 0.5, 0.5), EPS);
    }

    #[test]
    fn multiply_color_by_scalar() {
        let c = Color::new(0.2, 0.3, 0.4);
        assert_eq!(c * 2., Color::new(0.4, 0.6, 0.8));
    }

    #[test]
    fn multiply_colors() {
        let c1 = Color::new(1., 0.2, 0.4);
        let c2 = Color::new(0.9, 1., 0.1);
        assert_approx!(c1 * c2, Color::new(0.9, 0.2, 0.04), EPS);
    }
}
