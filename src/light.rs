// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::color::Color;
use crate::geometry::Point;

pub struct PointLight {
    pub position: Point,
    pub intensity: Color,
}

impl PointLight {
    pub fn new(position: Point, intensity: Color) -> Self {
        PointLight{position, intensity}
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn point_light_defaults() {
        let intensity = Color::white();
        let position = Point::new(0., 0., 0.);
        let light = PointLight::new(position, intensity);
        assert_eq!(light.position, position);
        assert_eq!(light.intensity, intensity);
    }
}
