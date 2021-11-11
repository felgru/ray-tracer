// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use crate::color::Color;
use std::path::Path;

pub struct Canvas {
    width: usize,
    height: usize,
    data: Vec<Color>,
}

impl Canvas {
    pub fn new(width: usize, height: usize) -> Self {
        Canvas {
            width,
            height,
            data: vec![Color::default(); width*height],
        }
    }

    pub fn width(&self) -> usize {
        self.width
    }

    pub fn height(&self) -> usize {
        self.height
    }

    pub fn write_pixel(&mut self, x: usize, y: usize, c: Color) {
        self.data[x + y * self.width] = c;
    }

    pub fn pixel_at(&self, x: usize, y: usize) -> Color {
        self.data[x + y * self.width]
    }

    pub fn write_image<P: AsRef<Path>>(&self, path: &P)
                                -> Result<(), image::ImageError> {
        let img = self.to_image();
        img.save(path)
    }

    fn to_image(&self) -> image::RgbImage {
        use image::{ImageBuffer, RgbImage};
        let mut img: RgbImage = ImageBuffer::new(self.width as u32,
                                                 self.height as u32);
        for (dest, src) in img.pixels_mut().zip(self.data.iter()) {
            *dest = src.to_rgb();
        }
        img
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canvas_construction() {
        let c = Canvas::new(10, 20);
        assert_eq!(c.width(), 10);
        assert_eq!(c.height(), 20);
        for color in c.data {
            assert_eq!(color, Color::default());
        }
    }

    #[test]
    fn writing_pixel_to_canvas() {
        let mut c = Canvas::new(10, 20);
        let red = Color::new(1., 0., 0.);
        c.write_pixel(2, 3, red);
        assert_eq!(c.pixel_at(2, 3), red);
    }
}
