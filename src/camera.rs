use crate::canvas::Canvas;
use crate::geometry::{identity, Point, Transform};
use crate::rays::Ray;
use crate::world::World;

pub struct Camera {
    hsize: usize,
    vsize: usize,
    pixel_size: f64,
    half_width: f64,
    half_height: f64,
    pub transform: Transform,
}

impl Camera {
    pub fn new(hsize: usize, vsize: usize, field_of_view: f64) -> Self {
        let half_view = (field_of_view / 2.).tan();
        let aspect = (hsize as f64) / (vsize as f64);

        let (half_width, half_height) = if aspect >= 1. {
            (half_view, half_view / aspect)
        } else {
            (half_view * aspect, half_view)
        };

        let pixel_size = (half_width * 2.) / (hsize as f64);

        Camera{
            hsize, vsize,
            pixel_size, half_width, half_height,
            transform: identity(),
        }
    }

    pub fn render(&self, world: &World) -> Canvas {
        let mut canvas = Canvas::new(self.hsize, self.vsize);

        // for each row
        for y in 0..self.vsize {
            for x in 0..self.hsize {
                let ray = self.ray_for_pixel(x, y);
                let color = world.color_at(&ray);
                canvas.write_pixel(x, y, color);
            }
        }

        canvas
    }

    fn ray_for_pixel(&self, x: usize, y: usize) -> Ray {
        let xoffset = (x as f64 + 0.5) * self.pixel_size;
        let yoffset = (y as f64 + 0.5) * self.pixel_size;

        // the camera looks toward -z, so +x is to the *left*.
        let world_x = self.half_width - xoffset;
        let world_y = self.half_height - yoffset;

        // the canvas is at z=-1
        let pixel = self.transform.inverse() * Point::new(world_x, world_y, -1.);
        let origin = self.transform.inverse() * Point::new(0., 0., 0.);
        let direction = (pixel - origin).normalize();
        Ray::new(origin, direction)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::color::Color;
    use crate::geometry::*;
    use crate::world::default_world;

    #[test]
    fn constructing_a_camera() {
        let hsize = 160;
        let vsize = 120;
        let field_of_view = std::f64::consts::PI/2.;
        let c = Camera::new(hsize, vsize, field_of_view);
        assert_eq!(c.hsize, 160);
        assert_eq!(c.vsize, 120);
        assert_eq!(c.transform, identity());
    }

    #[test]
    fn pixel_size_for_horizontal_canvas() {
        let field_of_view = std::f64::consts::PI/2.;
        let c = Camera::new(200, 125, field_of_view);
        assert_relative_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn pixel_size_for_vertical_canvas() {
        let field_of_view = std::f64::consts::PI/2.;
        let c = Camera::new(125, 200, field_of_view);
        assert_relative_eq!(c.pixel_size, 0.01);
    }

    #[test]
    fn constructing_a_ray_through_the_center_of_the_canvas() {
        let field_of_view = std::f64::consts::PI/2.;
        let c = Camera::new(201, 101, field_of_view);
        let r = c.ray_for_pixel(100, 50);
        assert_relative_eq!(r.origin, Point::new(0., 0., 0.));
        assert_relative_eq!(r.direction, Vector::new(0., 0., -1.));
    }

    #[test]
    fn constructing_a_ray_through_a_corner_of_the_canvas() {
        let field_of_view = std::f64::consts::PI/2.;
        let c = Camera::new(201, 101, field_of_view);
        let r = c.ray_for_pixel(0, 0);
        assert_relative_eq!(r.origin, Point::new(0., 0., 0.));
        assert_relative_eq!(r.direction,
                            Vector::new(0.66519, 0.33259, -0.66851),
                            max_relative = 1e-4);
    }

    #[test]
    fn constructing_a_ray_when_the_camera_is_transformed() {
        let pi = std::f64::consts::PI;
        let field_of_view = pi/2.;
        let mut c = Camera::new(201, 101, field_of_view);
        c.transform = rotation_y(pi/4.)
                    * translation(&Vector::new(0., -2., 5.));
        let r = c.ray_for_pixel(100, 50);
        assert_relative_eq!(r.origin, Point::new(0., 2., -5.));
        let x = 1./f64::sqrt(2.);
        assert_relative_eq!(r.direction, Vector::new(x, 0., -x),
                            max_relative = 1e-12);
    }

    #[test]
    fn rendering_a_world_with_a_camera() {
        let w = default_world();
        let field_of_view = std::f64::consts::PI/2.;
        let mut c = Camera::new(11, 11, field_of_view);
        let from = Point::new(0., 0., -5.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        c.transform = view_transform(&from, &to, &up);

        let image = c.render(&w);
        assert_relative_eq!(image.pixel_at(5, 5),
                            Color::new(0.38066, 0.47583, 0.2855),
                            max_relative = 1e-4);
    }
}
