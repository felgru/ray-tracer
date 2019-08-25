use ray_tracer::*;
use canvas::Canvas;
use color::Color;
use geometry::Point;
use light::PointLight;
use rays::Ray;
use shapes::sphere::Sphere;
use world::World;

fn main() {
    let ray_origin = Point::new(0., 0., -5.);
    let wall_z = 10.;
    let wall_size = 7.0;
    let canvas_pixels = 200;
    let pixel_size = wall_size / (canvas_pixels as f64);
    let half = wall_size / 2.;

    let mut canvas = Canvas::new(canvas_pixels, canvas_pixels);

    let world = define_world();

    // for each row
    for y in 0..canvas_pixels {
        // compute the world y coordinate (top = +half, bottom = -half)
        let world_y = half - pixel_size * (y as f64);

        // for each pixel in the row
        for x in 0..canvas_pixels {
            // compute the world x coordinate (left = -half, right = half)
            let world_x = -half + pixel_size * (x as f64);

            // describe the point on the wall that the ray will target
            let position = Point::new(world_x, world_y, wall_z);

            let ray = Ray::new(ray_origin, (position - ray_origin).normalize());
            let color = world.color_at(&ray);
            canvas.write_pixel(x, y, color);
        }
    }

    canvas.write_image("test.png").unwrap();
}

fn define_world() -> World {
    let mut sphere = Sphere::new();
    sphere.material.color = Color::new(1., 0.2, 1.);

    let light_position = Point::new(-10., 10., -10.);
    let light_color    = Color::new(1., 1., 1.);
    let light          = PointLight::new(light_position, light_color);

    let mut world = World::new();
    world.add_object(sphere);
    world.add_light(light);
    world
}
