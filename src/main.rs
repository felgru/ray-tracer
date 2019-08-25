use ray_tracer::*;
use camera::Camera;
use color::Color;
use geometry::{Point, Vector, view_transform};
use light::PointLight;
use shapes::sphere::Sphere;
use world::World;

fn main() {
    let field_of_view = std::f64::consts::PI/2.;
    let canvas_pixels = 200;
    let mut c = Camera::new(canvas_pixels, canvas_pixels, field_of_view);
    let from = Point::new(0., 0., -5.);
    let to = Point::new(0., 0., 0.);
    let up = Vector::new(0., 1., 0.);
    c.transform = view_transform(&from, &to, &up);

    let world = define_world();

    let canvas = c.render(&world);
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
