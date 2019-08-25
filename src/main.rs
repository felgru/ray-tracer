use ray_tracer::*;
use camera::Camera;
use color::Color;
use geometry::*;
use light::PointLight;
use material::Material;
use patterns::Pattern;
use shapes::Shape;
use world::World;

fn main() {
    let world = define_world();
    let camera = define_camera();

    let canvas = camera.render(&world);
    canvas.write_image("test.png").unwrap();
}

fn define_camera() -> Camera {
    let pi = std::f64::consts::PI;
    let mut camera = Camera::new(400, 200, pi/3.);
    camera.transform = view_transform(&Point::new(0., 1.5, -5.),
                                      &Point::new(0., 1., 0.),
                                      &Vector::new(0., 1., 0.));
    camera
}

fn define_world() -> World {
    let mut world = World::new();

    let mut floor = floor();
    let mut floor_material = floor_material();
    floor_material.shader.reflective = 0.5;
    floor.set_material(floor_material);
    world.add_object(floor);
    world.add_object(left_wall());
    world.add_object(right_wall());

    world.add_object(middle_sphere());
    world.add_object(right_sphere());
    world.add_object(left_sphere());

    let light_position = Point::new(-10., 10., -10.);
    let light_color    = Color::white();
    let light          = PointLight::new(light_position, light_color);

    world.add_light(light);
    world
}

fn floor_material() -> Material {
    let mut m = Material::new();
    m.pattern = Pattern::stripes(Color::new(1., 0.9, 0.9), Color::new(0.8, 0.4, 0.4));
    m.shader.specular = 0.;
    m
}

fn floor() -> Shape {
    let mut floor = Shape::plane();
    floor.set_material(floor_material());
    floor
}

fn left_wall() -> Shape {
    let mut wall = Shape::plane();
    let pi = std::f64::consts::PI;
    wall.set_transform(translation(&Vector::new(0., 0., 5.))
                       * rotation_y(-pi/4.) * rotation_x(pi/2.));
    wall.set_material(floor_material());
    wall
}

fn right_wall() -> Shape {
    let mut wall = Shape::plane();
    let pi = std::f64::consts::PI;
    wall.set_transform(translation(&Vector::new(0., 0., 5.))
                       * rotation_y(pi/4.) * rotation_x(pi/2.));
    wall.set_material(floor_material());
    wall
}

fn middle_sphere() -> Shape {
    let mut s = Shape::sphere();
    s.set_transform(translation(&Vector::new(-0.5, 1., 0.5)));
    s.set_material(sphere_material(Color::new(0.1, 1., 0.5)));
    s
}

fn right_sphere() -> Shape {
    let mut s = Shape::sphere();
    s.set_transform(translation(&Vector::new(1.5, 0.5, -0.5))
                    * scaling(0.5, 0.5, 0.5));
    s.set_material(sphere_material(Color::new(0.5, 1., 0.1)));
    s
}

fn left_sphere() -> Shape {
    let mut s = Shape::sphere();
    s.set_transform(translation(&Vector::new(-1.5, 0.33, -0.75))
                    * scaling(0.33, 0.33, 0.33));
    s.set_material(sphere_material(Color::new(1., 0.8, 0.1)));
    s
}

fn sphere_material(color: Color) -> Material {
    let mut m = Material::new();
    m.pattern = Pattern::stripes(color, color * 0.5);
    m.shader.diffuse = 0.7;
    m.shader.specular = 0.3;
    m
}
