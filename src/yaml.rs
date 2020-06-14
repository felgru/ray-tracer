use std::collections::HashMap;
use std::ops::Index;
use yaml_rust::{Yaml, YamlLoader};

use crate::camera::Camera;
use crate::color::Color;
use crate::csg::{CsgIndex, CsgOperator};
use crate::geometry::{identity, Point, rotation_x, rotation_y, rotation_z, scaling, Transform,
                      translation, Vector, view_transform};
use crate::group::GroupIndex;
use crate::light::PointLight;
use crate::material::Material;
use crate::object_store::ObjectIndex;
use crate::patterns::Pattern;
use crate::shapes::{Shape, ShapeIndex};
use crate::world::{World, WorldBuilder};

pub fn load_world_and_cameras_from_str(s: &str)
        -> (World, Vec<Camera>) {
    // TODO: add error handling
    let docs = YamlLoader::load_from_str(s).unwrap();
    let doc = &docs[0];
    let mut world = WorldBuilder::new();
    let mut cameras = Vec::<Camera>::new();
    let mut materials = MaterialStore::new();
    let mut groups = GroupStore::new();
    for entry in doc.as_vec().unwrap() {
        if entry.as_hash().unwrap().contains_key(&Yaml::from_str("add")) {
            match entry["add"].as_str().unwrap() {
                "camera" => cameras.push(load_camera(&entry)),
                "light" => {
                    let light = load_light(&entry);
                    world.add_light(light);
                },
                "plane" => {
                    load_toplevel_shape(Shape::plane(), &entry,
                                        &mut world, &materials, &groups);
                },
                "sphere" => {
                    load_toplevel_shape(Shape::sphere(), &entry,
                                        &mut world, &materials, &groups);
                },
                "cube" => {
                    load_toplevel_shape(Shape::cube(), &entry,
                                        &mut world, &materials, &groups);
                },
                "csg" => {
                    let csg = load_csg(&entry, &mut world, &materials);
                    world.add_object_to_scene(csg.into());
                },
                unknown => { // TODO: error
                    println!("trying to add unknown object: {}", unknown);
                }
            }
        } else if entry.as_hash().unwrap()
                       .contains_key(&Yaml::from_str("define")) {
            match entry["define"].as_str().unwrap() {
                "material" => materials.parse(&entry),
                "group" => groups.parse(&entry, &mut world),
                unknown => { // TODO: error
                    println!("trying to define unknown property: {}", unknown);
                }
            }
        } else {
            // TODO: error
        }
    }
    (world.build(), cameras)
}

fn load_camera(entry: &Yaml) -> Camera {
    // TODO: add error handling
    let width = entry["width"].as_i64().unwrap() as usize;
    let height = entry["height"].as_i64().unwrap() as usize;
    let field_of_view = parse_f64_expression(&entry["field-of-view"]);
    let mut camera = Camera::new(width, height, field_of_view);
    if entry.as_hash().unwrap().contains_key(&Yaml::from_str("up")) {
        camera.transform = view_transform(&parse_point(&entry["from"]),
                                          &parse_point(&entry["to"]),
                                          &parse_vector(&entry["up"]));
    }
    camera
}

fn load_light(entry: &Yaml) -> PointLight {
    let position = parse_point(&entry["at"]);
    let color = parse_color(&entry["intensity"]);
    PointLight::new(position, color)
}

fn load_toplevel_shape<'a>(shape: Shape, entry: &'a Yaml,
                           world: &mut WorldBuilder,
                           materials: &MaterialStore,
                           groups: &GroupStore) {
    let (shape, group) = load_shape_properties(shape, entry, world, materials);
    if let Some(group) = group {
        let group = groups[group];
        world.set_group_of(shape.into(), group);
    } else {
        world.add_object_to_scene(shape.into());
    }
}

fn load_shape_properties<'a>(mut shape: Shape, entry: &'a Yaml,
                             world: &mut WorldBuilder,
                             materials: &MaterialStore)
                            -> (ShapeIndex, Option<&'a str>) {
    let mut group = None;
    let mut transform = identity();
    for (key, entry) in entry.as_hash().unwrap().iter() {
        match key.as_str().unwrap() {
            "type" => {},
            "material" => {
                let material = entry.as_str().unwrap();
                let material = materials[material];
                shape.set_material(material);
            },
            "transform" => {
                transform = parse_transform(&entry);
            },
            "group" => {
                group = entry.as_str();
            },
            "add" => {
                // ignore add: {shape_type}
            },
            unknown => {
                println!("unknown shape parameter: {}", unknown);
            }
        }
    }
    (world.add_shape(shape, transform), group)
}

fn load_object(entry: &Yaml, mut world: &mut WorldBuilder,
               materials: &MaterialStore) -> ObjectIndex {
    if !entry.as_hash().unwrap().contains_key(&Yaml::from_str("add")) {
        panic!("Expected keyword 'add' not found.");
    }
    let (obj, group) = match entry["add"].as_str().unwrap() {
        "plane" => {
            let (s, group) = load_shape_properties(Shape::plane(), &entry,
                                                   &mut world, &materials);
            (s.into(), group)
        },
        "sphere" => {
            let (s, group) = load_shape_properties(Shape::sphere(), &entry,
                                                   &mut world, &materials);
            (s.into(), group)
        },
        "cube" => {
            let (s, group) = load_shape_properties(Shape::cube(), &entry,
                                                   &mut world, &materials);
            (s.into(), group)
        },
        "csg" => {
            let csg = load_csg(&entry, &mut world, &materials).into();
            (csg, None)
        },
        "group" => unimplemented!(),
        unknown => { // TODO: error
            panic!("trying to add unknown object: {}", unknown);
        }
    };
    assert_eq!(group, None);
    obj
}

fn load_csg(entry: &Yaml, world: &mut WorldBuilder,
            materials: &MaterialStore) -> CsgIndex {
    let entry = entry.as_hash().unwrap();
    let mut transform = identity();
    let mut operator: Option<CsgOperator> = None;
    let mut operands: Option<(ObjectIndex, ObjectIndex)> = None;
    for (key, entry) in entry.iter() {
        match key.as_str().unwrap() {
            "operator" => {
                use CsgOperator::*;
                operator = match entry.as_str().unwrap() {
                    "union" => Some(Union),
                    "intersection" => Some(Intersection),
                    "difference" => Some(Difference),
                    unknown => panic!("Unknown CSG operator: {}", unknown),
                }
            },
            "transform" => {
                transform = parse_transform(&entry);
            },
            "operands" => {
                let entries = entry.as_vec().unwrap();
                assert_eq!(entries.len(), 2);
                let l = load_object(&entries[0], world, materials);
                let r = load_object(&entries[1], world, materials);
                operands = Some((l, r));
            },
            "add" => {
                // ignore add: csg
            },
            unknown => {
                // TODO: unknown keyword
                println!("unknown csg parameter: {}", unknown);
            }
        }
    }
    if operator.is_none() {
        panic!("Missing operator of CSG");
    }
    if operands.is_none() {
        panic!("Missing operands of CSG");
    }
    world.add_csg(operator.unwrap(), operands.unwrap(), transform)
}

fn parse_point(entry: &Yaml) -> Point {
    // TODO: add error handling
    let (x, y, z) = parse_f64_tuple(&entry);
    Point::new(x, y, z)
}

fn parse_vector(entry: &Yaml) -> Vector {
    // TODO: add error handling
    let (x, y, z) = parse_f64_tuple(&entry);
    Vector::new(x, y, z)
}

fn parse_color(entry: &Yaml) -> Color {
    // TODO: add error handling
    let (r, g, b) = parse_f64_tuple(&entry);
    Color::new(r, g, b)
}

fn parse_f64_tuple(entry: &Yaml) -> (f64, f64, f64) {
    // TODO: add error handling
    let arr = entry.as_vec().unwrap();
    assert_eq!(arr.len(), 3);
    let x = parse_f64_expression(&arr[0]);
    let y = parse_f64_expression(&arr[1]);
    let z = parse_f64_expression(&arr[2]);
    (x, y, z)
}

fn parse_transform(entry: &Yaml) -> Transform {
    let mut transform = identity();
    for trans in entry.as_vec().unwrap() {
        let trans = trans.as_hash().unwrap();
        assert_eq!(trans.len(), 1);
        let (key, trans) = trans.iter().next().unwrap();
        match key.as_str().unwrap() {
            "scale" => {
                let (x, y, z) = parse_f64_tuple(trans);
                transform = scaling(x, y, z) * transform;
            },
            "translate" => {
                let dir = parse_vector(trans);
                transform = translation(&dir) * transform;
            },
            "rotate-x" => {
                let angle = parse_f64_expression(trans);
                transform = rotation_x(angle) * transform;
            },
            "rotate-y" => {
                let angle = parse_f64_expression(trans);
                transform = rotation_y(angle) * transform;
            },
            "rotate-z" => {
                let angle = parse_f64_expression(trans);
                transform = rotation_z(angle) * transform;
            },
            unknown => {
                println!("unknown transform: {}", unknown);
            }
        }
    }
    transform
}

struct GroupStore {
    groups: HashMap<String, GroupIndex>,
}

impl GroupStore {
    fn new() -> Self {
        GroupStore{ groups: HashMap::new() }
    }

    fn add(&mut self, name: &str, group: GroupIndex) {
        self.groups.insert(name.to_string(), group);
    }

    fn parse(&mut self, entry: &Yaml, world: &mut WorldBuilder) {
        let entry = entry.as_hash().unwrap();
        let mut name: Option<&str> = None;
        let mut parent: Option<GroupIndex> = None;
        let mut transform = identity();
        for (key, entry) in entry.iter() {
            match key.as_str().unwrap() {
                "name" => {
                    name = entry.as_str();
                },
                "parent" => {
                    // TODO: guard against lookup errors
                    let name = entry.as_str().unwrap();
                    parent = Some(self.groups[name]);
                },
                "transform" => {
                    transform = parse_transform(&entry);
                },
                "define" => {
                    // ignore define: group
                },
                unknown => {
                    // TODO: unknown keyword
                    println!("unknown group parameter: {}", unknown);
                }
            }
        }
        let group = world.add_group(transform);
        if let Some(parent) = parent {
            world.set_group_of(group.into(), parent);
        } else {
            world.add_object_to_scene(group.into());
        }
        if let Some(name) = name {
            self.add(name, group);
        } else {
            // TODO: group has no name
            panic!("Group without name.");
        }
    }
}

impl Index<&str> for GroupStore {
    type Output = GroupIndex;

    fn index(&self, idx: &str) -> &GroupIndex {
        &self.groups[idx]
    }
}

struct MaterialStore {
    materials: HashMap<String, Material>,
}

impl MaterialStore {
    fn new() -> Self {
        MaterialStore{ materials: HashMap::new() }
    }

    fn add(&mut self, name: &str, material: Material) {
        self.materials.insert(name.to_string(), material);
    }

    fn parse(&mut self, entry: &Yaml) {
        let mut material = Material::new();
        let mut name: Option<&str> = None;
        let entry = entry.as_hash().unwrap();
        for (key, entry) in entry.iter() {
            match key.as_str().unwrap() {
                "name" => {
                    name = entry.as_str();
                },
                "extends" => {
                    // TODO: guard against lookup errors
                    let name = entry.as_str().unwrap();
                    material = self.materials[name];
                }
                "color" => {
                    let color = parse_color(entry);
                    material.pattern = Pattern::uniform(color);
                }
                "pattern" => {
                    let pattern = parse_pattern(entry);
                    material.pattern = pattern;
                },
                // parse shader properties
                "ambient" => {
                    material.shader.ambient = parse_f64_expression(entry);
                },
                "diffuse" => {
                    material.shader.diffuse = parse_f64_expression(entry);
                },
                "specular" => {
                    material.shader.specular = parse_f64_expression(entry);
                },
                "shininess" => {
                    material.shader.shininess = parse_f64_expression(entry);
                },
                "reflective" => {
                    material.shader.reflective = parse_f64_expression(entry);
                },
                "transparency" => {
                    material.shader.transparency = parse_f64_expression(entry);
                },
                "refractive-index" => {
                    material.shader.refractive_index
                        = parse_f64_expression(entry);
                },
                "define" => {
                    // ignore define: material
                },
                unknown => {
                    // TODO: unknown keyword
                    println!("unknown material parameter: {}", unknown);
                }
            }
        }
        if let Some(name) = name {
            self.add(name, material);
        } else {
            // TODO: material has no name
        }
    }
}

impl Index<&str> for MaterialStore {
    type Output = Material;

    fn index(&self, idx: &str) -> &Material {
        &self.materials[idx]
    }
}

fn parse_pattern(entry: &Yaml) -> Pattern {
    // TODO: add error handling
    let mut pattern : Pattern = match entry["type"].as_str().unwrap() {
        "plain" | "uniform" => {
            let color = parse_color(&entry["color"]);
            Pattern::uniform(color)
        },
        "stripes" => {
            let color1 = parse_color(&entry["color1"]);
            let color2 = parse_color(&entry["color2"]);
            Pattern::stripes(color1, color2)
        },
        "rings" => {
            let color1 = parse_color(&entry["color1"]);
            let color2 = parse_color(&entry["color2"]);
            Pattern::rings(color1, color2)
        },
        "checkers" => {
            let color1 = parse_color(&entry["color1"]);
            let color2 = parse_color(&entry["color2"]);
            Pattern::checkers(color1, color2)
        },
        "gradient" => {
            let color1 = parse_color(&entry["color1"]);
            let color2 = parse_color(&entry["color2"]);
            Pattern::gradient(color1, color2)
        },
        unknown => {
            // TODO: raise error instead
            println!("unknown pattern type: {}", unknown);
            Pattern::uniform(Color::black())
        },
    };
    if entry.as_hash().unwrap().contains_key(&Yaml::from_str("transform")) {
        pattern.transform = parse_transform(&entry["transform"]);
    }
    pattern
}

fn parse_f64_expression(entry: &Yaml) -> f64 {
    // TODO: build a proper parsing tree
    if let Some(res) = entry.as_f64() {
        res
    } else if let Some(res) = entry.as_i64() {
        res as f64
    } else {
        let s = entry.as_str().unwrap();
        let mut op: Option<&str> = None;
        let mut words = s.split_whitespace();
        let mut res: f64 = parse_f64_const(words.next().unwrap());
        for w in words {
            if op.is_none() {
                op = Some(w);
            } else {
                let right = parse_f64_const(w);
                if let Some("+") = op {
                    res += right;
                } else if let Some("-") = op {
                    res -= right;
                } else if let Some("*") = op {
                    res *= right;
                } else if let Some("/") = op {
                    res /= right;
                }
            }
        }
        res
    }
}

fn parse_f64_const(w: &str) -> f64 {
    // TODO: error handling
    if let Ok(f) = w.parse::<f64>() {
        f
    } else if w == "pi" {
        std::f64::consts::PI
    } else if w == "-pi" {
        -std::f64::consts::PI
    } else {
        println!("Cannot parse as f64: {:?}", w);
        std::f64::NAN
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn loading_a_camera() {
        let s = "
- add: camera
  width: 160
  height: 120
  field-of-view: pi / 2
";
        let docs = YamlLoader::load_from_str(s).unwrap();
        let camera_yaml = &docs[0][0];
        let c = load_camera(camera_yaml);
        //assert_eq!(c.hsize, 160);
        //assert_eq!(c.vsize, 120);
        assert_eq!(c.transform, identity());
    }

    #[test]
    fn parsing_a_simple_f64() {
        let expr = &Yaml::from_str("1.");
        assert_eq!(parse_f64_expression(expr), 1.);
    }

    #[test]
    fn parsing_an_f64_expression() {
        let expr = &Yaml::from_str("1. / pi");
        let expected = 1. / std::f64::consts::PI;
        assert_relative_eq!(parse_f64_expression(expr), expected);
    }

    #[test]
    fn parsing_a_point() {
        let y = |s| Yaml::from_str(s);
        let yaml = &Yaml::Array(vec![y("pi"), y("1.5"), y("-5")]);
        let expected = Point::new(std::f64::consts::PI, 1.5, -5.);
        assert_relative_eq!(parse_point(yaml), expected);
    }
}
