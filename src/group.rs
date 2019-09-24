use crate::color::Color;
use crate::geometry::{identity, Point, Transform, Vector};
use crate::intersections::Intersections;
use crate::light::PointLight;
use crate::material::Material;
use crate::rays::Ray;
use crate::shapes::Shape;

pub struct Object {
    parent: usize,
    object: ObjectData,
}

impl Object {
    pub fn normal_at(&self, world_point: &Point, objects: &Vec<Object>)
                                                                    -> Vector {
        if let Some(shp) = self.as_shape() {
            let local_point = self.world_to_object(world_point, objects);
            let local_normal = shp.local_normal_at(&local_point);
            self.normal_to_world(&local_normal, objects)
        } else {
            panic!("Trying to take normal of non-shape object!");
        }
    }

    pub fn color_at(&self, world_point: &Point, objects: &Vec<Object>)
                                                                    -> Color {
        let pattern = match self.as_shape() {
            Some(shp) => shp.get_material().pattern,
            None => panic!("Trying to take color of non-shape object!"),
        };
        let object_point = self.world_to_object(world_point, objects);
        let pattern_point = pattern.transform.inverse() * object_point;
        pattern.local_color_at(&pattern_point)
    }

    pub fn lighting(&self, light: &PointLight,
                    position: &Point, eyev: &Vector, normalv: &Vector,
                    in_shadow: bool,
                    objects: &Vec<Object>) -> Color {
        let color = self.color_at(position, objects);
        let shader = self.get_material().shader;
        shader.lighting(color, light, position, eyev, normalv, in_shadow)
    }

    pub fn intersect(&self, ray: &Ray, obj_id: usize, objects: &Vec<Object>)
                                                            -> Intersections {
        use ObjectData::*;
        match &self.object {
            Shape(shp) => {
                let mut intersections = Intersections::new();
                intersections.add_intersections(obj_id, shp.intersect(ray));
                intersections
            }
            Group(grp) => grp.intersect(ray, objects),
        }
    }

    fn parent<'a>(&self, objects: &'a Vec<Object>) -> Option<&'a Object> {
        if self.parent != std::usize::MAX {
            Some(&objects[self.parent])
        } else {
            None
        }
    }

    pub fn get_transform(&self) -> &Transform {
        self.object.get_transform()
    }

    pub fn set_transform(&mut self, transform: Transform) {
        self.object.set_transform(transform);
    }

    pub fn get_material(&self) -> &Material {
        if let Some(shp) = self.as_shape() {
            shp.get_material()
        } else {
            panic!("Trying to take material of non-shape object!")
        }
    }

    pub fn set_material(&mut self, material: Material) {
        match &mut self.object {
            ObjectData::Shape(shp) => shp.set_material(material),
            _ => panic!("Trying to set material of non-shape object!"),
        }
    }

    fn as_shape(&self) -> Option<&Shape> {
        match &self.object {
            ObjectData::Shape(shp) => Some(&shp),
            _ => None,
        }
    }

    fn world_to_object(&self, p: &Point, objects: &Vec<Object>) -> Point {
        let p = match self.parent(objects) {
            Some(parent) => parent.world_to_object(&p, objects),
            None => p.clone(),
        };
        self.get_transform().inverse() * p
    }

    fn normal_to_world(&self, n: &Vector, objects: &Vec<Object>)
                                                                -> Vector {
        let inverse_transform = self.get_transform().inverse();
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        let world_normal = (inv_transp * n).normalize();
        match self.parent(objects) {
            Some(parent) => parent.normal_to_world(&world_normal, objects),
            None => world_normal,
        }
    }
}

impl From<Shape> for Object {
    fn from(s: Shape) -> Self {
        Object {
            parent: std::usize::MAX,
            object: ObjectData::Shape(s),
        }
    }
}

impl From<Group> for Object {
    fn from(g: Group) -> Self {
        Object {
            parent: std::usize::MAX,
            object: ObjectData::Group(g),
        }
    }
}

enum ObjectData {
    Group(Group),
    Shape(Shape),
}

impl ObjectData {
    fn get_transform(&self) -> &Transform {
        use ObjectData::*;
        match self {
            Group(grp) => &grp.transform,
            Shape(shp) => shp.get_transform(),
        }
    }

    fn set_transform(&mut self, transform: Transform) {
        use ObjectData::*;
        match self {
            Group(grp) => grp.transform = transform,
            Shape(shp) => shp.set_transform(transform),
        }
    }
}

pub struct Group {
    pub index: usize,
    pub children: Vec<usize>,
    pub transform: Transform,
}

impl Group {
    pub fn new(index: usize) -> Self {
        Group {
            index: index,
            children: Vec::new(),
            transform: identity(),
        }
    }

    pub fn add_child(&mut self, i: usize, objects: &mut Vec<Object>) {
        objects[i].parent = self.index;
        self.children.push(i);
    }

    pub fn get_child<'a>(&self, i: usize, objects: &'a Vec<Object>) -> &'a Object {
        &objects[self.children[i]]
    }

    pub fn intersect(&self, ray: &Ray, objects: &Vec<Object>)
                                                            -> Intersections {
        let local_ray = ray.transform(&self.transform.inverse());
        self.local_intersect(&local_ray, objects)
    }

    pub fn local_intersect(&self, ray: &Ray, objects: &Vec<Object>)
                                                            -> Intersections {
        let mut intersections = Intersections::new();
        for &i in self.children.iter() {
            let obj = &objects[i];
            match &obj.object {
                ObjectData::Shape(shp) => {
                    intersections.add_intersections(i, shp.intersect(ray));
                },
                ObjectData::Group(grp) => {
                    let is = grp.intersect(ray, objects);
                    intersections.add_intersections_from(is);
                }
            }
        }
        intersections
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::patterns::Pattern;
    use crate::shapes::Shape;

    #[test]
    fn creating_an_empty_group() {
        let g = Group::new(1);
        assert_eq!(g.index, 1);
        assert_relative_eq!(g.transform, identity());
        g.children.is_empty();
    }

    #[test]
    fn adding_a_child_to_a_group() {
        let mut objects: Vec<Object> = vec![Shape::sphere().into()];
        let mut g = Group::new(1);
        assert_eq!(g.index, 1);
        g.add_child(0, &mut objects);
        assert!(!g.children.is_empty());
        assert_eq!(g.children[0], 0);
        assert_eq!(g.get_child(0, &objects).parent, 1);
    }

    #[test]
    fn intersecting_a_ray_with_an_empty_group() {
        let objects: Vec<Object> = Vec::new();
        let g = Group::new(1);
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let is = g.local_intersect(&r, &objects);
        assert_eq!(is.len(), 0);
    }

    #[test]
    fn intersecting_a_ray_with_a_nonempty_group() {
        let mut g = Group::new(3);
        let mut objects: Vec<Object> = Vec::new();
        objects.push(Shape::sphere().into());
        let mut s1 = Shape::sphere();
        s1.set_transform(translation(&Vector::new(0., 0., -3.)));
        objects.push(s1.into());
        let mut s2 = Shape::sphere();
        s2.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s2.into());
        g.add_child(0, &mut objects);
        g.add_child(1, &mut objects);
        g.add_child(2, &mut objects);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = g.local_intersect(&r, &objects);
        assert_eq!(is.len(), 4);
        assert_eq!(is[0].object_index, 1);
        assert_eq!(is[1].object_index, 1);
        assert_eq!(is[2].object_index, 0);
        assert_eq!(is[3].object_index, 0);
    }

    #[test]
    fn intersecting_a_ray_with_a_transformed_group() {
        let mut g = Group::new(3);
        g.transform = scaling(2., 2., 2.);
        let mut objects: Vec<Object> = Vec::new();
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        g.add_child(0, &mut objects);
        let r = Ray::new(Point::new(10., 0., -10.), Vector::new(0., 0., 1.));
        let is = g.intersect(&r, &objects);
        assert_eq!(is.len(), 2);
    }

    #[test]
    fn converting_a_point_from_world_to_object_space() {
        let mut objects: Vec<Object> = Vec::new();
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        let mut g1 = Group::new(2);
        g1.transform = rotation_y(std::f64::consts::PI/2.);
        let mut g2 = Group::new(1);
        g2.transform = scaling(2., 2., 2.);
        g2.add_child(0, &mut objects);
        objects.push(g2.into());
        g1.add_child(1, &mut objects);
        objects.push(g1.into());
        let s = &objects[0];
        let p = s.world_to_object(&Point::new(-2., 0., -10.), &objects);
        assert_relative_eq!(p, Point::new(0., 0., -1.));
    }

    #[test]
    fn converting_a_normal_form_object_to_world_space() {
        let mut objects: Vec<Object> = Vec::new();
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        let mut g1 = Group::new(2);
        g1.transform = rotation_y(std::f64::consts::PI/2.);
        let mut g2 = Group::new(1);
        g2.transform = scaling(1., 2., 3.);
        g2.add_child(0, &mut objects);
        objects.push(g2.into());
        g1.add_child(1, &mut objects);
        objects.push(g1.into());
        let s = &objects[0];
        let x = f64::sqrt(3.) / 3.;
        let p = s.normal_to_world(&Vector::new(x, x, x), &objects);
        assert_relative_eq!(p, Vector::new(0.2857, 0.4286, -0.8571),
                            max_relative = 1e-4);
    }

    #[test]
    fn finding_the_normal_on_a_child_object() {
        let mut objects: Vec<Object> = Vec::new();
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        let mut g1 = Group::new(2);
        g1.transform = rotation_y(std::f64::consts::PI/2.);
        let mut g2 = Group::new(1);
        g2.transform = scaling(1., 2., 3.);
        g2.add_child(0, &mut objects);
        objects.push(g2.into());
        g1.add_child(1, &mut objects);
        objects.push(g1.into());
        let s = &objects[0];
        let sqrt3 = f64::sqrt(3.);
        let p = s.normal_at(&Point::new(sqrt3, 2./3.*sqrt3, -5. - sqrt3/3.),
                            &objects);
        assert_relative_eq!(p, Vector::new(0.2857, 0.4286, -0.8571),
                            max_relative = 1e-4);
    }

    #[test]
    fn normal_of_translated_sphere() {
        let mut s: Object = Shape::sphere().into();
        use crate::geometry::translation;
        s.set_transform(translation(&Vector::new(0., 1., 0.)));
        let n = s.normal_at(&Point::new(0., 1.70711, -0.70711), &vec![]);
        assert_relative_eq!(n, Vector::new(0., 0.70711, -0.70711), max_relative = 1e-5);
    }

    #[test]
    fn normal_of_transformed_sphere() {
        let mut s: Object = Shape::sphere().into();
        use crate::geometry::{scaling, rotation_z};
        s.set_transform(scaling(1., 0.5, 1.) * rotation_z(std::f64::consts::PI/5.));
        let x = 1./f64::sqrt(2.);
        let n = s.normal_at(&Point::new(0., x, -x), &vec![]);
        assert_relative_eq!(n, Vector::new(0., 0.97014, -0.24254), max_relative = 1e-4);
    }

    #[test]
    fn stripes_with_transformed_sphere() {
        let mut object: Object = Shape::sphere().into();
        object.set_transform(scaling(2., 2., 2.));
        let pattern = Pattern::stripes(Color::white(), Color::black());
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.), &vec![]);
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern() {
        let mut object: Object = Shape::sphere().into();
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = scaling(2., 2., 2.);
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.), &vec![]);
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern_and_object() {
        let mut object: Object = Shape::sphere().into();
        object.set_transform(scaling(2., 2., 2.));
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = translation(&Vector::new(0.5, 0., 0.));
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(2.5, 0., 0.), &vec![]);
        assert_eq!(c, Color::white());
    }

    fn set_pattern_of_object(object: &mut Object, pattern: Pattern) {
        let mut material = object.get_material().clone();
        material.pattern = pattern;
        object.set_material(material);
    }

    #[test]
    fn lighting_with_stripe_pattern() {
        use crate::patterns::Pattern;
        let mut m = Material::new();
        m.pattern = Pattern::stripes(Color::white(), Color::black());
        m.shader.ambient = 1.;
        m.shader.diffuse = 0.;
        m.shader.specular = 0.;
        let mut obj: Object = Shape::sphere().into();
        obj.set_material(m);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let c1 = obj.lighting(&light, &Point::new(0.9, 0., 0.),
                              &eyev, &normalv, false,
                              &vec![]);
        let c2 = obj.lighting(&light, &Point::new(1.1, 0., 0.),
                              &eyev, &normalv, false,
                              &vec![]);
        assert_relative_eq!(c1, Color::white());
        assert_relative_eq!(c2, Color::black());
    }
}
