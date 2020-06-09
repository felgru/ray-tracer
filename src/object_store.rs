use crate::color::Color;
use crate::geometry::{Point, Transform, Vector};
use crate::group::{GroupIndex, Groups};
use crate::intersections::Intersections;
use crate::light::PointLight;
use crate::material::Material;
use crate::rays::Ray;
use crate::shapes::bounds::Bounds;
use crate::shapes::{Shape, ShapeIndex};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum ObjectIndex {
    Group(GroupIndex),
    Shape(ShapeIndex),
}

impl ObjectIndex {
    pub fn new_group(index: usize) -> ObjectIndex {
        ObjectIndex::Group(GroupIndex::new(index))
    }

    pub fn new_shape(index: usize) -> ObjectIndex {
        ObjectIndex::Shape(ShapeIndex::new(index))
    }
}

impl From<ShapeIndex> for ObjectIndex {
    fn from(s: ShapeIndex) -> Self {
        ObjectIndex::Shape(s)
    }
}

impl From<GroupIndex> for ObjectIndex {
    fn from(g: GroupIndex) -> Self {
        ObjectIndex::Group(g)
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum Parent {
    Group(GroupIndex),
    None,
}

impl From<GroupIndex> for Parent {
    fn from(g: GroupIndex) -> Self {
        Parent::Group(g)
    }
}

struct Parents {
    group_parents: Vec<Parent>,
    shape_parents: Vec<Parent>,
}

impl Parents {
    fn new() -> Parents {
        Parents{
            group_parents: Vec::new(),
            shape_parents: Vec::new(),
        }
    }

    fn parent_of(&self, object: ObjectIndex) -> Parent {
        match object {
            ObjectIndex::Group(g) => self.group_parents[g.index],
            ObjectIndex::Shape(s) => self.shape_parents[s.index],
        }
    }

    fn add_group_parent(&mut self) {
        self.group_parents.push(Parent::None);
    }

    fn add_shape_parent(&mut self) {
        self.shape_parents.push(Parent::None);
    }

    fn set_parent_of(&mut self, obj: ObjectIndex, parent: Parent) {
        use ObjectIndex::*;
        match obj {
            Group(g) => {
                self.group_parents[g.index] = parent;
            },
            Shape(s) => {
                self.shape_parents[s.index] = parent;
            },
        }
    }
}

struct Transforms {
    group_transforms: Vec<Transform>,
    shape_transforms: Vec<Transform>,
}

impl Transforms {
    pub fn new() -> Self {
        Transforms{
            group_transforms: Vec::new(),
            shape_transforms: Vec::new(),
        }
    }

    pub fn add_group_transform(&mut self, transform: Transform) {
        self.group_transforms.push(transform);
    }

    pub fn add_shape_transform(&mut self, transform: Transform) {
        self.shape_transforms.push(transform);
    }

    pub fn get_transform_of(&self, i: ObjectIndex) -> &Transform {
        match i {
            ObjectIndex::Group(g) => {
                &self.group_transforms[g.index]
            }
            ObjectIndex::Shape(s) => {
                &self.shape_transforms[s.index]
            }
        }
    }

    pub fn get_transform_mut_of(&mut self, i: ObjectIndex) -> &mut Transform {
        match i {
            ObjectIndex::Group(g) => {
                &mut self.group_transforms[g.index]
            }
            ObjectIndex::Shape(s) => {
                &mut self.shape_transforms[s.index]
            }
        }
    }

    pub fn set_transform_of(&mut self, i: ObjectIndex, transform: Transform) {
        *self.get_transform_mut_of(i) = transform;
    }
}

pub struct ObjectStore {
    // TODO: adapt type of shapes?
    shapes: Vec<Shape>,
    groups: Groups,
    parents: Parents,
    transforms: Transforms,
}

impl ObjectStore {
    pub fn new() -> Self {
        ObjectStore{
            shapes: Vec::new(),
            groups: Groups::new(),
            parents: Parents::new(),
            transforms: Transforms::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.shapes.is_empty()
            && self.groups.is_empty()
    }

    pub fn add_shape(&mut self, shape: Shape, transform: Transform)
                                                    -> ShapeIndex {
        self.parents.add_shape_parent();
        self.transforms.add_shape_transform(transform);
        self.shapes.push(shape);
        ShapeIndex::new(self.shapes.len() - 1)
    }

    pub fn add_group(&mut self, transform: Transform) -> GroupIndex {
        self.parents.add_group_parent();
        self.transforms.add_group_transform(transform);
        self.groups.add_group()
    }

    pub fn add_subgroup(&mut self, transform: Transform, group: GroupIndex)
                                                                -> GroupIndex {
        let sub_group = self.add_group(transform);
        let i = ObjectIndex::Group(sub_group);
        self.add_child_to_group(i, group);
        self.parents.set_parent_of(i, Parent::Group(group));
        sub_group
    }

    pub fn add_shape_to_group(&mut self, shape: Shape, transform: Transform,
                              group: GroupIndex) -> ShapeIndex {
        let i = self.add_shape(shape, transform);
        let obj = ObjectIndex::Shape(i);
        self.add_child_to_group(obj, group);
        self.parents.set_parent_of(obj, group.into());
        i
    }

    fn add_child_to_group(&mut self, obj: ObjectIndex, group: GroupIndex) {
        self.groups.children_mut(group).push(obj);
        let child_transform = self.get_transform_of_object(obj);
        let child_bounds = Bounds::from_points(
            &self.local_bounds_of(obj)
                 .global_corners(child_transform));
        self.groups.bounds_mut(group).merge_bounds(&child_bounds);
        // TODO: merge the updated bounds into the parents bounds
    }

    pub fn intersect(&self, i: ObjectIndex, ray: &Ray) -> Intersections {
        let transform = self.get_transform_of_object(i);
        let local_ray = ray.transform(&transform.inverse());
        match i {
            ObjectIndex::Shape(shp) => {
                let mut intersections = Intersections::new();
                intersections.add_intersections(shp,
                    self.shapes[shp.index].local_intersect(&local_ray));
                intersections
            }
            ObjectIndex::Group(grp) => {
                self.groups.local_intersect(grp, &local_ray, &self)
            }
        }
    }

    pub fn get_bounds_of_group(&self, group: GroupIndex) -> &Bounds {
        self.groups.local_bounds_of(group)
    }

    pub fn parent_of_object<O: Into<ObjectIndex>>(&self, obj: O) -> Parent {
        self.parents.parent_of(obj.into())
    }

    #[cfg(test)]
    pub fn groups(&self) -> &Groups {
        &self.groups
    }

    pub fn get_transform_of_object(&self, i: ObjectIndex) -> &Transform {
        self.transforms.get_transform_of(i)
    }

    pub fn set_transform_of_object(&mut self, i: ObjectIndex,
                                   transform: Transform) {
        self.transforms.set_transform_of(i, transform);
    }

    pub fn get_material_of(&self, shp: ShapeIndex) -> &Material {
        self.shapes[shp.index].get_material()
    }

    pub fn set_material_of(&mut self, shp: ShapeIndex, material: Material) {
        self.shapes[shp.index].set_material(material)
    }

    fn color_at(&self, shp: ShapeIndex, world_point: &Point) -> Color {
        let pattern = self.get_material_of(shp).pattern;
        let shp_obj = ObjectIndex::Shape(shp);
        let object_point = self.world_to_object(shp_obj, world_point);
        let pattern_point = pattern.transform.inverse() * object_point;
        pattern.local_color_at(&pattern_point)
    }

    pub fn lighting(&self, shp: ShapeIndex, light: &PointLight,
                    position: &Point, eyev: &Vector, normalv: &Vector,
                    in_shadow: bool) -> Color {
        let color = self.color_at(shp, position);
        let shader = self.get_material_of(shp).shader;
        shader.lighting(color, light, position, eyev, normalv, in_shadow)
    }

    pub fn world_to_object(&self, obj: ObjectIndex, p: &Point) -> Point {
        let p = match self.parents.parent_of(obj) {
            Parent::Group(parent) => {
                self.world_to_object(ObjectIndex::Group(parent), &p)
            }
            Parent::None => *p,
        };
        self.get_transform_of_object(obj).inverse() * p
    }

    pub fn normal_at(&self, shp: ShapeIndex, world_point: &Point) -> Vector {
        let obj = ObjectIndex::Shape(shp);
        let shp = &self.shapes[shp.index];
        let local_point = self.world_to_object(obj,
                                               world_point);
        let local_normal = shp.local_normal_at(&local_point);
        self.normal_to_world(obj, &local_normal)
    }

    pub fn normal_to_world(&self, obj: ObjectIndex, n: &Vector) -> Vector {
        let inverse_transform = self.get_transform_of_object(obj).inverse();
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        let world_normal = (inv_transp * n).normalize();
        match self.parents.parent_of(obj) {
            Parent::Group(parent) => {
                self.normal_to_world(ObjectIndex::Group(parent),
                                     &world_normal)
            }
            Parent::None => world_normal,
        }
    }

    pub fn local_bounds_of(&self, obj: ObjectIndex) -> Bounds {
        use ObjectIndex::*;
        match obj {
            Shape(shp) => self.shapes[shp.index].local_bounds(),
            Group(grp) => *self.groups.local_bounds_of(grp),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::*;
    use crate::material::Material;
    use crate::patterns::Pattern;

    #[test]
    fn create_an_object_store() {
        let obs = ObjectStore::new();
        assert!(obs.shapes.is_empty());
        assert!(obs.groups.is_empty());
    }

    #[test]
    fn one_object() {
        let mut obs = ObjectStore::new();
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::test();
        s.set_material(m);
        let i = obs.add_shape(s, identity());
        assert_eq!(i, ShapeIndex::new(0));
        assert_eq!(obs.shapes.len(), 1);
        assert!(obs.groups.is_empty());
        assert_eq!(obs.get_material_of(i).pattern,
                   Pattern::test());
    }

    #[test]
    fn normal_of_translated_sphere() {
        let mut obs = ObjectStore::new();
        use crate::geometry::translation;
        let si = obs.add_shape(Shape::sphere(),
                               translation(&Vector::new(0., 1., 0.)));
        let n = obs.normal_at(si, &Point::new(0., 1.70711, -0.70711));
        assert_relative_eq!(n, Vector::new(0., 0.70711, -0.70711),
                            max_relative = 1e-5);
    }

    #[test]
    fn normal_of_transformed_sphere() {
        let mut obs = ObjectStore::new();
        use crate::geometry::{scaling, rotation_z};
        let transform = scaling(1., 0.5, 1.)
                        * rotation_z(std::f64::consts::PI/5.);
        let si = obs.add_shape(Shape::sphere(), transform);
        let x = 1./f64::sqrt(2.);
        let n = obs.normal_at(si, &Point::new(0., x, -x));
        assert_relative_eq!(n, Vector::new(0., 0.97014, -0.24254),
                            max_relative = 1e-4);
    }

    // Shape tests

    #[test]
    fn stripes_with_transformed_sphere() {
        let mut obs = ObjectStore::new();
        let si = obs.add_shape(Shape::sphere(), scaling(2., 2., 2.));
        let pattern = Pattern::stripes(Color::white(), Color::black());
        set_pattern_of_shape(&mut obs, si, pattern);
        let c = obs.color_at(si, &Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern() {
        let mut obs = ObjectStore::new();
        let si = obs.add_shape(Shape::sphere(), identity());
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = scaling(2., 2., 2.);
        set_pattern_of_shape(&mut obs, si, pattern);
        let c = obs.color_at(si, &Point::new(1.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern_and_object() {
        let mut obs = ObjectStore::new();
        let si = obs.add_shape(Shape::sphere(), scaling(2., 2., 2.));
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = translation(&Vector::new(0.5, 0., 0.));
        set_pattern_of_shape(&mut obs, si, pattern);
        let c = obs.color_at(si, &Point::new(2.5, 0., 0.));
        assert_eq!(c, Color::white());
    }

    fn set_pattern_of_shape(objects: &mut ObjectStore, si: ShapeIndex,
                            pattern: Pattern) {
        let mut material = objects.get_material_of(si).clone();
        material.pattern = pattern;
        objects.set_material_of(si, material);
    }

    #[test]
    fn lighting_with_stripe_pattern() {
        use crate::patterns::Pattern;
        let mut m = Material::new();
        m.pattern = Pattern::stripes(Color::white(), Color::black());
        m.shader.ambient = 1.;
        m.shader.diffuse = 0.;
        m.shader.specular = 0.;
        let mut s = Shape::sphere();
        s.set_material(m);
        let mut objects = ObjectStore::new();
        let si = objects.add_shape(s, identity());
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let c1 = objects.lighting(si, &light, &Point::new(0.9, 0., 0.),
                                  &eyev, &normalv, false);
        let c2 = objects.lighting(si, &light, &Point::new(1.1, 0., 0.),
                                  &eyev, &normalv, false);
        assert_relative_eq!(c1, Color::white());
        assert_relative_eq!(c2, Color::black());
    }
}
