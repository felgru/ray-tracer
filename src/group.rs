use crate::geometry::{identity, Point, Transform, Vector};
use crate::intersections::Intersections;
use crate::object::{Object, ObjectData};
use crate::object_store::ObjectStore;
use crate::rays::Ray;
use crate::shapes::bounds::Bounds;

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct Group {
    pub index: usize,
}

/// A collection of all groups
pub struct Groups {
    children: Vec<Vec<usize>>,
    parents: Vec<Group>,
    transforms: Vec<Transform>,
    bounds: Vec<Bounds>,
}

impl Groups {
    pub fn new() -> Self {
        Groups{
            children: Vec::new(),
            parents: Vec::new(),
            transforms: Vec::new(),
            bounds: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.children.len()
    }

    pub fn is_empty(&self) -> bool {
        self.children.is_empty()
    }

    pub fn add_group(&mut self) -> Group {
        Group::new(self)
    }

    pub fn add_child_to_group(&mut self, i: usize, group: Group,
                              objects: &mut [Object]) {
        self.children_mut(group).push(i);
        let mut object = &mut objects[i];
        object.parent = group;
        if let ObjectData::Group(g) = object.object {
            *self.parent_mut(g) = group;
        }
        let object_bounds = Bounds::from_points(
            &object.local_bounds(self)
                   .global_corners(object.get_transform(self)));
        self.bounds_mut(group).merge_bounds(&object_bounds);
    }

    pub fn children(&self, group: Group) -> &Vec<usize> {
        &self.children[group.index]
    }

    pub fn children_mut(&mut self, group: Group) -> &mut Vec<usize> {
        &mut self.children[group.index]
    }

    pub fn parent(&self, group: Group) -> Group {
        self.parents[group.index]
    }

    pub fn parent_mut(&mut self, group: Group) -> &mut Group {
        &mut self.parents[group.index]
    }

    pub fn transform(&self, group: Group) -> &Transform {
        &self.transforms[group.index]
    }

    pub fn transform_mut(&mut self, group: Group) -> &mut Transform {
        &mut self.transforms[group.index]
    }

    pub fn bounds(&self, group: Group) -> &Bounds {
        &self.bounds[group.index]
    }

    pub fn bounds_mut(&mut self, group: Group) -> &mut Bounds {
        &mut self.bounds[group.index]
    }
}

impl Group {
    fn new(groups: &mut Groups) -> Self {
        let nan = std::f64::NAN;
        groups.children.push(Vec::new());
        groups.parents.push(Group::with_index(std::usize::MAX));
        groups.transforms.push(identity());
        groups.bounds.push(Bounds::new(Point::new(nan, nan, nan),
                                       Point::new(nan, nan, nan)));
        let index = groups.bounds.len() - 1;
        Group {
            index,
        }
    }

    pub fn with_index(i: usize) -> Self {
        Group {index: i}
    }

    pub fn parent(&self, groups: &Groups) -> Option<Group> {
        let parent = groups.parent(*self);
        if parent.index != std::usize::MAX {
            Some(parent)
        } else {
            None
        }
    }

    pub fn intersect(&self, ray: &Ray, objects: &ObjectStore)
                                                            -> Intersections {
        let transform = objects.get_transform_of_group(*self);
        let local_ray = ray.transform(&transform.inverse());
        self.local_intersect(&local_ray, objects)
    }

    pub fn local_intersect(&self, ray: &Ray, objects: &ObjectStore)
                                                    -> Intersections {
        let mut intersections = Intersections::new();
        if !objects.get_bounds_of_group(*self).intersects(ray) {
            return intersections
        }
        for &i in objects.children_of_group(*self).iter() {
            let obj = objects.get_object(i);
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

    pub fn local_bounds(&self, groups: &Groups) -> Bounds {
        *groups.bounds(*self)
    }

    pub fn world_to_object(&self, p: &Point, groups: &Groups) -> Point {
        let p = match self.parent(groups) {
            Some(parent) => parent.world_to_object(&p, groups),
            None => *p,
        };
        groups.transform(*self).inverse() * p
    }

    pub fn normal_to_world(&self, n: &Vector, groups: &Groups) -> Vector {
        let inverse_transform = groups.transform(*self).inverse();
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        let world_normal = (inv_transp * n).normalize();
        match self.parent(groups) {
            Some(parent) => parent.normal_to_world(&world_normal, groups),
            None => world_normal,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::shapes::Shape;

    #[test]
    fn creating_an_empty_group() {
        let mut groups = Groups::new();
        let g = groups.add_group();
        assert_eq!(g.index, 0);
        assert_relative_eq!(groups.transform(g), &identity());
        groups.children(g).is_empty();
    }

    #[test]
    fn intersecting_a_ray_with_an_empty_group() {
        let mut objects = ObjectStore::new();
        let (g, _) = objects.add_group();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let is = g.local_intersect(&r, &objects);
        assert!(is.is_empty());
    }

    #[test]
    fn intersecting_a_ray_with_a_nonempty_group() {
        let mut objects = ObjectStore::new();
        let (g, _) = objects.add_group();
        objects.add_object_to_group(Shape::sphere(), g);
        let mut s1 = Shape::sphere();
        s1.set_transform(translation(&Vector::new(0., 0., -3.)));
        objects.add_object_to_group(s1, g);
        let mut s2 = Shape::sphere();
        s2.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.add_object_to_group(s2, g);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = g.local_intersect(&r, &objects);
        assert_eq!(is.len(), 4);
        assert_eq!(is[0].object_index, 2);
        assert_eq!(is[1].object_index, 2);
        assert_eq!(is[2].object_index, 1);
        assert_eq!(is[3].object_index, 1);
    }

    #[test]
    fn intersecting_a_ray_with_a_transformed_group() {
        let mut objects = ObjectStore::new();
        let (g, g_i) = objects.add_group();
        objects.set_transform_of_object(g_i, scaling(2., 2., 2.));
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.add_object_to_group(s, g);
        let r = Ray::new(Point::new(10., 0., -10.), Vector::new(0., 0., 1.));
        let is = g.intersect(&r, &objects);
        assert_eq!(is.len(), 2);
    }

    #[test]
    fn bounding_box_of_group_with_a_single_object() {
        let mut groups = Groups::new();
        let g = groups.add_group();
        *groups.transform_mut(g) = scaling(2., 2., 2.);
        let mut objects: Vec<Object> = Vec::new();
        objects.push(Shape::sphere().into());
        groups.add_child_to_group(0, g, &mut objects);
        let b = g.local_bounds(&groups);
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(1., 1., 1.)));
    }

    #[test]
    fn bounding_box_of_group_with_two_objects() {
        let mut groups = Groups::new();
        let g = groups.add_group();
        *groups.transform_mut(g) = scaling(2., 2., 2.);
        let mut objects: Vec<Object> = Vec::new();
        objects.push(Shape::sphere().into());
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        groups.add_child_to_group(0, g, &mut objects);
        groups.add_child_to_group(1, g, &mut objects);
        let b = g.local_bounds(&groups);
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(6., 1., 1.)));
    }
}
