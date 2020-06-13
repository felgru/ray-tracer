use crate::intersections::Intersections;
use crate::object_store::{ObjectIndex, ObjectStore};
use crate::rays::Ray;
use crate::shapes::bounds::Bounds;

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct GroupIndex {
    pub index: usize,
}

/// A collection of all groups
pub struct Groups {
    children: Vec<Vec<ObjectIndex>>,
    bounds: Vec<Bounds>,
}

impl Groups {
    pub fn new() -> Self {
        Groups{
            children: Vec::new(),
            bounds: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.children.len()
    }

    pub fn is_empty(&self) -> bool {
        self.children.is_empty()
    }

    pub fn add_group(&mut self) -> GroupIndex {
        self.children.push(Vec::new());
        self.bounds.push(Bounds::empty());
        let index = self.bounds.len() - 1;
        GroupIndex::new(index)
    }

    pub fn child_iter(&self, group: GroupIndex) -> GroupChildIter {
        GroupChildIter::new(group)
    }

    pub fn children(&self, group: GroupIndex) -> &Vec<ObjectIndex> {
        &self.children[group.index]
    }

    pub fn children_mut(&mut self, group: GroupIndex) -> &mut Vec<ObjectIndex> {
        &mut self.children[group.index]
    }

    pub fn local_bounds_of(&self, group: GroupIndex) -> &Bounds {
        &self.bounds[group.index]
    }

    pub fn bounds_mut(&mut self, group: GroupIndex) -> &mut Bounds {
        &mut self.bounds[group.index]
    }

    pub fn local_intersect(&self, group: GroupIndex, ray: &Ray,
                           objects: &ObjectStore) -> Intersections {
        let mut intersections = Intersections::new();
        if !self.local_bounds_of(group).intersects(ray) {
            return intersections;
        }
        for &child in self.children(group).iter() {
            let is = objects.intersect(child, ray);
            intersections.add_intersections_from(is);
        }
        intersections
    }
}

impl GroupIndex {
    pub fn new(index: usize) -> Self {
        GroupIndex{index}
    }
}

pub struct GroupChildIter {
    group: GroupIndex,
    child: usize,
}

impl GroupChildIter {
    fn new(group: GroupIndex) -> Self {
        Self{group, child: 0}
    }

    pub fn next(&mut self, groups: &Groups) -> Option<ObjectIndex> {
        let children = groups.children(self.group);
        if children.len() <= self.child { return None; }
        let obj = children[self.child];
        self.child += 1;
        Some(obj)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::object_store::Parent;
    use crate::shapes::{Shape, ShapeIndex};

    #[test]
    fn creating_an_empty_group() {
        let mut obs = ObjectStore::new();
        let g = obs.add_group(identity());
        assert_eq!(g.index, 0);
        let obj = ObjectIndex::Group(g);
        assert_relative_eq!(obs.get_transform_of_object(obj), &identity());
        obs.groups().children(g).is_empty();
    }

    #[test]
    fn adding_a_child_to_a_group() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(identity());
        assert_eq!(g.index, 0);
        let i = objects.add_shape_to_group(Shape::sphere(), identity(), g);
        assert!(!objects.groups().children(g).is_empty());
        assert_eq!(i, ShapeIndex::new(0));
        let i = objects.groups().children(g)[0];
        assert_eq!(i, ObjectIndex::Shape(ShapeIndex::new(0)));
        assert_eq!(objects.parent_of_object(i),
                   Parent::Group(g));
    }

    #[test]
    fn intersecting_a_ray_with_an_empty_group() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(identity());
        objects.set_bounds_of(g.into());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let is = objects.groups().local_intersect(g, &r, &objects);
        assert!(is.is_empty());
    }

    #[test]
    fn intersecting_a_ray_with_a_nonempty_group() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(identity());
        let s0 = objects.add_shape_to_group(Shape::sphere(), identity(), g);
        let s1_transform = translation(&Vector::new(0., 0., -3.));
        let s1 = objects.add_shape_to_group(Shape::sphere(), s1_transform, g);
        let s2_transform = translation(&Vector::new(5., 0., 0.));
        let _s2 = objects.add_shape_to_group(Shape::sphere(), s2_transform, g);
        objects.set_bounds_of(g.into());
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = objects.groups().local_intersect(g, &r, &objects);
        assert_eq!(is.len(), 4);
        assert_eq!(is[0].shape_index, s1);
        assert_eq!(is[1].shape_index, s1);
        assert_eq!(is[2].shape_index, s0);
        assert_eq!(is[3].shape_index, s0);
    }

    #[test]
    fn intersecting_a_ray_with_a_transformed_group() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(scaling(2., 2., 2.));
        let s = Shape::sphere();
        let s_transform = translation(&Vector::new(5., 0., 0.));
        objects.add_shape_to_group(s, s_transform, g);
        objects.set_bounds_of(g.into());
        let r = Ray::new(Point::new(10., 0., -10.), Vector::new(0., 0., 1.));
        let is = objects.intersect(g.into(), &r);
        assert_eq!(is.len(), 2);
    }

    #[test]
    fn bounding_box_of_group_with_a_single_object() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(scaling(2., 2., 2.));
        objects.add_shape_to_group(Shape::sphere(), identity(), g);
        objects.set_bounds_of(g.into());
        let b = *objects.get_bounds_of_group(g);
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(1., 1., 1.)));
    }

    #[test]
    fn bounding_box_of_group_with_two_objects() {
        let mut objects = ObjectStore::new();
        let g = objects.add_group(scaling(2., 2., 2.));
        objects.add_shape_to_group(Shape::sphere(), identity(), g);
        objects.add_shape_to_group(Shape::sphere(),
                                   translation(&Vector::new(5., 0., 0.)), g);
        objects.set_bounds_of(g.into());
        let b = *objects.get_bounds_of_group(g);
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(6., 1., 1.)));
    }
}
