use crate::geometry::{identity, Point, Transform};
use crate::intersections::Intersections;
use crate::object::{Object, ObjectData};
use crate::rays::Ray;
use crate::shapes::bounds::Bounds;

pub struct Group {
    pub index: usize,
    pub children: Vec<usize>,
    pub transform: Transform,
    bounds: Bounds,
}

impl Group {
    pub fn new(index: usize) -> Self {
        let nan = std::f64::NAN;
        Group {
            index: index,
            children: Vec::new(),
            transform: identity(),
            bounds: Bounds::new(Point::new(nan, nan, nan),
                                Point::new(nan, nan, nan)),
        }
    }

    pub fn add_child(&mut self, i: usize, objects: &mut Vec<Object>) {
        self.children.push(i);
        let mut object = &mut objects[i];
        object.parent = self.index;
        let object_bounds = Bounds::from_points(
            &object.local_bounds().global_corners(object.get_transform()));
        self.bounds.merge_bounds(&object_bounds);
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
        if !self.bounds.intersects(ray) {
            return intersections
        }
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

    pub fn local_bounds(&self) -> Bounds {
        self.bounds.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
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
    fn bounding_box_of_group_with_a_single_object() {
        let mut g = Group::new(3);
        g.transform = scaling(2., 2., 2.);
        let mut objects: Vec<Object> = Vec::new();
        objects.push(Shape::sphere().into());
        g.add_child(0, &mut objects);
        let b = g.local_bounds();
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(1., 1., 1.)));
    }

    #[test]
    fn bounding_box_of_group_with_two_objects() {
        let mut g = Group::new(3);
        g.transform = scaling(2., 2., 2.);
        let mut objects: Vec<Object> = Vec::new();
        objects.push(Shape::sphere().into());
        let mut s = Shape::sphere();
        s.set_transform(translation(&Vector::new(5., 0., 0.)));
        objects.push(s.into());
        g.add_child(0, &mut objects);
        g.add_child(1, &mut objects);
        let b = g.local_bounds();
        assert_relative_eq!(b, Bounds::new(Point::new(-1., -1., -1.),
                                           Point::new(6., 1., 1.)));
    }
}
