use crate::intersections::Intersections;
use crate::object_store::{ObjectIndex, ObjectStore};
use crate::rays::Ray;
use crate::shapes::{bounds::Bounds, ShapeIndex};

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct CsgIndex {
    pub index: usize,
}

impl CsgIndex {
    pub fn new(index: usize) -> Self {
        CsgIndex{index}
    }
}

#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub enum CsgOperator {
    Union,
    Intersection,
    Difference,
}

/// A collection of all CSGs
pub struct CSGs {
    operators: Vec<CsgOperator>,
    children: Vec<(ObjectIndex, ObjectIndex)>,
    bounds: Vec<Bounds>,
}

impl CSGs {
    pub fn new() -> Self {
        CSGs{
            operators: Vec::new(),
            children: Vec::new(),
            bounds: Vec::new(),
        }
    }

    pub fn len(&self) -> usize {
        self.operators.len()
    }

    pub fn is_empty(&self) -> bool {
        self.operators.is_empty()
    }

    pub fn add_csg(&mut self, operator: CsgOperator,
                   left: ObjectIndex, right: ObjectIndex) -> CsgIndex {
        self.operators.push(operator);
        self.children.push((left, right));
        self.bounds.push(Bounds::empty());
        let index = self.bounds.len() - 1;
        CsgIndex::new(index)
    }

    pub fn children(&self, csg: CsgIndex) -> (ObjectIndex, ObjectIndex) {
        self.children[csg.index]
    }

    pub fn local_bounds_of(&self, group: CsgIndex) -> &Bounds {
        &self.bounds[group.index]
    }

    pub fn bounds_mut(&mut self, csg: CsgIndex) -> &mut Bounds {
        &mut self.bounds[csg.index]
    }

    pub fn local_intersect(&self, csg: CsgIndex, ray: &Ray,
                           objects: &ObjectStore) -> Intersections {
        if !self.local_bounds_of(csg).intersects(ray) {
            return Intersections::new();
        }

        let (left, right) = self.children[csg.index];
        let mut intersections = objects.intersect(left, ray);
        let right_intersections = objects.intersect(right, ray);

        let mut left_shapes: Vec<ShapeIndex>
            = intersections.iter().map(|&i| i.shape_index).collect();
        left_shapes.sort();

        intersections.add_intersections_from(right_intersections);

        filter_intersections(self.operators[csg.index], &left_shapes,
                             intersections)
    }
}

fn filter_intersections(op: CsgOperator,
                        left_indices: &[ShapeIndex],
                        intersections: Intersections) -> Intersections {
    // begin outside of both children
    let mut inl = false;
    let mut inr = false;

    // prepare a list to receive the filtered intersections
    let mut result = Intersections::new();

    for i in intersections {
        // if i.shape_index is part of the left child, then lhit is true
        let lhit = left_indices.binary_search(&i.shape_index).is_ok();

        if intersection_allowed(op, lhit, inl, inr) {
            result.add_intersection(i);
        }

        // depending on which object was hit, toggle either inl or inr
        if lhit {
            inl = !inl;
        } else {
            inr = !inr;
        }
    }

    result
}

fn intersection_allowed(op: CsgOperator, lhit: bool, inl: bool, inr: bool)
                                                                    -> bool {
    match op {
        CsgOperator::Union => {
            (lhit && !inr) || (!lhit && !inl)
        },
        CsgOperator::Intersection => {
            (lhit && inr) || (!lhit && inl)
        },
        CsgOperator::Difference => {
            (lhit && !inr) || (!lhit && inl)
        },
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::object_store::Parent;
    use crate::shapes::Shape;

    #[test]
    fn constructing_a_csg() {
        let mut objects = ObjectStore::new();
        let l = objects.add_shape(Shape::sphere(), identity());
        let r = objects.add_shape(Shape::cube(), identity());
        let c = objects.add_csg(CsgOperator::Union, (l, r), identity());
        assert_eq!(c.index, 0);
        assert_eq!(objects.csgs().children[c.index],
                   (ObjectIndex::Shape(l), ObjectIndex::Shape(r)));
        assert_eq!(objects.parent_of_object(l), Parent::CSG(c));
        assert_eq!(objects.parent_of_object(r), Parent::CSG(c));
    }

    #[test]
    fn evaluating_the_rule_for_a_csg_operation() {
        use CsgOperator::*;
        assert_eq!(intersection_allowed(Union, true, true, true), false);
        assert_eq!(intersection_allowed(Union, true, true, false), true);
        assert_eq!(intersection_allowed(Union, true, false, true), false);
        assert_eq!(intersection_allowed(Union, true, false, false), true);
        assert_eq!(intersection_allowed(Union, false, true, true), false);
        assert_eq!(intersection_allowed(Union, false, true, false), false);
        assert_eq!(intersection_allowed(Union, false, false, true), true);
        assert_eq!(intersection_allowed(Union, false, false, false), true);

        assert_eq!(intersection_allowed(Intersection, true, true, true),
                   true);
        assert_eq!(intersection_allowed(Intersection, true, true, false),
                   false);
        assert_eq!(intersection_allowed(Intersection, true, false, true),
                   true);
        assert_eq!(intersection_allowed(Intersection, true, false, false),
                   false);
        assert_eq!(intersection_allowed(Intersection, false, true, true),
                   true);
        assert_eq!(intersection_allowed(Intersection, false, true, false),
                   true);
        assert_eq!(intersection_allowed(Intersection, false, false, true),
                   false);
        assert_eq!(intersection_allowed(Intersection, false, false, false),
                   false);

        assert_eq!(intersection_allowed(Difference, true, true, true),
                   false);
        assert_eq!(intersection_allowed(Difference, true, true, false),
                   true);
        assert_eq!(intersection_allowed(Difference, true, false, true),
                   false);
        assert_eq!(intersection_allowed(Difference, true, false, false),
                   true);
        assert_eq!(intersection_allowed(Difference, false, true, true),
                   true);
        assert_eq!(intersection_allowed(Difference, false, true, false),
                   true);
        assert_eq!(intersection_allowed(Difference, false, false, true),
                   false);
        assert_eq!(intersection_allowed(Difference, false, false, false),
                   false);
    }

    #[test]
    fn filtering_a_list_of_intersections() {
        use CsgOperator::*;
        check_filtering_of_intersections(Union, 0, 3);
        check_filtering_of_intersections(Intersection, 1, 2);
        check_filtering_of_intersections(Difference, 0, 1);
    }

    fn check_filtering_of_intersections(op: CsgOperator,
                                        x0: usize, x1: usize) {
        let mut objects = ObjectStore::new();
        let l = objects.add_shape(Shape::sphere(), identity());
        let r = objects.add_shape(Shape::cube(), identity());
        let c = objects.add_csg(op, (l, r), identity());
        objects.set_bounds_of(c.into());
        use crate::intersections::Intersection;
        let xs = vec![Intersection::new(1., l),
                      Intersection::new(2., r),
                      Intersection::new(3., l),
                      Intersection::new(4., r)];
        let left_indices = [l];
        let result = filter_intersections(op, &left_indices,
                                          Intersections::from_vec(xs.clone()));
        assert_eq!(result.len(), 2);
        assert_eq!(result[0], xs[x0]);
        assert_eq!(result[1], xs[x1]);
    }

    #[test]
    fn ray_misses_csg_object() {
        let mut objects = ObjectStore::new();
        let l = objects.add_shape(Shape::sphere(), identity());
        let r = objects.add_shape(Shape::cube(), identity());
        let c = objects.add_csg(CsgOperator::Union, (l, r), identity());
        let ray = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let xs = objects.csgs().local_intersect(c, &ray, &objects);
        assert!(xs.is_empty());
    }

    #[test]
    fn ray_hits_a_csg_object() {
        let mut objects = ObjectStore::new();
        let l = objects.add_shape(Shape::sphere(), identity());
        let r = objects.add_shape(Shape::sphere(),
                                  translation(&Vector::new(0., 0., 0.5)));
        let c = objects.add_csg(CsgOperator::Union, (l, r), identity());
        objects.set_bounds_of(c.into());
        let ray = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = objects.csgs().local_intersect(c, &ray, &objects);
        assert_eq!(xs.len(), 2);
        assert_relative_eq!(xs[0].t, 4.);
        assert_eq!(xs[0].shape_index, l);
        assert_relative_eq!(xs[1].t, 6.5);
        assert_eq!(xs[1].shape_index, r);
    }
}
