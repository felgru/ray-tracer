use std::ops::Index;
use crate::shapes::ShapeIndex;

#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Intersection {
    pub t: f64,
    pub shape_index: ShapeIndex,
}

impl Intersection {
    pub fn new(t: f64, shape_index: ShapeIndex) -> Self {
        Self{t, shape_index}
    }
}

pub struct Intersections {
    intersections: Vec<Intersection>,
}

impl Intersections {
    pub fn new() -> Self {
        Intersections { intersections: Vec::new() }
    }

    pub fn from_vec(intersections: Vec<Intersection>) -> Self {
        let mut intersections = intersections;
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        Intersections {
            intersections
        }
    }

    pub fn hit(&self) -> Option<Intersection> {
        self.intersections.iter().find(|i| i.t >= 0.).copied()
    }

    pub fn add_intersection(&mut self, intersection: Intersection) {
        let t = &intersection.t;
        let pos = self.intersections.binary_search_by(|i| i.t.partial_cmp(t)
                                                             .unwrap())
                                    .unwrap_or_else(|e| e);
        self.intersections.insert(pos, intersection);
    }

    pub fn add_intersections(&mut self, shp_index: ShapeIndex,
                             intersections: Vec<f64>) {
        self.intersections.extend(intersections.into_iter()
                                  .map(|t| Intersection::new(t, shp_index)));
        self.intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }

    pub fn add_intersections_from(&mut self, intersections: Self) {
        self.intersections.extend(intersections);
        self.intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }

    pub fn iter(&self) -> std::slice::Iter<Intersection> {
        self.intersections.iter()
    }

    pub fn len(&self) -> usize {
        self.intersections.len()
    }

    pub fn is_empty(&self) -> bool {
        self.intersections.is_empty()
    }
}

impl IntoIterator for Intersections {
    type Item = Intersection;
    type IntoIter = ::std::vec::IntoIter<Self::Item>;

    fn into_iter(self) -> Self::IntoIter {
        self.intersections.into_iter()
    }
}

impl Index<usize> for Intersections {
    type Output = Intersection;

    fn index(&self, i: usize) -> &Intersection {
        &self.intersections[i]
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersection_properties() {
        let i = Intersection::new(3.5, ShapeIndex::new(0));
        assert_eq!(i.t, 3.5);
        assert_eq!(i.shape_index, ShapeIndex::new(0));
    }

    #[test]
    fn aggregating_intersections() {
        let i1 = Intersection::new(1., ShapeIndex::new(0));
        let i2 = Intersection::new(2., ShapeIndex::new(0));
        let xs = Intersections::from_vec(vec![i1, i2]);
        assert_eq!(xs.intersections.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_when_all_intersections_have_positive_t() {
        let i1 = Intersection::new(1., ShapeIndex::new(0));
        let i2 = Intersection::new(2., ShapeIndex::new(0));
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.shape_index, ShapeIndex::new(0));
    }

    #[test]
    fn hit_when_some_intersections_have_negative_t() {
        let i1 = Intersection::new(-1., ShapeIndex::new(0));
        let i2 = Intersection::new(1., ShapeIndex::new(0));
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.shape_index, ShapeIndex::new(0));
    }

    #[test]
    fn hit_when_all_intersections_have_negative_t() {
        let i1 = Intersection::new(-2., ShapeIndex::new(0));
        let i2 = Intersection::new(-1., ShapeIndex::new(0));
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit();
        assert!(i.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let intersections: Vec<Intersection>
            = [5., 7., -3., 2.].iter()
                               .map(|&t| Intersection::new(t,
                                            ShapeIndex::new(0))).collect();
        let xs = Intersections::from_vec(intersections);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 2.);
    }
}
