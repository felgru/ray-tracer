use std::ops::Index;

#[derive(Copy, Clone, PartialEq)]
pub struct Intersection {
    pub t: f64,
    pub object_index: usize,
}

impl Intersection {
    pub fn new(t: f64, object_index: usize) -> Self {
        Self{t, object_index}
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
        self.intersections.iter().find(|i| i.t >= 0.).map(|&i| i)
    }

    pub fn add_intersections(&mut self, obj_index: usize, intersections: Vec<f64>) {
        self.intersections.extend(intersections.into_iter().map(|t| Intersection::new(t, obj_index)));
        self.intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
    }

    pub fn iter(&self) -> std::slice::Iter<Intersection> {
        self.intersections.iter()
    }

    pub fn len(&self) -> usize {
        self.intersections.len()
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
        let i = Intersection::new(3.5, 0);
        assert_eq!(i.t, 3.5);
        assert_eq!(i.object_index, 0);
    }

    #[test]
    fn aggregating_intersections() {
        let i1 = Intersection::new(1., 0);
        let i2 = Intersection::new(2., 0);
        let xs = Intersections::from_vec(vec![i1, i2]);
        assert_eq!(xs.intersections.len(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_when_all_intersections_have_positive_t() {
        let i1 = Intersection::new(1., 0);
        let i2 = Intersection::new(2., 0);
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.object_index, 0);
    }

    #[test]
    fn hit_when_some_intersections_have_negative_t() {
        let i1 = Intersection::new(-1., 0);
        let i2 = Intersection::new(1., 0);
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.object_index, 0);
    }

    #[test]
    fn hit_when_all_intersections_have_negative_t() {
        let i1 = Intersection::new(-2., 0);
        let i2 = Intersection::new(-1., 0);
        let xs = Intersections::from_vec(vec![i2, i1]);
        let i = xs.hit();
        assert!(i.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let intersections: Vec<Intersection>
            = [5., 7., -3., 2.].into_iter()
                               .map(|&t| Intersection::new(t, 0)).collect();
        let xs = Intersections::from_vec(intersections);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 2.);
    }
}
