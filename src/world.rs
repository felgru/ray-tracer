use std::ops::Index;

use crate::light::PointLight;
use crate::rays::Ray;
use crate::shapes::sphere::Sphere;

pub struct World {
    objects: Vec<Sphere>,
    lights: Vec<PointLight>,
}

impl World {
    pub fn new() -> Self {
        World{
            objects: Vec::new(),
            lights: Vec::new(),
        }
    }

    pub fn add_object(&mut self, obj: Sphere) -> usize {
        self.objects.push(obj);
        self.objects.len() - 1
    }

    pub fn add_light(&mut self, light: PointLight) -> usize {
        self.lights.push(light);
        self.lights.len() - 1
    }

    pub fn get_object(&self, i: usize) -> &Sphere {
        &self.objects[i]
    }

    pub fn get_light(&self, i: usize) -> &PointLight {
        &self.lights[i]
    }

    pub fn intersect(&self, ray: &Ray) -> Intersections {
        let mut intersections = Intersections::new(vec![]);
        for (i, obj) in self.objects.iter().enumerate() {
            intersections.add_intersections(i, obj.intersect(ray));
        }
        intersections
    }
}

#[derive(Copy, Clone)]
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
    pub fn new(intersections: Vec<Intersection>) -> Self {
        let mut intersections = intersections;
        intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
        Intersections {
            intersections
        }
    }

    pub fn count(&self) -> usize {
        self.intersections.len()
    }

    pub fn hit(&self) -> Option<Intersection> {
        self.intersections.iter().find(|i| i.t >= 0.).map(|&i| i)
    }

    pub fn add_intersections(&mut self, obj_index: usize, intersections: Vec<f64>) {
        self.intersections.extend(intersections.into_iter().map(|t| Intersection::new(t, obj_index)));
        self.intersections.sort_by(|a, b| a.t.partial_cmp(&b.t).unwrap());
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

    use crate::shapes::sphere::Sphere;
    use crate::color::Color;
    use crate::geometry::{Point, Vector};

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
        let xs = Intersections::new(vec![i1, i2]);
        assert_eq!(xs.count(), 2);
        assert_eq!(xs[0].t, 1.);
        assert_eq!(xs[1].t, 2.);
    }

    #[test]
    fn hit_when_all_intersections_have_positive_t() {
        let i1 = Intersection::new(1., 0);
        let i2 = Intersection::new(2., 0);
        let xs = Intersections::new(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.object_index, 0);
    }

    #[test]
    fn hit_when_some_intersections_have_negative_t() {
        let i1 = Intersection::new(-1., 0);
        let i2 = Intersection::new(1., 0);
        let xs = Intersections::new(vec![i2, i1]);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 1.);
        assert_eq!(i.object_index, 0);
    }

    #[test]
    fn hit_when_all_intersections_have_negative_t() {
        let i1 = Intersection::new(-2., 0);
        let i2 = Intersection::new(-1., 0);
        let xs = Intersections::new(vec![i2, i1]);
        let i = xs.hit();
        assert!(i.is_none());
    }

    #[test]
    fn hit_is_always_lowest_nonnegative_intersection() {
        let intersections: Vec<Intersection> = [5., 7., -3., 2.].into_iter().map(|&t| Intersection::new(t, 0)).collect();
        let xs = Intersections::new(intersections);
        let i = xs.hit().unwrap();
        assert_eq!(i.t, 2.);
    }

    #[test]
    fn create_a_world() {
        let w = World::new();
        assert!(w.objects.is_empty());
        assert!(w.lights.is_empty());
    }

    fn default_world() -> World {
        let intensity = Color::new(1., 1., 1.);
        let light = PointLight::new(Point::new(-10., 10., -10.),
                                    intensity);
        let mut s1 = Sphere::new();
        let c = Color::new(0.8, 1., 0.6);
        s1.material.color = c;
        s1.material.diffuse = 0.7;
        s1.material.specular = 0.2;
        let mut s2 = Sphere::new();
        use crate::geometry::scaling;
        let t = scaling(0.5, 0.5, 0.5);
        s2.transform = t;

        let mut w = World::new();
        w.add_light(light);
        w.add_object(s1);
        w.add_object(s2);
        w
    }

    #[test]
    fn let_there_be_light() {
        let intensity = Color::new(1., 1., 1.);
        let c = Color::new(0.8, 1., 0.6);
        use crate::geometry::scaling;
        let t = scaling(0.5, 0.5, 0.5);

        let w = default_world();
        assert_eq!(w.lights[0].intensity, intensity);
        assert_eq!(w.objects[0].material.color, c);
        assert_eq!(w.objects[1].transform, t);
    }

    #[test]
    fn intersect_a_world_with_a_ray() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = w.intersect(&r);
        assert_eq!(xs.count(), 4);
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.);
    }
}
