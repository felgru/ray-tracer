use std::ops::Index;

use crate::color::Color;
use crate::geometry::{Point, Vector};
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

    pub fn color_at(&self, ray: &Ray) -> Color {
        let is = self.intersect(ray);
        match is.hit() {
            Some(i) => {
                let comps = self.prepare_computations(&i, ray);
                self.shade_hit(&comps)
            },
            None => Color::black(),
        }
    }

    fn intersect(&self, ray: &Ray) -> Intersections {
        let mut intersections = Intersections::new(vec![]);
        for (i, obj) in self.objects.iter().enumerate() {
            intersections.add_intersections(i, obj.intersect(ray));
        }
        intersections
    }

    fn prepare_computations<'a>(&'a self, intersection: &'a Intersection,
                                ray: &'a Ray) -> PreparedComputations<'a> {
        PreparedComputations::new(self, intersection, ray)
    }

    fn shade_hit(&self, comps: &PreparedComputations) -> Color {
        let lighting = |l| comps.object.material.lighting(l,
                                                           &comps.point,
                                                           &comps.eyev,
                                                           &comps.normalv);
        self.lights.iter().map(lighting).sum::<Color>()
    }
}

#[derive(Copy, Clone)]
struct Intersection {
    pub t: f64,
    pub object_index: usize,
}

impl Intersection {
    pub fn new(t: f64, object_index: usize) -> Self {
        Self{t, object_index}
    }
}

struct Intersections {
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

struct PreparedComputations<'a> {
    pub t: f64,
    pub object: &'a Sphere,
    pub point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub inside: bool,
}

impl<'a> PreparedComputations<'a> {
    pub fn new(world: &'a World, intersection: &'a Intersection, ray: &'a Ray) -> Self {
        let point = ray.position(intersection.t);
        let object = world.get_object(intersection.object_index);
        let eyev = -ray.direction;
        let mut normalv = object.normal_at(&point);
        let inside;
        if normalv.dot(&eyev) < 0. {
            inside = true;
            normalv = -normalv;
        } else {
            inside = false;
        }
        PreparedComputations{
            t: intersection.t,
            object: object,
            point: point,
            eyev: eyev,
            normalv: normalv,
            inside: inside,
        }
    }
}

#[cfg(test)]
pub fn default_world() -> World {
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

#[cfg(test)]
mod tests {
    use super::*;

    use crate::shapes::sphere::Sphere;

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
        assert_eq!(xs.intersections.len(), 2);
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
        assert_eq!(xs.intersections.len(), 4);
        assert_eq!(xs[0].t, 4.);
        assert_eq!(xs[1].t, 4.5);
        assert_eq!(xs[2].t, 5.5);
        assert_eq!(xs[3].t, 6.);
    }

    #[test]
    fn precompute_state_of_an_outside_hit() {
        let mut w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut shape = Sphere::new();
        let red = Color::new(1., 0., 0.);
        shape.material.color = red;
        let indx = w.add_object(shape);
        let i = Intersection::new(4., indx);
        let comps = w.prepare_computations(&i, &r);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.material.color, red);
        assert_relative_eq!(comps.point, Point::new(0., 0., -1.));
        assert_relative_eq!(comps.eyev, Vector::new(0., 0., -1.));
        assert_relative_eq!(comps.normalv, Vector::new(0., 0., -1.));
        assert_eq!(comps.inside, false);
    }

    #[test]
    fn precompute_state_of_an_inside_hit() {
        let mut w = default_world();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let mut shape = Sphere::new();
        let red = Color::new(1., 0., 0.);
        shape.material.color = red;
        let indx = w.add_object(shape);
        let i = Intersection::new(1., indx);
        let comps = w.prepare_computations(&i, &r);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.material.color, red);
        assert_relative_eq!(comps.point, Point::new(0., 0., 1.));
        assert_relative_eq!(comps.eyev, Vector::new(0., 0., -1.));
        // Normal has been inverted.
        assert_relative_eq!(comps.normalv, Vector::new(0., 0., -1.));
        assert_eq!(comps.inside, true);
    }

    #[test]
    fn shading_an_intersection() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., 0);
        let comps = w.prepare_computations(&i, &r);
        let c = w.shade_hit(&comps);
        let expected = Color::new(0.38066119308103435,
                                  0.47582649135129296,
                                  0.28549589481077575);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn shading_an_intersection_from_the_inside() {
        let mut w = default_world();
        w.lights[0] = PointLight::new(Point::new(0., 0.25, 0.),
                                      Color::new(1., 1., 1.));
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(0.5, 1);
        let comps = w.prepare_computations(&i, &r);
        let c = w.shade_hit(&comps);
        let intensity = 0.9049844720832575;
        let expected = Color::new(intensity, intensity, intensity);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn color_when_ray_misses() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let c = w.color_at(&r);
        assert_eq!(c, Color::black());
    }

    #[test]
    fn color_when_ray_hits() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let c = w.color_at(&r);
        let expected = Color::new(0.38066119308103435,
                                  0.47582649135129296,
                                  0.28549589481077575);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn color_whith_intersection_behind_the_ray() {
        let w = {
            let mut w = default_world();
            let outer = &mut w.objects[0];
            outer.material.ambient = 1.;
            let inner = &mut w.objects[1];
            inner.material.ambient = 1.;
            w
        };
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(&r);
        assert_relative_eq!(c, w.get_object(1).material.color);
    }
}
