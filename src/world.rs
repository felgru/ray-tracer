use std::ops::Index;

use crate::color::Color;
use crate::geometry::{Point, reflect, Vector};
use crate::light::PointLight;
use crate::rays::Ray;
use crate::shapes::Shape;

pub struct World {
    objects: Vec<Shape>,
    lights: Vec<PointLight>,
}

impl World {
    pub fn new() -> Self {
        World{
            objects: Vec::new(),
            lights: Vec::new(),
        }
    }

    pub fn add_object(&mut self, obj: Shape) -> usize {
        self.objects.push(obj);
        self.objects.len() - 1
    }

    pub fn add_light(&mut self, light: PointLight) -> usize {
        self.lights.push(light);
        self.lights.len() - 1
    }

    pub fn get_object(&self, i: usize) -> &Shape {
        &self.objects[i]
    }

    pub fn get_light(&self, i: usize) -> &PointLight {
        &self.lights[i]
    }

    pub fn color_at(&self, ray: &Ray, remaining: usize) -> Color {
        let is = self.intersect(ray);
        match is.hit() {
            Some(i) => {
                let comps = self.prepare_computations(&i, ray);
                self.shade_hit(&comps, remaining)
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

    fn shade_hit(&self, comps: &PreparedComputations,
                 remaining: usize) -> Color {
        let lighting = |l| comps.object.get_material()
                                .lighting(l,
                                          &comps.point,
                                          &comps.eyev,
                                          &comps.normalv,
                                          self.is_shadowed(&comps.over_point, l));
        let surface = self.lights.iter().map(lighting).sum::<Color>();
        let reflected = self.reflected_color(comps, remaining);

        surface + reflected
    }

    fn reflected_color(&self, comps: &PreparedComputations,
                       remaining: usize) -> Color {
        let reflectivity = comps.object.get_material().reflective;
        if remaining <= 0 || reflectivity == 0. {
            return Color::black();
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(&reflect_ray, remaining - 1);

        color * reflectivity
    }

    fn is_shadowed(&self, point: &Point, light: &PointLight) -> bool {
        let v = light.position - point;
        let distance = v.norm();
        let direction = v.normalize();

        let ray = Ray::new(*point, direction);
        let intersections = self.intersect(&ray);

        match intersections.hit() {
            Some(i) => i.t < distance,
            None    => false,
        }
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
    pub object: &'a Shape,
    pub point: Point,
    pub over_point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub reflectv: Vector,
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
        let eps = Self::over_point_eps();
        let over_point = point + normalv * eps;
        PreparedComputations{
            t: intersection.t,
            object: object,
            point: point,
            over_point: over_point,
            eyev: eyev,
            normalv: normalv,
            reflectv: reflect(&ray.direction, &normalv),
            inside: inside,
        }
    }

    fn over_point_eps() -> f64 {
        1e-10
    }
}

#[cfg(test)]
pub fn default_world() -> World {
    let intensity = Color::white();
    let light = PointLight::new(Point::new(-10., 10., -10.),
                                intensity);
    let mut s1 = Shape::sphere();
    let m = {
        let c = Color::new(0.8, 1., 0.6);
        use crate::material::Material;
        let mut m = Material::new();
        m.color = c;
        m.diffuse = 0.7;
        m.specular = 0.2;
        m
    };
    s1.set_material(m);

    let mut s2 = Shape::sphere();
    use crate::geometry::scaling;
    s2.set_transform(scaling(0.5, 0.5, 0.5));

    let mut w = World::new();
    w.add_light(light);
    w.add_object(s1);
    w.add_object(s2);
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::material::Material;

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
        let intensity = Color::white();
        let c = Color::new(0.8, 1., 0.6);
        use crate::geometry::scaling;
        let t = scaling(0.5, 0.5, 0.5);

        let w = default_world();
        assert_eq!(w.lights[0].intensity, intensity);
        assert_eq!(w.objects[0].get_material().color, c);
        assert_eq!(w.objects[1].get_transform(), &t);
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
        let mut shape = Shape::sphere();
        let red = Color::new(1., 0., 0.);
        let mut m = Material::new();
        m.color = red;
        shape.set_material(m);
        let indx = w.add_object(shape);
        let i = Intersection::new(4., indx);
        let comps = w.prepare_computations(&i, &r);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.get_material().color, red);
        assert_relative_eq!(comps.point, Point::new(0., 0., -1.));
        assert_relative_eq!(comps.eyev, Vector::new(0., 0., -1.));
        assert_relative_eq!(comps.normalv, Vector::new(0., 0., -1.));
        assert_eq!(comps.inside, false);
    }

    #[test]
    fn precompute_state_of_an_inside_hit() {
        let mut w = default_world();
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let mut shape = Shape::sphere();
        let red = Color::new(1., 0., 0.);
        let mut m = shape.get_material().clone();
        m.color = red;
        shape.set_material(m);
        let indx = w.add_object(shape);
        let i = Intersection::new(1., indx);
        let comps = w.prepare_computations(&i, &r);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.get_material().color, red);
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
        let c = w.shade_hit(&comps, 4);
        let expected = Color::new(0.38066119308103435,
                                  0.47582649135129296,
                                  0.28549589481077575);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn shading_an_intersection_from_the_inside() {
        let mut w = default_world();
        w.lights[0] = PointLight::new(Point::new(0., 0.25, 0.),
                                      Color::white());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(0.5, 1);
        let comps = w.prepare_computations(&i, &r);
        let c = w.shade_hit(&comps, 4);
        let intensity = 0.9049844720832575;
        let expected = Color::new(intensity, intensity, intensity);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn color_when_ray_misses() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 1., 0.));
        let c = w.color_at(&r, 4);
        assert_eq!(c, Color::black());
    }

    #[test]
    fn color_when_ray_hits() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let c = w.color_at(&r, 4);
        let expected = Color::new(0.38066119308103435,
                                  0.47582649135129296,
                                  0.28549589481077575);
        assert_relative_eq!(c, expected);
    }

    #[test]
    fn color_whith_intersection_behind_the_ray() {
        let w = {
            let mut w = default_world();
            set_ambient_of_object(&mut w, 0, 1.);
            set_ambient_of_object(&mut w, 1, 1.);
            w
        };
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(&r, 4);
        assert_relative_eq!(c, w.get_object(1).get_material().color);
    }

    fn set_ambient_of_object(w: &mut World, i: usize, ambient: f64) {
        let obj = &mut w.objects[i];
        let mut m = obj.get_material().clone();
        m.ambient = ambient;
        obj.set_material(m);
    }

    #[test]
    fn no_shadow_when_nothing_is_between_point_and_light() {
        let w = default_world();
        let l = w.get_light(0);
        let p = Point::new(0., 10., 0.);
        assert!(!w.is_shadowed(&p, l));
    }

    #[test]
    fn shadow_when_object_between_point_and_light() {
        let w = default_world();
        let l = w.get_light(0);
        let p = Point::new(10., -10., 10.);
        assert!(w.is_shadowed(&p, l));
    }

    #[test]
    fn no_shadow_when_object_behind_light() {
        let w = default_world();
        let l = w.get_light(0);
        let p = Point::new(-20., 20., -20.);
        assert!(!w.is_shadowed(&p, l));
    }

    #[test]
    fn no_shadow_when_object_behind_point() {
        let w = default_world();
        let l = w.get_light(0);
        let p = Point::new(-2., 2., -2.);
        assert!(!w.is_shadowed(&p, l));
    }

    #[test]
    fn shade_hit_is_given_an_intersection_in_shadow() {
        let mut w = World::new();
        w.add_light(PointLight::new(Point::new(0., 0., -10.),
                                    Color::white()));
        let s1 = Shape::sphere();
        w.add_object(s1);
        let mut s2 = Shape::sphere();
        s2.set_transform(translation(&Vector::new(0., 0., 10.)));
        w.add_object(s2);
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., 1);
        let comps = w.prepare_computations(&i, &r);
        let c = w.shade_hit(&comps, 4);
        assert_relative_eq!(c, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn hit_should_offset_the_point() {
        let w = {
            let mut w = World::new();
            let mut shape = Shape::sphere();
            shape.set_transform(translation(&Vector::new(0., 0., 1.)));
            w.add_object(shape);
            w
        };
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(5., 0);
        let comps = w.prepare_computations(&i, &r);
        let eps = PreparedComputations::over_point_eps();
        assert!(comps.over_point[2] < -eps/2.);
        assert!(comps.point[2] > comps.over_point[2]);
    }

    #[test]
    fn precomputing_the_reflection_vector() {
        let mut w = World::new();
        w.add_object(Shape::plane());
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 1., -1.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 0);
        let comps = w.prepare_computations(&i, &r);
        assert_relative_eq!(comps.reflectv,
                            Vector::new(0., 1./sqrt2, 1./sqrt2));
    }

    #[test]
    fn reflected_color_of_a_nonreflective_material() {
        let mut w = default_world();
        set_ambient_of_object(&mut w, 1, 1.);
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., 1);
        let comps = w.prepare_computations(&i, &r);
        let color = w.reflected_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn reflected_color_of_a_reflective_material() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let comps = w.prepare_computations(&i, &r);
        let color = w.reflected_color(&comps, 4);
        let expected = Color::new(0.19032, 0.2379, 0.14274);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn reflected_color_at_maximum_recursion_depth() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let comps = w.prepare_computations(&i, &r);
        let color = w.reflected_color(&comps, 0);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn shade_hit_with_a_reflective_material() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let comps = w.prepare_computations(&i, &r);
        let color = w.shade_hit(&comps, 4);
        let expected = Color::new(0.87677, 0.92436, 0.82918);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn terminiation_in_case_of_infinite_reflection() {
        let mut w = World::new();
        let light = PointLight::new(Point::new(0., 0., 0.), Color::white());
        w.add_light(light);
        let mut m = Material::new();
        m.reflective = 1.;
        let mut lower = Shape::plane();
        lower.set_material(m.clone());
        lower.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(lower);
        let mut upper = Shape::plane();
        upper.set_material(m.clone());
        upper.set_transform(translation(&Vector::new(0., 1., 0.)));
        w.add_object(upper);
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        let _ = w.color_at(&r, 4); // should terminate successfully
    }
}
