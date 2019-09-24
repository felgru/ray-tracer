use crate::color::Color;
use crate::geometry::{Point, reflect, Vector};
use crate::intersections::{Intersection, Intersections};
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
                let comps = self.prepare_computations(&i, ray, &is);
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

    fn prepare_computations<'a, 'b>(&'a self, hit: &'b Intersection,
                                    ray: &'b Ray,
                                    intersections: &'b Intersections)
                                                -> PreparedComputations<'a> {
        PreparedComputations::new(self, hit, ray, intersections)
    }

    fn shade_hit(&self, comps: &PreparedComputations,
                 remaining: usize) -> Color {
        let lighting = |l| comps.object
                                .lighting(l,
                                          &comps.point,
                                          &comps.eyev,
                                          &comps.normalv,
                                          self.is_shadowed(&comps.over_point, l));
        let surface = self.lights.iter().map(lighting).sum::<Color>();

        let reflected = self.reflected_color(comps, remaining);
        let refracted = self.refracted_color(comps, remaining);

        let shader = comps.object.get_material().shader;
        if shader.reflective > 0. && shader.transparency > 0. {
            let reflectance = schlick(&comps);
            surface + reflected * reflectance
                    + refracted * (1. - reflectance)
        } else {
            surface + reflected + refracted
        }
    }

    fn reflected_color(&self, comps: &PreparedComputations,
                       remaining: usize) -> Color {
        let reflectivity = comps.object.get_material().shader.reflective;
        if remaining <= 0 || reflectivity == 0. {
            return Color::black();
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(&reflect_ray, remaining - 1);

        color * reflectivity
    }

    fn refracted_color(&self, comps: &PreparedComputations,
                       remaining: usize) -> Color {
        if remaining == 0
           || comps.object.get_material().shader.transparency == 0. {
            return Color::black();
        }

        // Snell's Law
        let n_ratio = comps.n1 / comps.n2;
        let cos_i = comps.eyev.dot(&comps.normalv);
        // Find sin(theta_t)^2 via trigonometric identity
        let sin2_t = n_ratio*n_ratio * (1. - cos_i*cos_i);
        if sin2_t > 1. {
            // Total internal reflection
            return Color::black();
        }

        // Find cos(theta_t) via trigonometric identity
        let cos_t = (1.0 - sin2_t).sqrt();

        let direction = comps.normalv * (n_ratio * cos_i - cos_t)
                        - comps.eyev * n_ratio;
        let refract_ray = Ray::new(comps.under_point, direction);

        self.color_at(&refract_ray, remaining - 1)
            * comps.object.get_material().shader.transparency
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

struct PreparedComputations<'a> {
    pub t: f64,
    pub object: &'a Shape,
    pub point: Point,
    pub over_point: Point,
    pub under_point: Point,
    pub eyev: Vector,
    pub normalv: Vector,
    pub reflectv: Vector,
    pub n1: f64,
    pub n2: f64,
    pub inside: bool,
}

impl<'a> PreparedComputations<'a> {
    pub fn new<'b>(world: &'a World, hit: &'b Intersection, ray: &'b Ray,
                   intersections: &'b Intersections) -> Self {
        let point = ray.position(hit.t);
        let object = world.get_object(hit.object_index);
        let eyev = -ray.direction;
        let mut normalv = object.normal_at(&point);
        let inside;
        if normalv.dot(&eyev) < 0. {
            inside = true;
            normalv = -normalv;
        } else {
            inside = false;
        }
        let eps = Self::collision_offset();
        let over_point = point + normalv * eps;
        let under_point = point - normalv * eps;
        let (n1, n2) = Self::get_refractive_indices_at(world, hit,
                                                       intersections);
        PreparedComputations{
            t: hit.t,
            object: object,
            point: point,
            over_point: over_point,
            under_point: under_point,
            eyev: eyev,
            normalv: normalv,
            reflectv: reflect(&ray.direction, &normalv),
            n1,
            n2,
            inside: inside,
        }
    }

    fn get_refractive_indices_at(world: &World, hit: &Intersection,
                                 xs: &Intersections) -> (f64, f64) {
        let mut visited_objects = Vec::<usize>::new();

        for i in xs.iter().take_while(|&i| i != hit) {
            let idx = i.object_index;
            if let Some(pos) = visited_objects.iter().position(|p| *p == idx) {
                visited_objects.remove(pos);
            } else {
                visited_objects.push(idx);
            }
        }

        let refractive_index_of = |i| {
            world.get_object(i).get_material()
                               .shader.refractive_index
        };

        let (n1, n2);

        let hit_idx = hit.object_index;
        if let Some(&idx) = visited_objects.last() {
            n1 = refractive_index_of(idx);
            if let Some(pos) = visited_objects.iter()
                                              .position(|p| *p == hit_idx) {
                visited_objects.remove(pos);
                if let Some(&idx) = visited_objects.last() {
                    n2 = refractive_index_of(idx);
                } else {
                    n2 = 1.;
                }
            } else {
                n2 = refractive_index_of(hit_idx);
            }
        } else {
            n1 = 1.;
            n2 = refractive_index_of(hit_idx);
        }
        (n1, n2)
    }

    fn collision_offset() -> f64 {
        1e-10
    }
}

/// Schlick's approximation to the Fresnel effect
fn schlick(comps: &PreparedComputations) -> f64 {
    let mut cos = comps.eyev.dot(&comps.normalv);

    // total internal reflection can only occur if n1 > n2
    if comps.n1 > comps.n2 {
        let n = comps.n1 / comps.n2;
        let sin2_t = n*n * (1. - cos*cos);
        if sin2_t > 1. {
            return 1.;
        }

        let cos_t = (1. - sin2_t).sqrt();

        // when n1 > n2, use cos(theta_t) instead
        cos = cos_t;
    }

    let r0 = ((comps.n1 - comps.n2) / (comps.n1 + comps.n2)).powi(2);
    r0 + (1. - r0) * (1. - cos).powi(5)
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
        use crate::patterns::Pattern;
        m.pattern = Pattern::uniform(c);
        m.shader.diffuse = 0.7;
        m.shader.specular = 0.2;
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
    use crate::patterns::Pattern;

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
        assert_eq!(w.objects[0].get_material().pattern, Pattern::uniform(c));
        assert_eq!(w.objects[1].get_transform(), &t);
    }

    #[test]
    fn intersect_a_world_with_a_ray() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let xs = w.intersect(&r);
        assert_eq!(xs.len(), 4);
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
        m.pattern = Pattern::uniform(red);
        shape.set_material(m);
        let indx = w.add_object(shape);
        let i = Intersection::new(4., indx);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.get_material().pattern, Pattern::uniform(red));
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
        m.pattern = Pattern::uniform(red);
        shape.set_material(m);
        let indx = w.add_object(shape);
        let i = Intersection::new(1., indx);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(comps.object.get_material().pattern, Pattern::uniform(red));
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
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
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
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
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
        assert_relative_eq!(c, Color::white());
    }

    fn set_ambient_of_object(w: &mut World, i: usize, ambient: f64) {
        let obj = &mut w.objects[i];
        let mut m = obj.get_material().clone();
        m.shader.ambient = ambient;
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
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let c = w.shade_hit(&comps, 4);
        assert_relative_eq!(c, Color::new(0.1, 0.1, 0.1));
    }

    #[test]
    fn over_point_accuracy() {
        let w = offset_test_world();
        let comps = offset_test_computations(&w);
        let eps = PreparedComputations::collision_offset();
        assert!(comps.over_point[2] < -eps/2., "assert {} < {}", comps.over_point[2], -eps/2.);
        assert!(comps.point[2] > comps.over_point[2], "assert {} > {}", comps.point[2], comps.over_point[2]);
    }

    #[test]
    fn under_point_accuracy() {
        let w = offset_test_world();
        let comps = offset_test_computations(&w);
        let eps = PreparedComputations::collision_offset();
        assert!(comps.under_point[2] > eps/2., "assert {} > {}", comps.under_point[2], eps/2.);
        assert!(comps.point[2] < comps.under_point[2], "assert {} < {}", comps.point[2], comps.under_point[2]);
    }

    fn offset_test_world() -> World {
        let mut w = World::new();
        let mut shape = Shape::sphere();
        shape.set_transform(translation(&Vector::new(0., 0., 1.)));
        w.add_object(shape);
        w
    }

    fn offset_test_computations<'a>(w: &World) -> PreparedComputations {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(5., 0);
        let is = Intersections::new(vec![i]);
        w.prepare_computations(&i, &r, &is)
    }

    #[test]
    fn precomputing_the_reflection_vector() {
        let mut w = World::new();
        w.add_object(Shape::plane());
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 1., -1.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 0);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.reflectv,
                            Vector::new(0., 1./sqrt2, 1./sqrt2));
    }

    #[test]
    fn reflected_color_of_a_nonreflective_material() {
        let mut w = default_world();
        set_ambient_of_object(&mut w, 1, 1.);
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., 1);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.reflected_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn reflected_color_of_a_reflective_material() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.shader.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.reflected_color(&comps, 4);
        let expected = Color::new(0.19032, 0.2379, 0.14274);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn reflected_color_at_maximal_recursion_depth() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.shader.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.reflected_color(&comps, 0);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn shade_hit_with_a_reflective_material() {
        let mut w = default_world();
        let mut shape = Shape::plane();
        let mut m = Material::new();
        m.shader.reflective = 0.5;
        shape.set_material(m);
        shape.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(shape);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
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
        m.shader.reflective = 1.;
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

    fn glass_sphere() -> Shape {
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        s.set_material(m);
        s
    }

    fn refraction_test_sphere(refractive_index: f64,
                              transform: Transform) -> Shape {
        let mut s = glass_sphere();
        s.set_transform(transform);
        let mut m = *s.get_material();
        m.shader.refractive_index = refractive_index;
        s.set_material(m);
        s
    }

    #[test]
    fn test_refraction() {
        refraction_test(0, 1.0, 1.5);
        refraction_test(1, 1.5, 2.0);
        refraction_test(2, 2.0, 2.5);
        refraction_test(3, 2.5, 2.5);
        refraction_test(4, 2.5, 1.5);
        refraction_test(5, 1.5, 1.0);
    }

    fn refraction_test(index: usize, n1: f64, n2: f64) {
        let a = refraction_test_sphere(1.5, scaling(2., 2., 2.));
        let b = refraction_test_sphere(2.0,
                        translation(&Vector::new(0., 0., -0.25)));
        let c = refraction_test_sphere(2.5,
                        translation(&Vector::new(0., 0., 0.25)));
        let mut w = World::new();
        w.add_object(a);
        w.add_object(b);
        w.add_object(c);
        let r = Ray::new(Point::new(0., 0., -4.), Vector::new(0., 0., 1.));
        let xs = w.intersect(&r);
        let comps = w.prepare_computations(&xs[index], &r, &xs);
        assert_eq!(comps.n1, n1);
        assert_eq!(comps.n2, n2);
    }

    #[test]
    fn refracted_color_of_opaque_surface() {
        let w = default_world();
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = Intersections::new(vec![Intersection::new(4., 0),
                                         Intersection::new(6., 0)]);
        let comps = w.prepare_computations(&is[0], &r, &is);
        let color = w.refracted_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_at_maximal_recursion_depth() {
        let mut w = default_world();
        let s = &mut w.objects[0];
        let mut m = s.get_material().clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        s.set_material(m);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let is = Intersections::new(vec![Intersection::new(4., 0),
                                         Intersection::new(6., 0)]);
        let comps = w.prepare_computations(&is[0], &r, &is);
        let color = w.refracted_color(&comps, 0);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let mut w = default_world();
        let s = &mut w.objects[0];
        let mut m = s.get_material().clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        s.set_material(m);
        let x = 1./f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., x), Vector::new(0., 1., 0.));
        let is = Intersections::new(vec![Intersection::new(-x, 0),
                                         Intersection::new( x, 0)]);
        // We're inside the sphere, so take second intersection.
        let comps = w.prepare_computations(&is[1], &r, &is);
        let color = w.refracted_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_with_a_refracted_ray() {
        let mut w = default_world();
        let a = &mut w.objects[0];
        let mut m = a.get_material().clone();
        m.shader.ambient = 1.;
        m.pattern = Pattern::test();
        a.set_material(m);
        let b = &mut w.objects[1];
        let mut m = b.get_material().clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        b.set_material(m);
        let r = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
        let is = Intersections::new(vec![Intersection::new(-0.9899, 0),
                                         Intersection::new(-0.4899, 1),
                                         Intersection::new( 0.4899, 1),
                                         Intersection::new( 0.9899, 0),
                                        ]);
        let comps = w.prepare_computations(&is[2], &r, &is);
        let color = w.refracted_color(&comps, 5);
        let expected = Color::new(0., 0.99888, 0.04722);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn shade_hit_with_a_transparent_material() {
        let mut w = default_world();
        let mut floor = Shape::plane();
        let mut m = Material::new();
        m.shader.transparency = 0.5;
        m.shader.refractive_index = 1.5;
        floor.set_material(m);
        floor.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(floor);
        let mut ball = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::uniform(Color::new(1., 0., 0.));
        m.shader.ambient = 0.5;
        ball.set_material(m);
        ball.set_transform(translation(&Vector::new(0., -3.5, -0.5)));
        w.add_object(ball);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.shade_hit(&comps, 5);
        let expected = Color::new(0.93642, 0.68642, 0.68642);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn schlick_approximation_under_total_internal_reflection() {
        let mut w = World::new();
        w.add_object(glass_sphere());
        let x = 1./f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., x), Vector::new(0., 1., 0.));
        let is = Intersections::new(vec![Intersection::new(-x, 0),
                                         Intersection::new( x, 0),
                                        ]);
        let comps = w.prepare_computations(&is[1], &r, &is);
        let reflectance = schlick(&comps);
        assert_relative_eq!(reflectance, 1.);
    }

    #[test]
    fn schlick_approximation_with_perpendicular_viewing_angle() {
        let mut w = World::new();
        w.add_object(glass_sphere());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        let is = Intersections::new(vec![Intersection::new(-1., 0),
                                         Intersection::new( 1., 0),
                                        ]);
        let comps = w.prepare_computations(&is[1], &r, &is);
        let reflectance = schlick(&comps);
        assert_relative_eq!(reflectance, 0.04);
    }

    #[test]
    fn schlick_approximation_with_small_angle_and_n2_larger_than_n1() {
        let mut w = World::new();
        w.add_object(glass_sphere());
        let r = Ray::new(Point::new(0., 0.99, -2.), Vector::new(0., 0., 1.));
        let is = Intersections::new(vec![Intersection::new(1.8589, 0)]);
        let comps = w.prepare_computations(&is[0], &r, &is);
        let reflectance = schlick(&comps);
        assert_relative_eq!(reflectance, 0.48873, max_relative = 1e-4);
    }

    #[test]
    fn shade_hit_with_reflective_transparent_material() {
        let mut w = default_world();
        let mut floor = Shape::plane();
        let mut m = Material::new();
        m.shader.reflective = 0.5;
        m.shader.transparency = 0.5;
        m.shader.refractive_index = 1.5;
        floor.set_material(m);
        floor.set_transform(translation(&Vector::new(0., -1., 0.)));
        w.add_object(floor);
        let mut ball = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::uniform(Color::new(1., 0., 0.));
        m.shader.ambient = 0.5;
        ball.set_material(m);
        ball.set_transform(translation(&Vector::new(0., -3.5, -0.5)));
        w.add_object(ball);
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, 2);
        let is = Intersections::new(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.shade_hit(&comps, 5);
        let expected = Color::new(0.93391, 0.69643, 0.69243);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }
}
