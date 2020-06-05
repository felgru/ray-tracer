use crate::color::Color;
use crate::geometry::{Point, reflect, Transform, Vector};
use crate::group::GroupIndex;
use crate::intersections::{Intersection, Intersections};
use crate::light::PointLight;
use crate::rays::Ray;
use crate::object_store::{ObjectIndex, ObjectStore};
use crate::shapes::{Shape, ShapeIndex};

pub struct World {
    scene: Vec<ObjectIndex>,
    objects: ObjectStore,
    lights: Vec<PointLight>,
}

impl World {
    pub fn new() -> Self {
        World{
            scene: Vec::new(),
            objects: ObjectStore::new(),
            lights: Vec::new(),
        }
    }

    pub fn add_shape(&mut self, shape: Shape, transform: Transform)
                                                    -> ShapeIndex {
        let i = self.objects.add_shape(shape, transform);
        self.scene.push(ObjectIndex::Shape(i));
        i
    }

    pub fn add_group(&mut self, transform: Transform) -> GroupIndex {
        let group = self.objects.add_group(transform);
        self.scene.push(ObjectIndex::Group(group));
        group
    }

    pub fn add_subgroup(&mut self, transform: Transform, g: GroupIndex)
                                                        -> GroupIndex {
        self.objects.add_subgroup(transform, g)
    }

    pub fn add_shape_to_group(&mut self, shape: Shape, transform: Transform,
                              group: GroupIndex) -> ShapeIndex {
        self.objects.add_shape_to_group(shape, transform, group)
    }

    pub fn add_light(&mut self, light: PointLight) -> usize {
        self.lights.push(light);
        self.lights.len() - 1
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
        let mut intersections = Intersections::new();
        for i in self.scene.iter().copied() {
            let is = self.objects.intersect(i, ray);
            intersections.add_intersections_from(is);
        }
        intersections
    }

    fn prepare_computations(&self, hit: &Intersection,
                            ray: &Ray,
                            intersections: &Intersections)
                                        -> PreparedComputations {
        PreparedComputations::new(self, hit, ray, intersections)
    }

    fn shade_hit(&self, comps: &PreparedComputations,
                 remaining: usize) -> Color {
        let lighting = |l| self.objects
                               .lighting(comps.shape_index,
                                         l,
                                         &comps.point,
                                         &comps.eyev,
                                         &comps.normalv,
                                         self.is_shadowed(&comps.over_point, l));
        let surface = self.lights.iter().map(lighting).sum::<Color>();

        let reflected = self.reflected_color(comps, remaining);
        let refracted = self.refracted_color(comps, remaining);

        let shader = self.objects.get_material_of(comps.shape_index).shader;
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
        let reflectivity = self.objects.get_material_of(comps.shape_index)
                                       .shader.reflective;
        if remaining == 0 || reflectivity == 0. {
            return Color::black();
        }

        let reflect_ray = Ray::new(comps.over_point, comps.reflectv);
        let color = self.color_at(&reflect_ray, remaining - 1);

        color * reflectivity
    }

    fn refracted_color(&self, comps: &PreparedComputations,
                       remaining: usize) -> Color {
        if remaining == 0
           || self.objects.get_material_of(comps.shape_index)
                          .shader.transparency == 0. {
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
            * self.objects.get_material_of(comps.shape_index).shader
                                                             .transparency
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

    pub fn set_transform_of(&mut self, i: ObjectIndex, t: Transform) {
        self.objects.set_transform_of_object(i, t);
    }
}

struct PreparedComputations {
    pub t: f64,
    pub shape_index: ShapeIndex,
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

impl PreparedComputations {
    pub fn new(world: &World, hit: &Intersection, ray: &Ray,
               intersections: &Intersections) -> Self {
        let point = ray.position(hit.t);
        let shape_index = hit.shape_index;
        let eyev = -ray.direction;
        let mut normalv = world.objects.normal_at(shape_index, &point);
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
            shape_index,
            point,
            over_point,
            under_point,
            eyev,
            normalv,
            reflectv: reflect(&ray.direction, &normalv),
            n1,
            n2,
            inside,
        }
    }

    fn get_refractive_indices_at(world: &World, hit: &Intersection,
                                 xs: &Intersections) -> (f64, f64) {
        let mut visited_objects = Vec::<ShapeIndex>::new();

        for i in xs.iter().take_while(|&i| i != hit) {
            let idx = i.shape_index;
            if let Some(pos) = visited_objects.iter().position(|p| *p == idx) {
                visited_objects.remove(pos);
            } else {
                visited_objects.push(idx);
            }
        }

        let refractive_index_of = |i| {
            world.objects.get_material_of(i).shader.refractive_index
        };

        let (n1, n2);

        let hit_idx = hit.shape_index;
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

    let s2 = Shape::sphere();

    let mut w = World::new();
    w.add_light(light);
    use crate::geometry::identity;
    w.add_shape(s1, identity());
    use crate::geometry::scaling;
    w.add_shape(s2, scaling(0.5, 0.5, 0.5));
    w
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::material::Material;
    use crate::object_store::Parent;
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
        assert_eq!(w.objects.get_material_of(ShapeIndex::new(0)).pattern,
                   Pattern::uniform(c));
        let s1_idx = ObjectIndex::Shape(ShapeIndex::new(1));
        assert_eq!(w.objects.get_transform_of_object(s1_idx), &t);
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
        let indx = w.add_shape(shape, identity());
        let i = Intersection::new(4., indx);
        let is = Intersections::from_vec(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(w.objects.get_material_of(comps.shape_index).pattern,
                   Pattern::uniform(red));
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
        let indx = w.add_shape(shape, identity());
        let i = Intersection::new(1., indx);
        let is = Intersections::from_vec(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.t, i.t);
        assert_eq!(w.objects.get_material_of(comps.shape_index).pattern,
                   Pattern::uniform(red));
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
        let i = Intersection::new(4., ShapeIndex::new(0));
        let is = Intersections::from_vec(vec![i]);
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
        let i = Intersection::new(0.5, ShapeIndex::new(1));
        let is = Intersections::from_vec(vec![i]);
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
            set_ambient_of_shape(&mut w, 0, 1.);
            set_ambient_of_shape(&mut w, 1, 1.);
            w
        };
        let r = Ray::new(Point::new(0., 0., 0.75), Vector::new(0., 0., -1.));
        let c = w.color_at(&r, 4);
        assert_relative_eq!(c, Color::white());
    }

    fn set_ambient_of_shape(w: &mut World, i: usize, ambient: f64) {
        let si = ShapeIndex::new(i);
        let mut m = w.objects.get_material_of(si).clone();
        m.shader.ambient = ambient;
        w.objects.set_material_of(si, m);
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
        w.add_shape(Shape::sphere(), identity());
        w.add_shape(Shape::sphere(), translation(&Vector::new(0., 0., 10.)));
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(4., ShapeIndex::new(1));
        let is = Intersections::from_vec(vec![i]);
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
        w.add_shape(Shape::sphere(), translation(&Vector::new(0., 0., 1.)));
        w
    }

    fn offset_test_computations(w: &World) -> PreparedComputations {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let i = Intersection::new(5., ShapeIndex::new(0));
        let is = Intersections::from_vec(vec![i]);
        w.prepare_computations(&i, &r, &is)
    }

    #[test]
    fn precomputing_the_reflection_vector() {
        let mut w = World::new();
        let s = w.add_shape(Shape::plane(), identity());
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 1., -1.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, s);
        let is = Intersections::from_vec(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        assert_relative_eq!(comps.reflectv,
                            Vector::new(0., 1./sqrt2, 1./sqrt2));
    }

    #[test]
    fn reflected_color_of_a_nonreflective_material() {
        let mut w = default_world();
        set_ambient_of_shape(&mut w, 1, 1.);
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let i = Intersection::new(1., ShapeIndex::new(1));
        let is = Intersections::from_vec(vec![i]);
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
        w.add_shape(shape, translation(&Vector::new(0., -1., 0.)));
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, ShapeIndex::new(2));
        let is = Intersections::from_vec(vec![i]);
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
        w.add_shape(shape, translation(&Vector::new(0., -1., 0.)));
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, ShapeIndex::new(2));
        let is = Intersections::from_vec(vec![i]);
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
        w.add_shape(shape, translation(&Vector::new(0., -1., 0.)));
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, ShapeIndex::new(2));
        let is = Intersections::from_vec(vec![i]);
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
        w.add_shape(lower, translation(&Vector::new(0., -1., 0.)));
        let mut upper = Shape::plane();
        upper.set_material(m.clone());
        w.add_shape(upper, translation(&Vector::new(0., 1., 0.)));
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

    fn add_refraction_test_sphere(w: &mut World, refractive_index: f64,
                                  transform: Transform) -> ShapeIndex {
        let mut s = glass_sphere();
        let mut m = *s.get_material();
        m.shader.refractive_index = refractive_index;
        s.set_material(m);
        w.add_shape(s, transform)
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
        let mut w = World::new();
        add_refraction_test_sphere(&mut w, 1.5, scaling(2., 2., 2.));
        add_refraction_test_sphere(&mut w, 2.0,
                translation(&Vector::new(0., 0., -0.25)));
        add_refraction_test_sphere(&mut w, 2.5,
                translation(&Vector::new(0., 0., 0.25)));
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
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(4., s),
                                              Intersection::new(6., s)]);
        let comps = w.prepare_computations(&is[0], &r, &is);
        let color = w.refracted_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_at_maximal_recursion_depth() {
        let mut w = default_world();
        let si = ShapeIndex::new(0);
        let mut m = w.objects.get_material_of(si).clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        w.objects.set_material_of(si, m);
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(4., s),
                                              Intersection::new(6., s)]);
        let comps = w.prepare_computations(&is[0], &r, &is);
        let color = w.refracted_color(&comps, 0);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_under_total_internal_reflection() {
        let mut w = default_world();
        let si = ShapeIndex::new(0);
        let mut m = w.objects.get_material_of(si).clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        w.objects.set_material_of(si, m);
        let x = 1./f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., x), Vector::new(0., 1., 0.));
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(-x, s),
                                              Intersection::new( x, s)]);
        // We're inside the sphere, so take second intersection.
        let comps = w.prepare_computations(&is[1], &r, &is);
        let color = w.refracted_color(&comps, 4);
        assert_eq!(color, Color::black());
    }

    #[test]
    fn refracted_color_with_a_refracted_ray() {
        let mut w = default_world();
        let ai = ShapeIndex::new(0);
        let mut m = w.objects.get_material_of(ai).clone();
        m.shader.ambient = 1.;
        m.pattern = Pattern::test();
        w.objects.set_material_of(ai, m);
        let bi = ShapeIndex::new(1);
        let mut m = w.objects.get_material_of(bi).clone();
        m.shader.transparency = 1.;
        m.shader.refractive_index = 1.5;
        w.objects.set_material_of(bi, m);
        let r = Ray::new(Point::new(0., 0., 0.1), Vector::new(0., 1., 0.));
        let s0 = ShapeIndex::new(0);
        let s1 = ShapeIndex::new(1);
        let is = Intersections::from_vec(vec![Intersection::new(-0.9899, s0),
                                              Intersection::new(-0.4899, s1),
                                              Intersection::new( 0.4899, s1),
                                              Intersection::new( 0.9899, s0),
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
        let floor = w.add_shape(floor, translation(&Vector::new(0., -1., 0.)));
        let mut ball = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::uniform(Color::new(1., 0., 0.));
        m.shader.ambient = 0.5;
        ball.set_material(m);
        w.add_shape(ball, translation(&Vector::new(0., -3.5, -0.5)));
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, floor);
        let is = Intersections::from_vec(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.shade_hit(&comps, 5);
        let expected = Color::new(0.93642, 0.68642, 0.68642);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn schlick_approximation_under_total_internal_reflection() {
        let mut w = World::new();
        w.add_shape(glass_sphere(), identity());
        let x = 1./f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., x), Vector::new(0., 1., 0.));
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(-x, s),
                                              Intersection::new( x, s),
                                             ]);
        let comps = w.prepare_computations(&is[1], &r, &is);
        let reflectance = schlick(&comps);
        assert_relative_eq!(reflectance, 1.);
    }

    #[test]
    fn schlick_approximation_with_perpendicular_viewing_angle() {
        let mut w = World::new();
        w.add_shape(glass_sphere(), identity());
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 1., 0.));
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(-1., s),
                                              Intersection::new( 1., s),
                                             ]);
        let comps = w.prepare_computations(&is[1], &r, &is);
        let reflectance = schlick(&comps);
        assert_relative_eq!(reflectance, 0.04);
    }

    #[test]
    fn schlick_approximation_with_small_angle_and_n2_larger_than_n1() {
        let mut w = World::new();
        w.add_shape(glass_sphere(), identity());
        let r = Ray::new(Point::new(0., 0.99, -2.), Vector::new(0., 0., 1.));
        let s = ShapeIndex::new(0);
        let is = Intersections::from_vec(vec![Intersection::new(1.8589, s)]);
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
        let floor = w.add_shape(floor, translation(&Vector::new(0., -1., 0.)));
        let mut ball = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::uniform(Color::new(1., 0., 0.));
        m.shader.ambient = 0.5;
        ball.set_material(m);
        w.add_shape(ball, translation(&Vector::new(0., -3.5, -0.5)));
        let sqrt2 = f64::sqrt(2.);
        let r = Ray::new(Point::new(0., 0., -3.),
                         Vector::new(0., -1./sqrt2, 1./sqrt2));
        let i = Intersection::new(sqrt2, floor);
        let is = Intersections::from_vec(vec![i]);
        let comps = w.prepare_computations(&i, &r, &is);
        let color = w.shade_hit(&comps, 5);
        let expected = Color::new(0.93391, 0.69643, 0.69243);
        assert_relative_eq!(color, expected, max_relative = 1e-4);
    }

    #[test]
    fn parents_in_group_hierarchy() {
        let mut world = World::new();
        let g1 = world.add_group(rotation_y(std::f64::consts::PI/2.));
        let g2 = world.add_subgroup(scaling(2., 2., 2.), g1);
        let si = world.add_shape_to_group(Shape::sphere(),
                            translation(&Vector::new(5., 0., 0.)), g2);
        assert_eq!(world.objects.parent_of_object(si), Parent::Group(g2));
        assert_eq!(world.objects.parent_of_object(g2), Parent::Group(g1));
        assert_eq!(world.objects.parent_of_object(g1), Parent::None);
    }

    #[test]
    fn converting_a_point_from_world_to_object_space() {
        let mut world = World::new();
        let g1 = world.add_group(rotation_y(std::f64::consts::PI/2.));
        let g2 = world.add_subgroup(scaling(2., 2., 2.), g1);
        let si = world.add_shape_to_group(Shape::sphere(),
                        translation(&Vector::new(5., 0., 0.)), g2);
        let p = world.objects.world_to_object(si.into(),
                                              &Point::new(-2., 0., -10.));
        assert_relative_eq!(p, Point::new(0., 0., -1.));
    }

    #[test]
    fn converting_a_normal_from_object_to_world_space() {
        let mut world = World::new();
        let g1 = world.add_group(rotation_y(std::f64::consts::PI/2.));
        let g2 = world.add_subgroup(scaling(1., 2., 3.), g1);
        let si = world.add_shape_to_group(Shape::sphere(),
                        translation(&Vector::new(5., 0., 0.)), g2);
        let x = f64::sqrt(3.) / 3.;
        let p = world.objects.normal_to_world(si.into(), &Vector::new(x, x, x));
        assert_relative_eq!(p, Vector::new(0.2857, 0.4286, -0.8571),
                            max_relative = 1e-4);
    }

    #[test]
    fn finding_the_normal_on_a_child_object() {
        let mut world = World::new();
        let g1 = world.add_group(rotation_y(std::f64::consts::PI/2.));
        let g2 = world.add_subgroup(scaling(1., 2., 3.), g1);
        let si = world.add_shape_to_group(Shape::sphere(),
                        translation(&Vector::new(5., 0., 0.)), g2);
        let sqrt3 = f64::sqrt(3.);
        let p = world.objects.normal_at(si,
                    &Point::new(sqrt3, 2./3.*sqrt3, -5. - sqrt3/3.));
        assert_relative_eq!(p, Vector::new(0.2857, 0.4286, -0.8571),
                            max_relative = 1e-4);
    }

    #[test]
    fn object_gets_transformed_with_group() {
        let mut world = World::new();
        use crate::geometry::translation;
        let group = world.add_group(translation(&Vector::new(0., 2., 0.)));
        use crate::shapes::Shape;
        world.add_shape_to_group(Shape::sphere(), identity(), group);

        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let xs = world.intersect(&r);
        // since the transformation of the group moved the sphere above the ray
        // we do not get any intersections.
        assert!(xs.is_empty());
    }
}
