use crate::geometry::{identity, Point, Transform, Vector};
use crate::material::Material;
use crate::rays::Ray;

pub struct Sphere {
    pub transform: Transform,
    pub material: Material,
}

impl Sphere {
    pub fn new() -> Self {
        Sphere{transform: identity(), material: Material::new()}
    }

    pub fn intersect(&self, ray: &Ray) -> Vec<f64> {
        let local_ray = ray.transform(&self.transform.inverse());
        self.local_intersect(&local_ray)
    }

    pub fn local_intersect(&self, ray: &Ray) -> Vec<f64> {
        let sphere_to_ray: Vector = ray.origin - Point::new(0., 0., 0.);

        let a: f64 = ray.direction.dot(&ray.direction);
        let b: f64 = 2. * ray.direction.dot(&sphere_to_ray);
        let c: f64 = sphere_to_ray.dot(&sphere_to_ray) - 1.;

        let discriminant = b*b - 4. * a * c;

        if discriminant < 0. {
            Vec::new()
        } else {
            let t1 = (-b - discriminant.sqrt()) / (2. * a);
            let t2 = (-b + discriminant.sqrt()) / (2. * a);
            vec![t1, t2]
        }
    }

    pub fn normal_at(&self, point: &Point) -> Vector {
        let inverse_transform = self.transform.inverse();
        let local_point = inverse_transform * point;
        let local_normal = self.local_normal_at(&local_point);
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        (inv_transp * local_normal).normalize()
    }

    pub fn local_normal_at(&self, point: &Point) -> Vector {
        point - Point::new(0., 0., 0.)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_approx;
    use crate::float::Approx;

    const EPS: f64 = 1e-8;

    #[test]
    fn ray_intersecting_sphere_at_two_points() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_approx!(xs[0], 4.0, EPS);
        assert_approx!(xs[1], 6.0, EPS);
    }

    #[test]
    fn ray_intersecting_sphere_at_a_tangent() {
        let r = Ray::new(Point::new(0., 1., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_approx!(xs[0], 5.0, EPS);
        assert_approx!(xs[1], 5.0, EPS);
    }

    #[test]
    fn ray_missing_sphere() {
        let r = Ray::new(Point::new(0., 2., -5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn ray_originating_inside_sphere() {
        let r = Ray::new(Point::new(0., 0., 0.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_approx!(xs[0], -1.0, EPS);
        assert_approx!(xs[1],  1.0, EPS);
    }

    #[test]
    fn ray_originating_behind_sphere() {
        let r = Ray::new(Point::new(0., 0., 5.), Vector::new(0., 0., 1.));
        let s = Sphere::new();
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_approx!(xs[0], -6.0, EPS);
        assert_approx!(xs[1], -4.0, EPS);
    }

    #[test]
    fn default_transformation_of_a_sphere() {
        let s = Sphere::new();
        assert_approx!(s.transform, identity(), EPS);
    }

    #[test]
    fn changing_transformation_of_a_sphere() {
        let mut s = Sphere::new();
        use crate::geometry::translation;
        let m = translation(&Vector::new(2., 3., 4.));
        s.transform = m;
        assert_approx!(s.transform, m, EPS);
    }

    #[test]
    fn intersecting_scaled_sphere_with_a_ray() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Sphere::new();
        use crate::geometry::scaling;
        s.transform = scaling(2., 2., 2.);
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 2);
        assert_approx!(xs[0], 3.0, EPS);
        assert_approx!(xs[1], 7.0, EPS);
    }

    #[test]
    fn intersecting_translated_sphere_with_a_ray() {
        let r = Ray::new(Point::new(0., 0., -5.), Vector::new(0., 0., 1.));
        let mut s = Sphere::new();
        use crate::geometry::translation;
        s.transform = translation(&Vector::new(5., 0., 0.));
        let xs = s.intersect(&r);
        assert_eq!(xs.len(), 0);
    }

    #[test]
    fn normal_of_sphere_on_x_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Point::new(1., 0., 0.));
        assert_approx!(n, Vector::new(1., 0., 0.), EPS);
    }

    #[test]
    fn normal_of_sphere_on_y_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Point::new(0., 1., 0.));
        assert_approx!(n, Vector::new(0., 1., 0.), EPS);
    }

    #[test]
    fn normal_of_sphere_on_z_axis() {
        let s = Sphere::new();
        let n = s.normal_at(&Point::new(0., 0., 1.));
        assert_approx!(n, Vector::new(0., 0., 1.), EPS);
    }

    #[test]
    fn normal_of_sphere_at_nonaxial_point() {
        let s = Sphere::new();
        let x = 1./f64::sqrt(3.);
        let n = s.normal_at(&Point::new(x, x, x));
        assert_approx!(n, Vector::new(x, x, x), EPS);
    }

    #[test]
    fn normal_of_translated_sphere() {
        let mut s = Sphere::new();
        use crate::geometry::translation;
        s.transform = translation(&Vector::new(0., 1., 0.));
        let n = s.normal_at(&Point::new(0., 1.70711, -0.70711));
        assert_approx!(n, Vector::new(0., 0.70711, -0.70711), 1e-5);
    }

    #[test]
    fn normal_of_transformed_sphere() {
        let mut s = Sphere::new();
        use crate::geometry::{scaling, rotation_z};
        s.transform = scaling(1., 0.5, 1.) * rotation_z(std::f64::consts::PI/5.);
        let x = 1./f64::sqrt(2.);
        let n = s.normal_at(&Point::new(0., x, -x));
        assert_approx!(n, Vector::new(0., 0.97014, -0.24254), 1e-5);
    }

    #[test]
    fn default_sphere_material() {
        let s = Sphere::new();
        assert_eq!(s.material, Material::new());
    }

    #[test]
    fn assign_material_to_sphere() {
        let mut s = Sphere::new();
        let mut m = Material::new();
        m.ambient = 1.;
        s.material = m;
        assert_eq!(s.material, m);
    }
}
