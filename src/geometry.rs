use na::{Matrix4, Point3, Projective3, Rotation3, U1, U3, Unit, Vector3};

pub type Point = Point3<f64>;
pub type Vector = Vector3<f64>;

pub type Transform = Projective3<f64>;

pub fn identity() -> Transform {
    Transform::identity()
}

pub fn translation(vec: &Vector) -> Transform {
    let mut m = Matrix4::new(1., 0., 0., 0.,
                             0., 1., 0., 0.,
                             0., 0., 1., 0.,
                             0., 0., 0., 1.);
    m.fixed_slice_mut::<U3, U1>(0, 3).copy_from(vec);
    Transform::from_matrix_unchecked(m)
}

pub fn scaling(x: f64, y: f64, z: f64) -> Transform {
    let m = Matrix4::new( x, 0., 0., 0.,
                         0.,  y, 0., 0.,
                         0., 0.,  z, 0.,
                         0., 0., 0., 1.);
    Transform::from_matrix_unchecked(m)
}

pub fn rotation_x(angle: f64) -> Transform {
    rotation_around_axis(&Vector::x_axis(), angle)
}

pub fn rotation_y(angle: f64) -> Transform {
    rotation_around_axis(&Vector::y_axis(), angle)
}

pub fn rotation_z(angle: f64) -> Transform {
    rotation_around_axis(&Vector::z_axis(), angle)
}

pub fn rotation_around_axis(axis: &Unit<Vector>, angle: f64) -> Transform {
    na::convert(Rotation3::from_axis_angle(axis, angle))
}

pub fn reflect(v: &Vector, normal: &Vector) -> Vector {
    v - normal * 2. * v.dot(normal)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reflection_of_45degree_approach() {
        let v = Vector::new(1., -1., 0.);
        let n = Vector::new(0., 1., 0.);
        let r = reflect(&v, &n);
        assert_relative_eq!(r, Vector::new(1., 1., 0.));
    }

    #[test]
    fn reflecting_a_vector_off_a_slanted_surface() {
        let v = Vector::new(0., -1., 0.);
        let x = 1./f64::sqrt(2.);
        let n = Vector::new(x, x, 0.);
        let r = reflect(&v, &n);
        assert_relative_eq!(r, Vector::new(1., 0., 0.));
    }
}
