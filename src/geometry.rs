// SPDX-FileCopyrightText: 2019–2021 Felix Gruber
//
// SPDX-License-Identifier: GPL-3.0-or-later

use na::{Matrix4, Point3, Projective3, Rotation3, Unit, Vector3};

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
    m.fixed_slice_mut::<3, 1>(0, 3).copy_from(vec);
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

pub fn view_transform(from: &Point, to: &Point, up: &Vector) -> Transform {
    let forward = (to - from).normalize();
    let upn = up.normalize();
    let left = forward.cross(&upn);
    let true_up = left.cross(&forward);

    let orientation = Matrix4::new(
        left[0],     left[1],     left[2],     0.,
        true_up[0],  true_up[1],  true_up[2],  0.,
        -forward[0], -forward[1], -forward[2], 0.,
        0.,          0.,          0.,          1.);
    let orientation = Transform::from_matrix_unchecked(orientation);

    orientation * translation(&(Point::new(0., 0., 0.) - from))
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

    #[test]
    fn view_transform_for_default_orientation() {
        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., -1.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(&from, &to, &up);
        assert_relative_eq!(t, identity());
    }

    #[test]
    fn view_transform_looking_in_positive_z_direction() {
        let from = Point::new(0., 0., 0.);
        let to = Point::new(0., 0., 1.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(&from, &to, &up);
        assert_relative_eq!(t, scaling(-1., 1., -1.));
    }

    #[test]
    fn view_transform_moves_the_world() {
        let from = Point::new(0., 0., 8.);
        let to = Point::new(0., 0., 0.);
        let up = Vector::new(0., 1., 0.);
        let t = view_transform(&from, &to, &up);
        assert_relative_eq!(t, translation(&Vector::new(0., 0., -8.)));
    }
    #[test]
    fn arbitrary_view_transformation() {
        let from = Point::new(1., 3., 2.);
        let to = Point::new(4., -2., 8.);
        let up = Vector::new(1., 1., 0.);
        let t = view_transform(&from, &to, &up);
        let expected = Transform::from_matrix_unchecked(Matrix4::<f64>::new(
            -0.50709, 0.50709,  0.67612, -2.36643,
             0.76772, 0.60609,  0.12122, -2.82843,
            -0.35857, 0.59761, -0.71714,  0.00000,
             0.00000, 0.00000,  0.00000,  1.00000,
        ));
        assert_relative_eq!(t, expected, max_relative = 1e-4);
    }
}
