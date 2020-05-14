use crate::color::Color;
use crate::geometry::{Point, Transform, Vector};
use crate::group::{Group, Groups};
use crate::intersections::Intersections;
use crate::light::PointLight;
use crate::material::Material;
use crate::rays::Ray;
use crate::shapes::bounds::Bounds;
use crate::shapes::Shape;

pub struct Object {
    pub parent: Group,
    pub object: ObjectData,
}

impl Object {
    pub fn normal_at(&self, world_point: &Point, groups: &Groups) -> Vector {
        if let Some(shp) = self.as_shape() {
            let local_point = self.world_to_object(world_point, groups);
            let local_normal = shp.local_normal_at(&local_point);
            self.normal_to_world(&local_normal, groups)
        } else {
            panic!("Trying to take normal of non-shape object!");
        }
    }

    pub fn color_at(&self, world_point: &Point, groups: &Groups) -> Color {
        let pattern = match self.as_shape() {
            Some(shp) => shp.get_material().pattern,
            None => panic!("Trying to take color of non-shape object!"),
        };
        let object_point = self.world_to_object(world_point, groups);
        let pattern_point = pattern.transform.inverse() * object_point;
        pattern.local_color_at(&pattern_point)
    }

    pub fn lighting(&self, light: &PointLight,
                    position: &Point, eyev: &Vector, normalv: &Vector,
                    in_shadow: bool,
                    groups: &Groups) -> Color {
        let color = self.color_at(position, groups);
        let shader = self.get_material().shader;
        shader.lighting(color, light, position, eyev, normalv, in_shadow)
    }

    pub fn intersect(&self, ray: &Ray, obj_id: usize,
                     objects: &[Object], groups: &Groups) -> Intersections {
        use ObjectData::*;
        match &self.object {
            Shape(shp) => {
                let mut intersections = Intersections::new();
                intersections.add_intersections(obj_id, shp.intersect(ray));
                intersections
            }
            Group(grp) => grp.intersect(ray, objects, groups),
        }
    }

    pub fn local_bounds(&self, groups: &Groups) -> Bounds {
        use ObjectData::*;
        match &self.object {
            Shape(shp) => shp.local_bounds(),
            Group(grp) => grp.local_bounds(groups),
        }
    }

    pub fn parent(&self) -> Option<Group> {
        if self.parent.index != std::usize::MAX {
            Some(self.parent)
        } else {
            None
        }
    }

    pub fn get_transform<'a>(&'a self, groups: &'a Groups) -> &'a Transform {
        self.object.get_transform(groups)
    }

    pub fn set_transform(&mut self, transform: Transform,
                         groups: &mut Groups) {
        self.object.set_transform(transform, groups);
    }

    pub fn get_material(&self) -> &Material {
        if let Some(shp) = self.as_shape() {
            shp.get_material()
        } else {
            panic!("Trying to take material of non-shape object!")
        }
    }

    pub fn set_material(&mut self, material: Material) {
        match &mut self.object {
            ObjectData::Shape(shp) => shp.set_material(material),
            _ => panic!("Trying to set material of non-shape object!"),
        }
    }

    pub fn is_group(&self) -> bool {
        match &self.object {
            ObjectData::Group(_) => true,
            _ => false,
        }
    }

    fn as_shape(&self) -> Option<&Shape> {
        match &self.object {
            ObjectData::Shape(shp) => Some(&shp),
            _ => None,
        }
    }

    pub fn world_to_object(&self, p: &Point, groups: &Groups) -> Point {
        let p = match self.parent() {
            Some(parent) => parent.world_to_object(&p, groups),
            None => *p,
        };
        self.get_transform(groups).inverse() * p
    }

    pub fn normal_to_world(&self, n: &Vector, groups: &Groups) -> Vector {
        let inverse_transform = self.get_transform(groups).inverse();
        use na::{U3};
        let inv_transp = inverse_transform.to_homogeneous()
                                          .fixed_slice::<U3,U3>(0, 0)
                                          .transpose().to_homogeneous();
        let inv_transp = Transform::from_matrix_unchecked(inv_transp);
        let world_normal = (inv_transp * n).normalize();
        match self.parent() {
            Some(parent) => parent.normal_to_world(&world_normal, groups),
            None => world_normal,
        }
    }
}

impl From<Shape> for Object {
    fn from(s: Shape) -> Self {
        Object {
            parent: Group::with_index(std::usize::MAX),
            object: ObjectData::Shape(s),
        }
    }
}

impl From<Group> for Object {
    fn from(g: Group) -> Self {
        Object {
            parent: Group::with_index(std::usize::MAX),
            object: ObjectData::Group(g),
        }
    }
}

pub enum ObjectData {
    Group(Group),
    Shape(Shape),
}

impl ObjectData {
    pub fn get_transform<'a>(&'a self, groups: &'a Groups) -> &'a Transform {
        use ObjectData::*;
        match self {
            Group(grp) => groups.transform(*grp),
            Shape(shp) => shp.get_transform(),
        }
    }

    pub fn set_transform(&mut self, transform: Transform, groups: &mut Groups) {
        use ObjectData::*;
        match self {
            Group(grp) => *groups.transform_mut(*grp) = transform,
            Shape(shp) => shp.set_transform(transform),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    use crate::geometry::*;
    use crate::patterns::Pattern;
    use crate::shapes::Shape;

    #[test]
    fn normal_of_translated_sphere() {
        let mut s: Object = Shape::sphere().into();
        use crate::geometry::translation;
        s.set_transform(translation(&Vector::new(0., 1., 0.)),
                        &mut Groups::new());
        let n = s.normal_at(&Point::new(0., 1.70711, -0.70711),
                            &Groups::new());
        assert_relative_eq!(n, Vector::new(0., 0.70711, -0.70711),
                            max_relative = 1e-5);
    }

    #[test]
    fn normal_of_transformed_sphere() {
        let mut s: Object = Shape::sphere().into();
        use crate::geometry::{scaling, rotation_z};
        s.set_transform(scaling(1., 0.5, 1.)
                        * rotation_z(std::f64::consts::PI/5.),
                        &mut Groups::new());
        let x = 1./f64::sqrt(2.);
        let n = s.normal_at(&Point::new(0., x, -x), &Groups::new());
        assert_relative_eq!(n, Vector::new(0., 0.97014, -0.24254),
                            max_relative = 1e-4);
    }

    #[test]
    fn stripes_with_transformed_sphere() {
        let mut object: Object = Shape::sphere().into();
        object.set_transform(scaling(2., 2., 2.), &mut Groups::new());
        let pattern = Pattern::stripes(Color::white(), Color::black());
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.), &Groups::new());
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern() {
        let mut object: Object = Shape::sphere().into();
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = scaling(2., 2., 2.);
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(1.5, 0., 0.), &Groups::new());
        assert_eq!(c, Color::white());
    }

    #[test]
    fn stripes_with_transformed_pattern_and_object() {
        let mut object: Object = Shape::sphere().into();
        object.set_transform(scaling(2., 2., 2.), &mut Groups::new());
        let mut pattern = Pattern::stripes(Color::white(), Color::black());
        pattern.transform = translation(&Vector::new(0.5, 0., 0.));
        set_pattern_of_object(&mut object, pattern);
        let c = object.color_at(&Point::new(2.5, 0., 0.), &Groups::new());
        assert_eq!(c, Color::white());
    }

    fn set_pattern_of_object(object: &mut Object, pattern: Pattern) {
        let mut material = object.get_material().clone();
        material.pattern = pattern;
        object.set_material(material);
    }

    #[test]
    fn lighting_with_stripe_pattern() {
        use crate::patterns::Pattern;
        let mut m = Material::new();
        m.pattern = Pattern::stripes(Color::white(), Color::black());
        m.shader.ambient = 1.;
        m.shader.diffuse = 0.;
        m.shader.specular = 0.;
        let mut obj: Object = Shape::sphere().into();
        obj.set_material(m);
        let eyev = Vector::new(0., 0., -1.);
        let normalv = Vector::new(0., 0., -1.);
        let light = PointLight::new(Point::new(0., 0., -10.), Color::white());
        let c1 = obj.lighting(&light, &Point::new(0.9, 0., 0.),
                              &eyev, &normalv, false,
                              &Groups::new());
        let c2 = obj.lighting(&light, &Point::new(1.1, 0., 0.),
                              &eyev, &normalv, false,
                              &Groups::new());
        assert_relative_eq!(c1, Color::white());
        assert_relative_eq!(c2, Color::black());
    }
}
