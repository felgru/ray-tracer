use crate::geometry::Transform;
use crate::group::{Group, Groups};
use crate::object::Object;
use crate::shapes::bounds::Bounds;

pub struct ObjectStore {
    objects: Vec<Object>,
    groups: Groups,
}

impl ObjectStore {
    pub fn new() -> Self {
        ObjectStore{
            objects: Vec::new(),
            groups: Groups::new(),
        }
    }

    pub fn is_empty(&self) -> bool {
        self.objects.is_empty()
    }

    pub fn add_object<O: Into<Object>>(&mut self, obj: O) -> usize {
        self.objects.push(obj.into());
        self.objects.len() - 1
    }

    pub fn add_group(&mut self) -> (Group, usize) {
        let group = self.groups.add_group();
        let i = self.add_object(group);
        (group, i)
    }

    pub fn add_subgroup(&mut self, g: Group) -> (Group, usize) {
        let sub_group = self.groups.add_group();
        let i = self.add_object_to_group(sub_group, g);
        (sub_group, i)
    }

    pub fn add_object_to_group<O: Into<Object>>(&mut self, obj: O,
                                                group: Group) -> usize {
        let i = self.add_object(obj.into());
        self.groups.add_child_to_group(i, group, &mut self.objects);
        i
    }

    pub fn get_object(&self, i: usize) -> &Object {
        &self.objects[i]
    }

    #[cfg(test)]
    pub fn get_object_mut(&mut self, i: usize) -> &mut Object {
        &mut self.objects[i]
    }

    pub fn get_group(&self, i: usize) -> Group {
        assert!(i < self.groups.len());
        Group::with_index(i)
    }

    pub fn groups(&self) -> &Groups {
        &self.groups
    }

    pub fn get_transform_of_group(&self, group: Group) -> &Transform {
        self.groups.transform(group)
    }

    pub fn get_bounds_of_group(&self, group: Group) -> &Bounds {
        self.groups.bounds(group)
    }

    pub fn children_of_group(&self, group: Group) -> &[usize] {
        &self.groups.children(group)
    }

    pub fn get_transform_of_object(&self, i: usize) -> &Transform {
        self.get_object(i).get_transform(&self.groups)
    }

    pub fn set_transform_of_object(&mut self, i: usize, transform: Transform) {
        let obj = &mut self.objects[i];
        obj.set_transform(transform, &mut self.groups);
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::material::Material;
    use crate::patterns::Pattern;
    use crate::shapes::Shape;

    #[test]
    fn create_an_object_store() {
        let obs = ObjectStore::new();
        assert!(obs.objects.is_empty());
        assert!(obs.groups.is_empty());
    }

    #[test]
    fn one_object() {
        let mut obs = ObjectStore::new();
        use crate::shapes::Shape;
        let mut s = Shape::sphere();
        let mut m = Material::new();
        m.pattern = Pattern::test();
        s.set_material(m);
        let i = obs.add_object(s);
        assert_eq!(i, 0);
        assert_eq!(obs.objects.len(), 1);
        assert!(obs.groups.is_empty());
        assert_eq!(obs.get_object(i).get_material().pattern,
                   Pattern::test());
    }

    #[test]
    fn adding_a_child_to_a_group() {
        let mut objects = ObjectStore::new();
        let (g, _) = objects.add_group();
        assert_eq!(g.index, 0);
        let i = objects.add_object_to_group(Shape::sphere(), g);
        assert!(!objects.children_of_group(g).is_empty());
        assert_eq!(i, 1);
        let i = objects.children_of_group(g)[0];
        assert_eq!(i, 1);
        assert_eq!(objects.get_object(i).parent.index, 0);
    }
}
