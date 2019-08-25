use crate::color::Color;
use crate::geometry::{identity, Point, Transform};

#[derive(Debug, PartialEq, Copy, Clone)]
pub struct Pattern {
    data: PatternData,
    pub transform: Transform,
}

impl Pattern {
    pub fn uniform(color: Color) -> Self {
        Self::with_data(PatternData::PlainColor(PlainColor{color}))
    }

    pub fn stripes(a: Color, b: Color) -> Self {
        Self::with_data(PatternData::StripePattern(StripePattern{a, b}))
    }

    fn with_data(data: PatternData) -> Self {
        Pattern{data, transform: identity()}
    }

    /// Give color at point in Pattern space
    pub fn local_color_at(&self, point: &Point) -> Color {
        self.get_pattern_data().local_color_at(point)
    }

    pub fn set_transform(&mut self, transform: Transform) {
        self.transform = transform;
    }

    fn get_pattern_data<'a>(&'a self) -> &'a dyn PatternType {
        use PatternData::*;
        match &self.data {
            PlainColor(p) => p,
            StripePattern(p) => p,
        }
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
enum PatternData {
    PlainColor(PlainColor),
    StripePattern(StripePattern),
}

trait PatternType {
    fn local_color_at(&self, point: &Point) -> Color;
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct PlainColor {
    color: Color,
}

impl PatternType for PlainColor {
    fn local_color_at(&self, _point: &Point) -> Color {
        self.color
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct StripePattern {
    a: Color,
    b: Color,
}

impl PatternType for StripePattern {
    fn local_color_at(&self, point: &Point) -> Color {
        if point[0].floor() as i64 % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stripe_pattern_constant_in_y() {
        let pattern = Pattern::stripes(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 1., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 2., 0.)), Color::white());
    }

    #[test]
    fn stripe_pattern_constant_in_z() {
        let pattern = Pattern::stripes(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 1.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 2.)), Color::white());
    }

    #[test]
    fn stripe_pattern_alternates_in_x() {
        let pattern = Pattern::stripes(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0.9, 0., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(1., 0., 0.)), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(-0.1, 0., 0.)), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(-1., 0., 0.)), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(-1.1, 0., 0.)), Color::white());
    }
}
