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

    pub fn rings(a: Color, b: Color) -> Self {
        Self::with_data(PatternData::RingPattern(RingPattern{a, b}))
    }

    pub fn checkers(a: Color, b: Color) -> Self {
        Self::with_data(PatternData::CheckerPattern(CheckerPattern{a, b}))
    }

    pub fn gradient(a: Color, b: Color) -> Self {
        Self::with_data(PatternData::GradientPattern(GradientPattern{a, b}))
    }

    pub fn test() -> Self {
        Self::with_data(PatternData::TestPattern(TestPattern{}))
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
            RingPattern(p) => p,
            CheckerPattern(p) => p,
            GradientPattern(p) => p,
            TestPattern(p) => p,
        }
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
enum PatternData {
    PlainColor(PlainColor),
    StripePattern(StripePattern),
    CheckerPattern(CheckerPattern),
    RingPattern(RingPattern),
    GradientPattern(GradientPattern),
    TestPattern(TestPattern),
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

#[derive(Debug, PartialEq, Copy, Clone)]
struct RingPattern {
    a: Color,
    b: Color,
}

impl PatternType for RingPattern {
    fn local_color_at(&self, point: &Point) -> Color {
        let (x, z) = (point[0], point[2]);
        if (x*x + z*z).sqrt().floor() as i64 % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}
#[derive(Debug, PartialEq, Copy, Clone)]
struct CheckerPattern {
    a: Color,
    b: Color,
}

impl PatternType for CheckerPattern {
    fn local_color_at(&self, point: &Point) -> Color {
        if point.iter().map(|x| x.floor() as i64).sum::<i64>() % 2 == 0 {
            self.a
        } else {
            self.b
        }
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct GradientPattern {
    a: Color,
    b: Color,
}

impl PatternType for GradientPattern {
    fn local_color_at(&self, point: &Point) -> Color {
        let distance = self.b - self.a;
        let fraction = point[0] - point[0].floor();

        self.a + distance * fraction
    }
}

#[derive(Debug, PartialEq, Copy, Clone)]
struct TestPattern {}

impl PatternType for TestPattern {
    fn local_color_at(&self, point: &Point) -> Color {
        Color::new(point[0], point[1], point[2])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn default_pattern_transform() {
        let pattern = Pattern::uniform(Color::white());
        assert_eq!(pattern.transform, identity());
    }

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

    #[test]
    fn ring_pattern_alternates_in_x_and_z() {
        let pattern = Pattern::rings(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)), Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(1., 0., 0.)), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 1.)), Color::black());
        let x = 0.708; // just slightly more than 1/√2
        assert_eq!(pattern.local_color_at(&Point::new(x, 0., x)), Color::black());
    }

    #[test]
    fn checker_pattern_alternates_in_x() {
        let pattern = Pattern::checkers(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0.99, 0., 0.)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(1.01, 0., 0.)),
                                          Color::black());
    }

    #[test]
    fn checker_pattern_alternates_in_y() {
        let pattern = Pattern::checkers(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0.99, 0.)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 1.01, 0.)),
                                          Color::black());
    }

    #[test]
    fn checker_pattern_alternates_in_z() {
        let pattern = Pattern::checkers(Color::white(), Color::black());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 0.99)),
                                          Color::white());
        assert_eq!(pattern.local_color_at(&Point::new(0., 0., 1.01)),
                                          Color::black());
    }

    #[test]
    fn gradient_linearly_interpolates() {
        let pattern = Pattern::gradient(Color::white(), Color::black());
        assert_relative_eq!(pattern.local_color_at(&Point::new(0., 0., 0.)),
                            Color::white());
        assert_relative_eq!(pattern.local_color_at(&Point::new(0.25, 0., 0.)),
                            Color::white() * 0.75);
        assert_relative_eq!(pattern.local_color_at(&Point::new(0.5, 0., 0.)),
                            Color::white() * 0.5);
        assert_relative_eq!(pattern.local_color_at(&Point::new(0.75, 0., 0.)),
                            Color::white() * 0.25);
    }
}
