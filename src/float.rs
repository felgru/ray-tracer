pub trait Approx<Precision=f64> {
    fn approx(self, rhs: Self, eps: Precision) -> bool;
}

impl Approx for f64 {
    fn approx(self, rhs: Self, eps: f64) -> bool {
        (self - rhs).abs() < eps
    }
}

// based on assert_eq
#[macro_export]
macro_rules! assert_approx {
    ($left:expr, $right:expr, $eps: expr) => ({
        match (&$left, &$right, &$eps) {
            (left_val, right_val, eps_val) => {
                if !left_val.approx(*right_val, *eps_val) {
                    // The reborrows below are intentional. Without them, the stack slot for the
                    // borrow is initialized even before the values are compared, leading to a
                    // noticeable slow down.
                    panic!(r#"assertion failed: `left.approx(right, eps)`
  left: `{:?}`,
 right: `{:?}`,
   eps: `{:?}`"#, &*left_val, &*right_val, &*eps_val)
                }
            }
        }
    });
    ($left:expr, $right:expr, $eps: expr,) => ({
        $crate::assert_approx!($left, $right, $eps)
    });
    ($left:expr, $right:expr, $eps: expr, $($arg:tt)+) => ({
        match (&($left), &($right), &($eps)) {
            (left_val, right_val, eps_val) => {
                if !left_val.approx(*right_val, *eps_val) {
                    // The reborrows below are intentional. Without them, the stack slot for the
                    // borrow is initialized even before the values are compared, leading to a
                    // noticeable slow down.
                    panic!(r#"assertion failed: `left.approx(right, eps)`
  left: `{:?}`,
 right: `{:?}`,
   eps: `{:?}`: {}"#, &*left_val, &*right_val, &*eps_val,
                           format_args!($($arg)+))
                }
            }
        }
    });
}
