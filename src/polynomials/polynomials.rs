const N: usize = 4;   // in real 256
const Q: u32 = 7681;    // in real 3329
const W: u32 = 3383;

const NTT: [[u32; N]; N] = [[1, 1   , 1   , 1   ],
                            [1, 3383, 7680, 4298],
                            [1, 7680, 1   , 7680],
                            [1, 4298, 7680, 3383]]; // well its bigger in reality 256 long.



#[derive(Clone, Debug)]
struct Poly {
    coeffs: [u32; N],
}

#[derive(Clone, Debug)]
struct PolyNTT {
    coeffs: [u32; N],
}

struct Matrix {
    data: [[u32; N]; N],
}


impl Poly {
    fn new(coeffs: [u32; N]) -> Self {
        Self { coeffs }
    }

    fn add(&mut self, other: &Poly){
        for i in 0..N {
            self.coeffs[i] = (self.coeffs[i] + other.coeffs[i]) % Q;
        }
    }

    fn to_ntt(&self, matrix: &Matrix) -> PolyNTT {
        let NTT_matrix = Matrix::new(NTT);
        let result = NTT_matrix.multiply_poly(self);
        PolyNTT::new(result.coeffs)
    }
}

impl PolyNTT {
    fn new(coeffs: [u32; N]) -> Self {
        Self { coeffs }
    }
}

//[rows][columns]
impl Matrix {
    fn new(data: [[u32; N]; N]) -> Self {
        Self { data }
    }

    fn multiply_poly(&self, poly: &Poly) -> Poly {
        let mut result = [0u32; N];
        for i in 0..N {
            for j in 0..N {
                result[i] = result[i] + self.data[i][j] * poly.coeffs[j];
            }
            result[i] %= Q;
        }
        Poly::new(result)
    }
}

// just for checking
pub fn NTT_matrix() -> Matrix {
    let mut data = [[0u32; N]; N];
    for i in 0..N {
        for j in 0..N {
            data[i][j] = mod_pow(W, (i * j) as u32, Q);
        }
    }
    Matrix::new(data)
}

pub fn mod_pow(base: u32, exp: u32, modulus: u32) -> u32 {
    let mut result = 1;
    let mut base = base % modulus;
    let mut exp = exp;

    while exp > 0 {
        if exp % 2 == 1 {
            result = (result * base) % modulus;
        }
        base = (base * base) % modulus;
        exp >>= 1;
    }
    result
}



#[cfg(test)]
mod tests{
    use super::*;

    #[test]
    fn test_poly_add() {
        let mut Poly = Poly::new([1, 2, 3, 4]);
        let Poly2 = Poly::new([4, 3, 2, 1]);
        Poly.add(&Poly2);
        assert_eq!(Poly.coeffs, [5, 5, 5, 5]);
    }

    #[test]
    fn test_poly_add_mod() {
        let mut Poly = Poly::new([1 + 7680, 2 + 7680, 3 + 7680, 4 + 7680]);
        let Poly2 = Poly::new([4, 3, 2, 1]);
        Poly.add(&Poly2);
        assert_eq!(Poly.coeffs, [4, 4, 4, 4]);
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 5, 13), 9);
        assert_eq!(mod_pow(5, 0, 7), 1);
    }

    #[test]
    fn test_eval_NTT_matrix() {
        let matrix = NTT_matrix();
        assert_eq!(matrix.data, NTT);
    }

    #[test]
    fn test_matrix_multiply_poly() {
        let matrix = Matrix::new([[1, 1   , 1   , 1   ],
                                  [1, 3383, 7680, 4298],
                                  [1, 7680, 1   , 7680],
                                  [1, 4298, 7680, 3383]]);
        let poly = Poly::new([1, 2, 3, 4]);
        let result = matrix.multiply_poly(&poly);
        assert_eq!(result.coeffs, [10, 913, 7679, 6764]);
    }

    #[test]
    fn test_poly_to_ntt() {
        let poly = Poly::new([1, 2, 3, 4]);
        let matrix = Matrix::new(NTT);
        let poly_ntt = poly.to_ntt(&matrix);
        assert_eq!(poly_ntt.coeffs, [10, 913, 7679, 6764]);
    }

}