// based on: https://eprint.iacr.org/2024/585.pdf
use std::marker::PhantomData;

const N: usize = 4;   // in real 256
const IN: u32 = 5761;
const Q: u32 = 7681;    // in real 3329

// first root
const W: u32 = 3383;
const IW: u32 = 4298; // inverse of W mod Q
const NTT: [[u32; N]; N] = [[1, 1   , 1   , 1   ],
                            [1, 3383, 7680, 4298],
                            [1, 7680, 1   , 7680],
                            [1, 4298, 7680, 3383]]; // well its bigger in reality 256 long.

const INTT: [[u32; N]; N] = [[1, 1   , 1   , 1   ],
                             [1, 4298, 7680, 3383],
                             [1, 7680, 1   , 7680],
                             [1, 3383, 7680, 4298]]; // well its bigger in reality 256 long.

//second root
const W2: u32 = 1925;

const NWNTT: [[u32; N]; N] = [[1, 1925, 3383, 6468],
                              [1, 6468, 4298, 1925],
                              [1, 5756, 3383, 1213],
                              [1, 1213, 4298, 5756]]; // well its bigger in reality 256 long.


trait Polynomial {
    fn new(coeffs: [u32; N]) -> Self;
    fn get_coeffs(&self) -> &[u32; N];
    fn add(&mut self, other: &Self);
    fn scalar_mul(&mut self, scalar: u32);
    fn el_wise_mul(&self, other: &Self) -> Self;
}

#[derive(Clone, Debug)]
struct Poly_norm;
#[derive(Clone, Debug)]
struct Poly_ntt;

#[derive(Clone, Debug)]
struct Poly<Domain> {
    coeffs: [u32; N],
    _domain: PhantomData<Domain>,
}

type Poly_normal = Poly<Poly_norm>;

type Poly_NTT = Poly<Poly_ntt>;

struct Matrix {
    data: [[u32; N]; N],
}

impl <T:Clone> Polynomial for Poly<T> {
    fn new(coeffs: [u32; N]) -> Self {
        Self { coeffs, _domain: PhantomData }
    }

    fn get_coeffs(&self) -> &[u32; N] {
        &self.coeffs
    }

    fn add(&mut self, other: &Self) {
        for i in 0..N {
            self.coeffs[i] = (self.coeffs[i] + other.coeffs[i]) % Q;
        }
    }

    fn scalar_mul(&mut self, scalar: u32) {
        for i in 0..N {
            self.coeffs[i] = (self.coeffs[i] * scalar) % Q;
        }
    }

    fn el_wise_mul(&self, other: &Self) -> Self {
        let mut result = self.clone();
        for i in 0..N {
            result.coeffs[i] = (result.coeffs[i] * other.coeffs[i]) % Q;
        }
        result
    }
}


impl Poly_normal {
        fn convolution(&self, other: &Poly_normal) -> Poly_normal {
        let a = self.to_ntt();
        let b = other.to_ntt();
        let c = a.el_wise_mul(&b);
        c.to_normal()
    }

    fn to_ntt(&self) -> Poly_NTT {
        Poly_NTT::new(Matrix::new(NTT).multiply_poly(self.get_coeffs()))
    }
}

impl Poly_NTT {
    fn to_normal(&self) -> Poly_normal {
        let mut poly = Poly_normal::new(Matrix::new(INTT).multiply_poly(self.get_coeffs()));
        poly.scalar_mul(IN);
        poly
    }
}


//[rows][columns]
impl Matrix {
    fn new(data: [[u32; N]; N]) -> Self {
        Self { data }
    }

    fn multiply_poly(&self, poly: &[u32; N]) -> [u32; N] {
        let mut result = [0u32; N];
        for i in 0..N {
            for j in 0..N {
                result[i] = result[i] + self.data[i][j] * poly[j];
            }
            result[i] %= Q;
        }
        result
    }
}

// just for checking
pub fn ntt_matrix() -> Matrix {
    let mut data = [[0u32; N]; N];
    for i in 0..N {
        for j in 0..N {
            data[i][j] = mod_pow(W, (i * j) as u32, Q);
        }
    }
    Matrix::new(data)
}

// just for checking
pub fn intt_matrix() -> Matrix {
    let mut data = [[0u32; N]; N];
    for i in 0..N {
        for j in 0..N {
            data[i][j] = mod_pow(IW, (i * j) as u32, Q);
        }
    }
    Matrix::new(data)  
}

// just for checking
// in docs the exp is (2*i*j + i), but example is (2*i*j + j)
pub fn nwtt_matrix() -> Matrix {
    let mut data = [[0u32; N]; N];
    for i in 0..N {
        for j in 0..N {
            data[i][j] = mod_pow(W2, ((2*i * j) + j) as u32, Q);
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
        let mut poly = Poly_normal::new([1, 2, 3, 4]);
        let poly2 = Poly_normal::new([4, 3, 2, 1]);
        poly.add(&poly2);
        assert_eq!(poly.get_coeffs(), &[5, 5, 5, 5]);
    }

    #[test]
    fn test_poly_add_mod() {
        let mut poly = Poly_normal::new([1 + 7680, 2 + 7680, 3 + 7680, 4 + 7680]);
        let poly2 = Poly_normal::new([4, 3, 2, 1]);
        poly.add(&poly2);
        assert_eq!(poly.get_coeffs(), &[4, 4, 4, 4]);
    }

    #[test]
    fn test_mod_pow() {
        assert_eq!(mod_pow(2, 10, 1000), 24);
        assert_eq!(mod_pow(3, 5, 13), 9);
        assert_eq!(mod_pow(5, 0, 7), 1);
    }

    #[test]
    fn test_eval_ntt_matrix() {
        let matrix = ntt_matrix();
        assert_eq!(matrix.data, NTT);
    }

    #[test]
    fn test_eval_intt_matrix(){
        let matrix  = intt_matrix();
        assert_eq!(matrix.data, INTT);
    }

    #[test]
    fn test_matrix_multiply_poly() {
        let matrix = Matrix::new([[1, 1   , 1   , 1   ],
                                  [1, 3383, 7680, 4298],
                                  [1, 7680, 1   , 7680],
                                  [1, 4298, 7680, 3383]]);
        let poly = [1, 2, 3, 4];
        let result = matrix.multiply_poly(&poly);
        assert_eq!(result, [10, 913, 7679, 6764]);
    }

    #[test]
    fn test_poly_to_ntt() {
        let poly = Poly_normal::new([1, 2, 3, 4]);
        let poly_ntt = poly.to_ntt();
        assert_eq!(poly_ntt.get_coeffs(), &[10, 913, 7679, 6764]);
    }

    #[test]
    fn test_poly_ntt_to_normal() {
        let poly = Poly_NTT::new([10, 913, 7679, 6764]);
        let poly_normal = poly.to_normal();
        assert_eq!(poly_normal.get_coeffs(), &[1, 2, 3, 4]);
    }

    #[test]
    fn test_poly_transforma_reduction() {
        let poly = Poly_normal::new([1, 2, 3, 4]);
        let poly_ntt = poly.to_ntt();
        let poly_normal = poly_ntt.to_normal();
        assert_eq!(poly_normal.get_coeffs(), &[1, 2, 3, 4]);
    }

    // evaluating convolution of two polynomials
    #[test]
    fn test_convolution_with_ntt(){
        let g = Poly_normal::new([1, 2, 3, 4]);
        let h = Poly_normal::new([5, 6, 7, 8]);
        let conv = g.convolution(&h);
        assert_eq!(conv.get_coeffs(), &[66, 68, 66, 60]);
    }


    #[test]
    fn test_generation_nwntt_matrix(){
        let matrix = nwtt_matrix();
        assert_eq!(matrix.data, NWNTT);
    }
}