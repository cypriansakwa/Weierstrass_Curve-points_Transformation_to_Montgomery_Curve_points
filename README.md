# Weierstrass to Montgomery Curve Transformation

This Rust project implements a transformation from Weierstrass curves to Montgomery curves over a finite field. The code leverages the `num-bigint`, `num-integer`, `num-traits`, and `rand` crates for mathematical operations and random number generation.

## Overview

The main functionality includes:
- Computing the modular inverse using the extended Euclidean algorithm.
- Finding the modular square root using the Tonelli-Shanks algorithm.
- Transforming a point on a Weierstrass curve to its equivalent on a Montgomery curve.

## Dependencies

This project requires the following external crates:
- `num-bigint`
- `num-integer`
- `num-traits`
- `rand`

To add these dependencies, include the following in your `Cargo.toml`:

 >```
 >[dependencies]
 >num-bigint = "0.4"
 >num-integer = "0.1"
 >num-traits = "0.2"
 >rand = "0.8"
## Functions
- `mod_inverse(value: &BigInt, modulus: &BigInt) -> Option<BigInt>`
  - Computes the modular inverse of `value` modulo `modulus` using the extended Euclidean algorithm. Returns `None` if the inverse does not exist.
- `extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt)`
  - Implements the extended Euclidean algorithm. Returns a tuple containing the greatest common divisor (gcd), and the coefficients $x$ and $y$ such that:
    $\gcd=a\cdot x+b\cdot y$
- `mod_sqrt(value: &BigInt, p: &BigInt) -> Option<BigInt>`
  - Calculates the modular square root of `value` modulo `p` using the Tonelli-Shanks algorithm. Returns `None` if no square root exists.
- `transform_to_montgomery(x: &BigInt, y: &BigInt, a: &BigInt, b: &BigInt, p: &BigInt) -> Option<(BigInt, BigInt, BigInt, BigInt)>`
  - Transforms a point $(x,y)$ on a Weierstrass curve defined by parameters $a$ and $b$ over the field $\mathbb{F}\_p$ to its equivalent point on a Montgomery curve. Returns the transformed coordinates and curve parameters 
$a_{\text{montgomery}}$ and $b_{\text{montgomery}}$.
## Usage
To use the transformation, instantiate the parameters of your Weierstrass curve and the point you wish to transform. Below is an example of how to use the transformation function:
>```
>fn main() {
 >   // Example values for a Weierstrass curve over F_p
  >  let a = BigInt::from_str("8").unwrap();
   > let b = BigInt::from_str("2").unwrap();
    >let p = BigInt::from_str("17").unwrap(); // Example prime modulus

    >let x = BigInt::from_str("14").unwrap();
    >let y = BigInt::from_str("6").unwrap();
>
 >   match transform_to_montgomery(&x, &y, &a, &b, &p) {
  >      Some((x_montgomery, y_montgomery, a_montgomery, b_montgomery)) => {
   >         println!("x_montgomery: {}", x_montgomery);
    >        println!("y_montgomery: {}", y_montgomery);
     >       println!("a_montgomery: {}", a_montgomery);
      >      println!("b_montgomery: {}", b_montgomery);
       > }
        >None => println!("No valid transformation found."),
    >}
>}
