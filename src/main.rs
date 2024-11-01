extern crate num_bigint;
extern crate num_integer;
extern crate num_traits;
extern crate rand;

use num_bigint::{BigInt, RandBigInt};
use num_integer::Integer;
use num_traits::{One, Zero};
use std::str::FromStr;

/// Computes the modular inverse of `value` modulo `modulus` using the extended Euclidean algorithm.
fn mod_inverse(value: &BigInt, modulus: &BigInt) -> Option<BigInt> {
    let (gcd, x, _) = extended_gcd(value, modulus);
    if gcd != BigInt::one() {
        None
    } else {
        Some((x % modulus + modulus) % modulus)
    }
}

/// Computes the extended Euclidean algorithm, returning (gcd, x, y) such that gcd = value * x + modulus * y.
fn extended_gcd(a: &BigInt, b: &BigInt) -> (BigInt, BigInt, BigInt) {
    if b.is_zero() {
        (a.clone(), BigInt::one(), BigInt::zero())
    } else {
        let (gcd, x1, y1) = extended_gcd(b, &(a % b));
        (gcd, y1.clone(), x1 - (a / b) * y1)
    }
}

/// Computes the modular square root using the Tonelli-Shanks algorithm.
/// Returns `None` if no square root exists.
fn mod_sqrt(value: &BigInt, p: &BigInt) -> Option<BigInt> {
    if value.is_zero() {
        return Some(BigInt::zero());
    }
    if p == &BigInt::from(2) {
        return Some(value.clone());
    }
    if value.modpow(&((p - 1u32) / 2u32), p) != BigInt::one() {
        return None; // No square root exists
    }
    
    let mut q = p - 1u32;
    let mut s = 0;
    while q.is_even() {
        q /= 2u32;
        s += 1;
    }
    
    let z = (2u32..)
        .map(BigInt::from)
        .find(|n| n.modpow(&((p - 1u32) / 2u32), p) == p - 1u32)
        .unwrap();
    
    let mut m = s;
    let mut c = z.modpow(&q, p);
    let mut t = value.modpow(&q, p);
    let mut r = value.modpow(&((q + 1u32) / 2u32), p);

    while t != BigInt::one() {
        let mut i = 0;
        let mut t2i = t.clone();
        while t2i != BigInt::one() {
            t2i = t2i.modpow(&BigInt::from(2), p);
            i += 1;
            if i == m {
                return None;
            }
        }
        
        let b = c.modpow(&BigInt::from(1u32 << (m - i - 1)), p);
        m = i;
        c = b.modpow(&BigInt::from(2), p);
        t = (t * &c) % p;
        r = (r * b) % p;
    }
    
    Some(r)
}

/// Transformation function from Weierstrass to Montgomery curve.
fn transform_to_montgomery(
    x: &BigInt,
    y: &BigInt,
    a: &BigInt,
    b: &BigInt,
    p: &BigInt,
) -> Option<(BigInt, BigInt, BigInt, BigInt)> {
    // Find a root z0 of the polynomial z^3 + az + b in the field F_p
    let mut rng = rand::thread_rng();
    let z0 = loop {
        let candidate = rng.gen_bigint_range(&BigInt::zero(), p);
        if (&candidate.pow(3) + a * &candidate + b).mod_floor(p).is_zero() {
            //println!("Found z0: {}", candidate); // Debug print
            break candidate;
        }
    };

    // Compute s = (sqrt(3 * z0^2 + a))^{-1} modulo p
    let s_squared = (BigInt::from(3) * &z0 * &z0 + a).mod_floor(p);
    //println!("s_squared: {}", s_squared); // Debug print

    let s = mod_sqrt(&s_squared, p)?;
    //println!("s: {}", s); // Debug print

    let s_inv = mod_inverse(&s, p)?;
    //println!("s_inv: {}", s_inv); // Debug print

    // Compute the new parameters a and b
    let a_montgomery = (BigInt::from(3) * &z0 * &s_inv).mod_floor(p);
    let b_montgomery = s_inv.mod_floor(p);

    // Map (x, y) to (x_montgomery, y_montgomery) on the Montgomery curve
    let x_montgomery = (s_inv.clone() * (x - &z0)).mod_floor(p); // Clone s_inv for the first usage
    let y_montgomery = (s_inv * y).mod_floor(p);                 // Use s_inv directly for the second usage
    Some((x_montgomery, y_montgomery, a_montgomery, b_montgomery))
}

fn main() {
    // Example values for a Weierstrass curve over F_p
    let a = BigInt::from_str("8").unwrap();
    let b = BigInt::from_str("2").unwrap();
    let p = BigInt::from_str("17").unwrap(); // Example prime modulus

    let x = BigInt::from_str("14").unwrap();
    let y = BigInt::from_str("6").unwrap();

    match transform_to_montgomery(&x, &y, &a, &b, &p) {
        Some((x_montgomery, y_montgomery, a_montgomery, b_montgomery)) => {
            println!("x_montgomery: {}", x_montgomery);
            println!("y_montgomery: {}", y_montgomery);
            println!("a_montgomery: {}", a_montgomery);
           println!("b_montgomery: {}", b_montgomery);
        }
        None => println!("No valid transformation found."),
    }
}
