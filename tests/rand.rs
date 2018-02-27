#![cfg(feature = "rand")]

extern crate num_bigint;
extern crate num_traits;
extern crate rand;

mod biguint {
    use rand::thread_rng;
    use num_bigint::{BigUint, RandBigInt};
    use num_traits::{FromPrimitive, Zero};

    #[test]
    fn test_rand() {
        let mut rng = thread_rng();
        let _n: BigUint = rng.gen_biguint(137);
        assert!(rng.gen_biguint(0).is_zero());
    }

    #[test]
    fn test_rand_range() {
        let mut rng = thread_rng();

        for _ in 0..10 {
            assert_eq!(rng.gen_bigint_range(&FromPrimitive::from_usize(236).unwrap(),
                                            &FromPrimitive::from_usize(237).unwrap()),
                       FromPrimitive::from_usize(236).unwrap());
        }

        let l = FromPrimitive::from_usize(403469000 + 2352).unwrap();
        let u = FromPrimitive::from_usize(403469000 + 3513).unwrap();
        for _ in 0..1000 {
            let n: BigUint = rng.gen_biguint_below(&u);
            assert!(n < u);

            let n: BigUint = rng.gen_biguint_range(&l, &u);
            assert!(n >= l);
            assert!(n < u);
        }
    }

    #[test]
    #[should_panic]
    fn test_zero_rand_range() {
        thread_rng().gen_biguint_range(&FromPrimitive::from_usize(54).unwrap(),
                                       &FromPrimitive::from_usize(54).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_negative_rand_range() {
        let mut rng = thread_rng();
        let l = FromPrimitive::from_usize(2352).unwrap();
        let u = FromPrimitive::from_usize(3513).unwrap();
        // Switching u and l should fail:
        let _n: BigUint = rng.gen_biguint_range(&u, &l);
    }
}

mod bigint {
    use rand::thread_rng;
    use num_bigint::{BigInt, RandBigInt};
    use num_traits::{FromPrimitive, Zero};

    #[test]
    fn test_rand() {
        let mut rng = thread_rng();
        let _n: BigInt = rng.gen_bigint(137);
        assert!(rng.gen_bigint(0).is_zero());
    }

    #[test]
    fn test_rand_range() {
        let mut rng = thread_rng();

        for _ in 0..10 {
            assert_eq!(rng.gen_bigint_range(&FromPrimitive::from_usize(236).unwrap(),
            &FromPrimitive::from_usize(237).unwrap()),
            FromPrimitive::from_usize(236).unwrap());
        }

        fn check(l: BigInt, u: BigInt) {
            let mut rng = thread_rng();
            for _ in 0..1000 {
                let n: BigInt = rng.gen_bigint_range(&l, &u);
                assert!(n >= l);
                assert!(n < u);
            }
        }
        let l: BigInt = FromPrimitive::from_usize(403469000 + 2352).unwrap();
        let u: BigInt = FromPrimitive::from_usize(403469000 + 3513).unwrap();
        check(l.clone(), u.clone());
        check(-l.clone(), u.clone());
        check(-u.clone(), -l.clone());
    }

    #[test]
    #[should_panic]
    fn test_zero_rand_range() {
        thread_rng().gen_bigint_range(&FromPrimitive::from_isize(54).unwrap(),
        &FromPrimitive::from_isize(54).unwrap());
    }

    #[test]
    #[should_panic]
    fn test_negative_rand_range() {
        let mut rng = thread_rng();
        let l = FromPrimitive::from_usize(2352).unwrap();
        let u = FromPrimitive::from_usize(3513).unwrap();
        // Switching u and l should fail:
        let _n: BigInt = rng.gen_bigint_range(&u, &l);
    }
}
