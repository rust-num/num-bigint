//! Test randomization of `BigUint` and `BigInt`
//!
//! This test is in a completely separate crate so `rand::rng()`
//! can be available without "infecting" the rest of the build with
//! `rand`'s default features, especially not `rand/std`.

#![cfg(test)]

mod torture;

mod biguint {
    use num_bigint::{BigUint, RandBigInt, RandomBits};
    use num_traits::Zero;
    use rand::distr::{Distribution, Uniform};
    use rand::{RngExt, SeedableRng};

    #[test]
    fn test_rand() {
        let mut rng = rand::rng();
        let n: BigUint = rng.random_biguint(137);
        assert!(n.bits() <= 137);
        assert!(rng.random_biguint(0).is_zero());
    }

    #[test]
    fn test_rand_bits() {
        let mut rng = rand::rng();
        let n: BigUint = rng.sample(&RandomBits::new(137));
        assert!(n.bits() <= 137);
        let z: BigUint = rng.sample(&RandomBits::new(0));
        assert!(z.is_zero());
    }

    #[test]
    fn test_rand_range() {
        let mut rng = rand::rng();

        for _ in 0..10 {
            assert_eq!(
                rng.random_biguint_range(&BigUint::from(236u32), &BigUint::from(237u32)),
                BigUint::from(236u32)
            );
        }

        let l = BigUint::from(403469000u32 + 2352);
        let u = BigUint::from(403469000u32 + 3513);
        for _ in 0..1000 {
            let n: BigUint = rng.random_biguint_below(&u);
            assert!(n < u);

            let n: BigUint = rng.random_biguint_range(&l, &u);
            assert!(n >= l);
            assert!(n < u);
        }
    }

    #[test]
    #[should_panic]
    fn test_zero_rand_range() {
        rand::rng().random_biguint_range(&BigUint::from(54u32), &BigUint::from(54u32));
    }

    #[test]
    #[should_panic]
    fn test_negative_rand_range() {
        let mut rng = rand::rng();
        let l = BigUint::from(2352u32);
        let u = BigUint::from(3513u32);
        // Switching u and l should fail:
        let _n: BigUint = rng.random_biguint_range(&u, &l);
    }

    #[test]
    fn test_rand_uniform() {
        let mut rng = rand::rng();

        let tiny = Uniform::new(BigUint::from(236u32), BigUint::from(237u32)).unwrap();
        for _ in 0..10 {
            assert_eq!(rng.sample(&tiny), BigUint::from(236u32));
        }

        let l = BigUint::from(403469000u32 + 2352);
        let u = BigUint::from(403469000u32 + 3513);
        let below = Uniform::new(BigUint::zero(), u.clone()).unwrap();
        let range = Uniform::new(l.clone(), u.clone()).unwrap();
        for _ in 0..1000 {
            let n: BigUint = rng.sample(&below);
            assert!(n < u);

            let n: BigUint = rng.sample(&range);
            assert!(n >= l);
            assert!(n < u);
        }
    }

    fn seeded_value_stability<R: SeedableRng + RandBigInt>(expected: &[&str]) {
        let mut seed = <R::Seed>::default();
        for (i, x) in seed.as_mut().iter_mut().enumerate() {
            *x = (i as u8).wrapping_mul(191);
        }
        let mut rng = R::from_seed(seed);
        for (i, &s) in expected.iter().enumerate() {
            let n: BigUint = s.parse().unwrap();
            let r = rng.random_biguint((1 << i) + i as u64);
            assert_eq!(n, r);
        }
    }

    #[test]
    fn test_chacha_value_stability() {
        const EXPECTED: &[&str] = &[
            "0",
            "5",
            "9",
            "1737",
            "937307",
            "65028028536",
            "279291567688057625466",
            "32684601755397150610585259058469172136086",
            "3318554627982869824783889169605932904164279211594869435303134228496248029647136",
            "2709999098147667596661525260287323466371098577491054736003573524257229599510813\
             950331137383286571179078064518241667831227807123801978731742492319680324618975",
        ];
        use rand_chacha::ChaChaRng;
        seeded_value_stability::<ChaChaRng>(EXPECTED);
    }

    #[test]
    fn test_isaac_value_stability() {
        const EXPECTED: &[&str] = &[
            "1",
            "7",
            "25",
            "1683",
            "117995",
            "67859584168",
            "518835239072611148581",
            "16943196501667364713622917423027819778745",
            "3128439869237800760464479228021661628664566245028432430864688183042066505636910",
            "1300089857736899601863677638636111069172421900222879296808895380997042087658914\
             2976856687389999376070549940495143742173349057217984875028392967823670091580",
        ];
        use rand_isaac::IsaacRng;
        seeded_value_stability::<IsaacRng>(EXPECTED);
    }

    #[test]
    fn test_xorshift_value_stability() {
        const EXPECTED: &[&str] = &[
            "1",
            "5",
            "4",
            "1412",
            "508011",
            "45408819789",
            "35449504201620804505",
            "1704901123711211464150080787386475139719",
            "3602478137269873610186753678058899426147621466964026295406154931693450025738094",
            "2676219517644961904040299206849473763923504769692259108134502736982105190584472\
             145963969095305013973645948123528298553311404505627077619646814679735036805715",
        ];
        use rand_xorshift::XorShiftRng;
        seeded_value_stability::<XorShiftRng>(EXPECTED);
    }

    #[test]
    fn test_roots_rand() {
        fn check<T: Into<BigUint>>(x: T, n: u32) {
            let x: BigUint = x.into();
            let root = x.nth_root(n);
            println!("check {}.nth_root({}) = {}", x, n, root);

            if n == 2 {
                assert_eq!(root, x.sqrt())
            } else if n == 3 {
                assert_eq!(root, x.cbrt())
            }

            let lo = root.pow(n);
            assert!(lo <= x);
            assert_eq!(lo.nth_root(n), root);
            if !lo.is_zero() {
                assert_eq!((&lo - 1u32).nth_root(n), &root - 1u32);
            }

            let hi = (&root + 1u32).pow(n);
            assert!(hi > x);
            assert_eq!(hi.nth_root(n), &root + 1u32);
            assert_eq!((&hi - 1u32).nth_root(n), root);
        }

        let mut rng = rand::rng();
        let bit_range = Uniform::new(0, 2048).unwrap();
        let sample_bits: Vec<_> = bit_range.sample_iter(&mut rng).take(100).collect();
        for bits in sample_bits {
            let x = rng.random_biguint(bits);
            for n in 2..11 {
                check(x.clone(), n);
            }
            check(x.clone(), 100);
        }
    }
}

mod bigint {
    use num_bigint::{BigInt, RandBigInt, RandomBits};
    use num_traits::Zero;
    use rand::distr::Uniform;
    use rand::{RngExt, SeedableRng};

    #[test]
    fn test_rand() {
        let mut rng = rand::rng();
        let n: BigInt = rng.random_bigint(137);
        assert!(n.bits() <= 137);
        assert!(rng.random_bigint(0).is_zero());
    }

    #[test]
    fn test_rand_bits() {
        let mut rng = rand::rng();
        let n: BigInt = rng.sample(&RandomBits::new(137));
        assert!(n.bits() <= 137);
        let z: BigInt = rng.sample(&RandomBits::new(0));
        assert!(z.is_zero());
    }

    #[test]
    fn test_rand_range() {
        let mut rng = rand::rng();

        for _ in 0..10 {
            assert_eq!(
                rng.random_bigint_range(&BigInt::from(236), &BigInt::from(237)),
                BigInt::from(236)
            );
        }

        fn check(l: BigInt, u: BigInt) {
            let mut rng = rand::rng();
            for _ in 0..1000 {
                let n: BigInt = rng.random_bigint_range(&l, &u);
                assert!(n >= l);
                assert!(n < u);
            }
        }
        let l: BigInt = BigInt::from(403469000 + 2352);
        let u: BigInt = BigInt::from(403469000 + 3513);
        check(l.clone(), u.clone());
        check(-l.clone(), u.clone());
        check(-u, -l);
    }

    #[test]
    #[should_panic]
    fn test_zero_rand_range() {
        rand::rng().random_bigint_range(&BigInt::from(54), &BigInt::from(54));
    }

    #[test]
    #[should_panic]
    fn test_negative_rand_range() {
        let mut rng = rand::rng();
        let l = BigInt::from(2352);
        let u = BigInt::from(3513);
        // Switching u and l should fail:
        let _n: BigInt = rng.random_bigint_range(&u, &l);
    }

    #[test]
    fn test_rand_uniform() {
        let mut rng = rand::rng();

        let tiny = Uniform::new(BigInt::from(236u32), BigInt::from(237u32)).unwrap();
        for _ in 0..10 {
            assert_eq!(rng.sample(&tiny), BigInt::from(236u32));
        }

        fn check(l: BigInt, u: BigInt) {
            let mut rng = rand::rng();
            let range = Uniform::new(l.clone(), u.clone()).unwrap();
            for _ in 0..1000 {
                let n: BigInt = rng.sample(&range);
                assert!(n >= l);
                assert!(n < u);
            }
        }
        let l: BigInt = BigInt::from(403469000 + 2352);
        let u: BigInt = BigInt::from(403469000 + 3513);
        check(l.clone(), u.clone());
        check(-l.clone(), u.clone());
        check(-u, -l);
    }

    fn seeded_value_stability<R: SeedableRng + RandBigInt>(expected: &[&str]) {
        let mut seed = <R::Seed>::default();
        for (i, x) in seed.as_mut().iter_mut().enumerate() {
            *x = (i as u8).wrapping_mul(191);
        }
        let mut rng = R::from_seed(seed);
        for (i, &s) in expected.iter().enumerate() {
            let n: BigInt = s.parse().unwrap();
            let r = rng.random_bigint((1 << i) + i as u64);
            assert_eq!(n, r);
        }
    }

    #[test]
    fn test_chacha_value_stability() {
        const EXPECTED: &[&str] = &[
            "0",
            "1",
            "27",
            "-2031",
            "194831",
            "-107270621334",
            "604778049616971649784",
            "-13279906003936626507900036341129035843658",
            "9598537977256701954366948527195592617018633007593014285361020521520881912967196",
            "-264136370194186516137263361513195811753545618309986722532270577136898598348583\
             0633637046727891746431626435333970024327119474183088397436173223513443625647336",
        ];
        use rand_chacha::ChaChaRng;
        seeded_value_stability::<ChaChaRng>(EXPECTED);
    }

    #[test]
    fn test_isaac_value_stability() {
        const EXPECTED: &[&str] = &[
            "-1",
            "1",
            "-43",
            "-1167",
            "593159",
            "87888368313",
            "863752686698588599066",
            "25660353483558192035558343375232549481936",
            "28977159681084521074061473674568372684301464954366590951863406150909603461154587",
            "87978171425457850961822054897752872385195765083062159849110972253115624397595266\
             8338609763273154793791639422990440190830083511986328340851006045387434229689",
        ];
        use rand_isaac::IsaacRng;
        seeded_value_stability::<IsaacRng>(EXPECTED);
    }

    #[test]
    fn test_xorshift_value_stability() {
        const EXPECTED: &[&str] = &[
            "1",
            "4",
            "-43",
            "810",
            "-389184",
            "127016555143",
            "860536200479125799740",
            "9835480991915491206245016261774321284293",
            "-14103609324431110452261709176070922521075146602747109938160234322089999252159519",
            "-22299954340459823432328915829626283467388738618069360141644037252844203947132116\
             47324956509902805306783167903408806315469788533784856707400727767026773441070",
        ];
        use rand_xorshift::XorShiftRng;
        seeded_value_stability::<XorShiftRng>(EXPECTED);
    }

    #[test]
    fn test_random_shr() {
        use rand::distr::StandardUniform;
        let rng = rand::rng();

        for p in rng.sample_iter::<i64, _>(&StandardUniform).take(1000) {
            let big = BigInt::from(p);
            let bigger = &big << 1000;
            assert_eq!(&bigger >> 1000, big);
            for i in 0..64 {
                let answer = BigInt::from(p >> i);
                assert_eq!(&big >> i, answer);
                assert_eq!(&bigger >> (1000 + i), answer);
            }
        }
    }
}
