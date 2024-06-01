use num_bigint::BigUint;

//addition to num_bigint by https://github.com/d3crvpt3d
//using this is going to create a bigger number for each, you should use "reduce" at the end for memory/performance reasons

#[allow(unused)]

impl BigURational {
  
  //reducing the term by using GCD (Greatest common divisor)
  pub fn reduce(&mut self){

    let gcd: BigUint = gcd(&self.numerator, &self.denominator);

    self.numerator = &self.numerator / &gcd;
    self.denominator = &self.denominator / &gcd;
  }

  //creating a new instance of BigURational
  pub fn new(numerator: BigUint, denominator: BigUint) -> BigURational{
    BigURational{
      numerator,
      denominator,
    }
  }

  //returns as equal denominator for comparison (should be checked sometime)
  pub fn as_equal_denominator(&self, rhs: &Self) -> (BigURational, BigURational){

    let mut self_num: BigUint = self.numerator.to_owned();
    let mut rhs_num: BigUint = rhs.numerator.to_owned();
    let mut self_dom: BigUint = self.denominator.to_owned();
    let mut rhs_dom: BigUint = rhs.denominator.to_owned();

    if self_dom == rhs_dom {
      
      return (
        BigURational{
        numerator: self_num,
        denominator: self_dom,
      },
      BigURational{
        numerator: rhs_num,
        denominator: rhs_dom,
      });
    }else if self_dom > rhs_dom{
      self_num = &self_num * &rhs_dom;
      
    }else{
      rhs_num = &self_dom * &rhs_num;

    }

    let denom: BigUint = &self_dom * &rhs_dom;
    
    self_dom = denom.clone();
    rhs_dom = denom;

    (
      BigURational{
      numerator: self_num,
      denominator: self_dom,
    },
    BigURational{
      numerator: rhs_num,
      denominator: rhs_dom,
    })

    //TODO
  }

  //
  //comparison functions start
  //

  //comparisons (has to be mut, because of auto reduce for performance/memory reasons)
  //TODO: change mutability for comparisons
  //just checking if both values are equal, use eq_mut for greedy version
  //use as_equal_denominator and compare numerators 
  

  pub fn eq(&self, rhs: &Self) -> bool{
    let (a, b) = self.as_equal_denominator(rhs);

    a.numerator.eq(&b.numerator)
  }

  pub fn lt(&self, rhs: &Self) -> bool{
    let (a, b) = self.as_equal_denominator(rhs);

    a.numerator.lt(&b.numerator)
  }

  pub fn gt(&self, rhs: &Self) -> bool{
    let (a, b) = self.as_equal_denominator(rhs);

    a.numerator.gt(&b.numerator)
  }

  //comparison functions end
  //
  //mut comparison functions start

  pub fn eq_reduce(&mut self, rhs: &mut Self) -> bool{
    self.to_equal_denominator(rhs);
    self.reduce();

    self.numerator.eq(&rhs.numerator)
  }

  //less than reduce
  pub fn lt_reduce(&mut self, rhs: &mut Self) -> bool{
    self.to_equal_denominator(rhs);
    self.reduce();

    self.numerator.lt(&rhs.numerator)
  }

  //short for reduce then gt
  pub fn gt_reduce(&mut self, rhs: &mut Self) -> bool{
    self.to_equal_denominator(rhs);
    self.reduce();

    self.numerator.gt(&rhs.numerator)
  }

  //
  //mut comparison functions end
  //

  //changes numerator of BigURational with bigger denominator
  pub fn to_equal_denominator(&mut self, rhs: &mut Self){

    if self.denominator == rhs.denominator {
      
      return;
    }else if self.denominator > rhs.denominator{
      self.numerator = &self.numerator * &rhs.denominator;
      
    }else{
      rhs.numerator = &self.denominator * &rhs.numerator;

    }

    let denom: BigUint = &self.denominator * &rhs.denominator;
    
    self.denominator = denom.clone();
    rhs.denominator = denom;
  
  }

}


//greatestCommonDivisor for BigUInt using Euclidean algorythm
fn gcd(in_a: &BigUint, in_b: &BigUint) -> BigUint{

  let mut a: BigUint = in_a.clone();
  let mut b: BigUint = in_b.clone();

  while b.eq(&BigUint::ZERO) {
      let t: BigUint = b.clone();
      b = a % b;
      a = t;
  }

  a.to_owned()
}


//derived it myself
//look in Desmos: \frac{a_{0}}{b_{0}}+\frac{a_{1}}{b_{1}}=\frac{\left(a_{0}\cdot b_{1}\right)+\left(a_{1}\cdot b_{0}\right)}{b_{0}\cdot b_{1}}
impl std::ops::Add for BigURational {
    type Output = Self;
    
    fn add(self, rhs: Self) -> Self::Output {

      Self{
        numerator: (&self.numerator*&rhs.denominator)+(&rhs.numerator*&self.denominator),
        denominator: &self.denominator * &rhs.denominator,
      }

    }
}

//derived it myself too
//look in Desmos: \frac{a_{0}}{b_{0}}\cdot\frac{a_{1}}{b_{1}}=\frac{a_{0}\cdot a_{1}}{b_{0}\cdot b_{1}}
impl std::ops::Mul for BigURational {
    type Output = Self;

    fn mul(self, rhs: Self) -> Self::Output {
        
      Self{
        numerator: self.numerator * rhs.numerator,
        denominator: self.denominator * rhs.denominator,
      }

    }
}


//same as multiply, but rhs is swapped ( a/b = a*(1/b) )
impl std::ops::Div for BigURational {
    type Output = Self;

    fn div(self, rhs: Self) -> Self::Output {
        
        Self{
          numerator: self.numerator * rhs.denominator,
          denominator: self.denominator * rhs.numerator,
        }

    }
}

pub struct BigURational{
  numerator: BigUint,
  denominator: BigUint,
}