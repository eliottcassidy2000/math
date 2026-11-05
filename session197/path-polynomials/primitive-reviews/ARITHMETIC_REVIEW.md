# Independent review: primitive target inverse and variable leaf

**PASS, with no correction.** This review covers the two addenda's arithmetic and all-order arguments. Their physical family identity is the explicitly named dependency `PRIMITIVE_CANCELLATION.md`; its independent graph/profile audit belongs to the separate primitive-family review. The conclusions here do not assert anything about tournaments outside that language.

## Frozen sources

In `session197-odd-cycle-path-polynomials/`:

- `PRIMITIVE_TARGET_INVERSE.md`: `7d2378f8dbc4d1fcf76af6e30f84485a1796a7d9f80210a5c9d6825dd034d614`.
- `verify_primitive_target_inverse.py`: `caaff817de21796d5486f0b9c574353dcf889b72ecd402bea7a5019178ad40eb`.
- `primitive_target_normal.stdout`: `82484ee2f8d4edc967f874c064a65d0e63c06258f8c00b5fb4d25281e5368f27`.
- `PRIMITIVE_LEAF_RIGIDITY.md`: `cac16f4a0a4fd82370c2abb18817e3e3c4f1881ec7e9aff7820a3477f90bed8b`.
- `verify_primitive_leaf_rigidity.py`: `f5fb267b7c2c056de3910407ab86caad2381370bb3676616524c79b0778c2ebc`.
- `primitive_leaf_normal.stdout`: `0a512c956f9ac7f7a665a419858dde1fd15e38210dbb284d15b09919ad902f4b`.

The dependency proof and verifier retain hashes `6c577db04d939196b04cbec9ba322f7f5abff7ae95e2f0fa19f7365c3acfd6ea` and `535d68803ccce6f1253e7ebf27c797bd9097665c2edc6e72cf897994c72ac0ce`.

## Target inverse

Write W=4^t bx^(2s) by^(2j) 3^(m-1) and V=hx^(2s) hy^(2j). The domains t,s,j>=0 and m>=1, and bx>hx, by>hy imply W>=V. Thus N=3A0 W-D0 V >= (3A0-D0)W. The positive factor H0 makes the finite s,j cutoff valid, including equality at the minimum recipe. After those exponents are fixed, the positive integer quotient (N+d)/a has a unique representation 4^t3^m exactly under the stated valuation test. No zero quotient or nontermination boundary can occur. The retained X in the left module and Y in the right make every recipe nonempty. Distinct tuples need not be presumed to have distinct counts.

For H0+6 the cutoff forces s=j=0 and then 4^t3^(m-1)<2 forces t=0,m=1. This contradicts the target value, independently of the implementation. The absence is all-order inside the language. The primality statement is separately checked.

## Variable leaf

The ratio group has generators with valuation vectors (at 2,11,31,37,127) equal to (1,0,0,0,0), (4,1,-1,-1,1), (5,1,0,-1,0). The two displayed linear relations and prime support are therefore necessary and sufficient, with the explicit integer reconstruction in the proof. Reducing the Mersenne numerator and denominator removes exactly gcd 1 or 3; their odd reduced parts remain coprime.

The order-36 claim is exact: every proper divisor of 36 divides 18 or 12, and neither checked power equals 1. A factor 37 in either reduced odd part forces its exponent divisible by 36. Its unreduced numerator is then divisible by 9, while reduction removes only a single factor 3. The other odd part is coprime, so this surviving factor 3 cannot cancel. It is outside the bank support. Consequently both 37 and, by the valuation relation, 11 disappear. The remaining odd ratio is (127/31)^a. Its strict interval (4,7] forces a=1 and then r=7. Cases r=1,2 have the correctly excluded cancellation signs. The converse uses actual nonempty X and Y. All leaf lengths are covered without factoring arbitrary large Mersenne numbers.

## Independent controls and replays

`check.py` independently enumerates the complete count table up to 10^32 by stopping on **actual increasing counts**, rather than using the source H0*W bound. It finds 2,388 recipes and values and compares the supplied inverse on 13,940 members, neighbours and deterministic additional targets. It independently sieves and tests the 7,604 primes through 77,395 to certify the omitted prime 5,990,037,979. A separate rational-factor implementation checks all 8,190 leaf lengths 3 through 8,192 and reconstructs the unique accepted exponent recipe (r,zK,zX,zY)=(7,0,1,-1).

- Independent checker: `fab6b9ed2c11e622867d691925ada4e721676f7fdda8242aaad2c562f283c66e`.
- Independent normal/optimized receipt: `cc9b6d425ca549e7d7196321fd4fa4b34ab0a1a7acb6ee669ec6b9168ab658f5`.

Both source verifiers were copied with their named dependencies into this audit directory before replay. Normal and optimized output is byte-identical to the frozen source receipts. No source was edited. Finite controls support the implementation; the arguments above, rather than the finite range, justify the unbounded claims.
