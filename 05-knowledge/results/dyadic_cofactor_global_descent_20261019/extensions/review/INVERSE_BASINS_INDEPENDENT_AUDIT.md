# Independent audit of inverse-basin divisor formulas

Date: 19 October 2026.

## Verdict

**Pass. No repair required.** The general inverse formula, fixed-point
divisor-index characterization, power-of-two divisor-index characterization,
native-domain restrictions, numerical example, and two resonance connections
are all correct.

Audited source:

`INVERSE_BASIN_DIVISOR_LATTICES.md` before its documented status-only publication update

SHA-256:

`75692a6c84baf430c7a5c39ac571585e2b33dc91feedbabcef8ac42101021c8a`

The conclusions describe the auxiliary cofactor map, not arbitrary
aliquot inverse basins.

## General inverse compiler

Direct rearrangement gives

    F(cp)-1=(c-1)(p+2).

For a specified exponent `a>=1`, an even target `m` therefore has exactly one
possible odd part, `(m-1)/(2^a-1)-2`. It is a positive odd integer precisely
when the quotient is an integer at least 3. The quotient is odd because both
numerator and denominator are odd. Consequently `2^a` is the exact dyadic
factor of the constructed preimage. Distinct exponents cannot yield the same
positive integer, since its exact dyadic factor is unique.

## Fixed-point preimages

For `m_A=2^A(2^(A+1)-1)`, independent expansion confirms

    m_A-1=(2^A-1)(2^(A+1)+1).

If `a|A`, divisibility follows and the quotient is at least `2^(A+1)+1>=5`.
Thus all positive divisors give valid preimages; `a=A` gives the fixed point.

For necessity, let `r=A mod a`, and `M=2^a-1`. Reduction of the exponents
gives the claimed congruence

    m_A-1 ≡2^(2r+1)-2^r-1 (mod M).

If `r>0`, then `a>=2`. In the case `2r+1<a`, the displayed representative is
strictly positive and below `2^(a-1)<M`.

Otherwise put `s=2r+1-a`. We have `0<=s<=r<=a-1`, so after reducing once the
representative satisfies

    -2^r <=2^s-2^r-1 <=-1.

In particular it lies strictly between `-M` and 0, since
`2^r<=2^(a-1)<2^a-1`. Thus neither branch can vanish modulo `M`.
Divisibility forces `r=0`. The separate `a=1` case already has `r=0`.

This proves necessity even for candidate exponents `a>A`; none were omitted.
There are exactly `tau(A)` immediate preimages, with no unproved prime-divisor
theorem or factorization hypothesis needed.

## Prime-domain controls and example

The formula

    p_1=m_A-3=(2^A+1)(2^(A+1)-3)

is exact, and both factors exceed one when `A>=2`. For `A=1`, the formula
instead gives the fixed point 6. This excludes the core-2 alternative as a
new native predecessor, but does not exclude other divisor indices.

For `A=4`, the three inverse odd parts are indeed `493=17*29`, 163, and 31.
The latter two are prime. Independent divisor enumeration gives

    s(652)=496, s(496)=496, s(986)=634.

The auxiliary map sends all three inputs 652, 496, and 986 to 496. Thus the
example accurately distinguishes native inputs from a composite-odd-part
formal preimage. The manuscript correctly limits the divisor indexing to
immediate fixed-point preimages rather than claiming that the same indices
classify every later level of the inverse tree.

## Power-of-two preimages

Euclidean division `b=qa+r`, `0<=r<a`, gives

    2^b-1 ≡2^r-1 (mod 2^a-1).

For `r>0` this remainder lies strictly between zero and the modulus, proving
the standard divisibility equivalence. The case `a=b` gives quotient one
and odd part `-1`, so it is correctly excluded. For a proper divisor,
`b/a>=2`, and the geometric quotient is at least `1+2^a>=3`.
The exact count `tau(b)-1` follows, including the empty case `b=1`.

Writing `d=2^t`, the `b=2t,a=t` preimage is `d(d-1)`. For `b=3t,a=t` it
is `d(d²+d-1)`, which is exactly the first image of `d²(d-1)`. Direct
substitution verifies both resonance links and preserves exact cofactors.

## Independent exact controls

Files:

- `audit_inverse_basins.py`
- `inverse_basins_audit_results.json`
- `inverse_basins_audit.out`
- `inverse_basins_audit.optimized.out`

All checks use exact integers, including primality checks for the numerical
example. Explicit errors remain enabled under Python optimization. Normal
and optimized runs give identical PASS output.

- Complete inverse sets agree with independent forward enumeration for all
  50,000 positive even targets through 100,000. Completeness of the finite
  forward search follows from `F(n)>=n/2+3`, hence `n<=2m-6` whenever
  `F(n)=m`.
- The fixed-point preimage indices agree with the positive divisors of `A`
  for every `1<=A<=512`, including targets up to 1,025 bits.
- Power-target preimage indices agree with the proper positive divisors of
  `b` for every `1<=b<=1024`.
- 141,072 modular-residue cases verify both strict sign intervals and the
  exact divisibility equivalence, including seeded exponents through 4,096.
- Both resonance connections are checked through `t=512`.
- The 496 example is independently verified using proper-divisor sums and
  exact primality tests.

These controls are supplementary; the casewise modular proof provides the
universal result. This review does not establish novelty or priority.
