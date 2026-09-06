# Independent audit of odd four-jet precision and the Deuring residue packet

**PASS — complete written proof audit and separate exact controls.**
The four claims in
[the producer report](overnight8_20260906_jets_residue.md)
are accepted on their distinct stated domains. No repair was required.
This audit does not promote a general full Smith-partition claim.

The audited producer source has SHA256
`102225c5476d59b3aec02c5b479093d105298000b9411fdb4526f9423efd3b65`;
its frozen LF output has SHA256
`f77f84da627c5c1f33468326baf24cc93c1bcf222028a935a956728f1d878061`.
The producer's later audit-status prose may change without changing this
mathematical version. No producer files were edited by the referee.

## Accepted statements and their boundaries

1. For three integer nodes, four complete Hasse jets at each node, and
   every odd prime p, the sharp largest Smith exponent is
   `11e+7d-[p=5]` in the close case `d>=1`. In the equilateral case it
   is zero at e=0 and `11e-beta` at e>=1, where beta is one at p3,
   p5, or p7 with normalized residues in arithmetic progression, and
   zero otherwise.
2. At p7, only the equilateral four-jet family has the asserted complete
   list: `(0^4,e,2e,4e,5e+1,7e,8e,10e-1+kappa,11e-kappa)` for e>=1,
   with the separate all-zero list at e=0. The AP bit kappa is intrinsic.
3. For n>=2 nodes with every pair at depth e>=1, distinct normalized
   residues in F_p, and uniform multiplicity m=(p+1)/2, the coefficient
   packet is nonzero if and only if the largest exponent attains
   `(nm-1)e`. Packet zero supplies only a strict drop from that ceiling.
4. At n=3 the packet is, up to a global nonzero sign, the classical
   Deuring polynomial. The cited supersingular interpretation is correct;
   it is not needed in the elementary proof of the precision criterion.

Neither the general packet nor its scalar three-node specialization
determines the amount of a further drop or the intermediate Smith ideals.
The report does not assert a close-pair full p7 partition, an all-odd
four-jet full partition, or a fixed-multiplicity metric failure at every
odd prime. Its uniform prime family lets multiplicity vary with p.

The minimality statement is also correctly scoped: four is the least
uniform jet count permitting a three-node odd-prime metric counterexample.
This uses the inherited values-only and two-jet laws and the independently
audited seventh-round three-jet law. No integer-diameter or universal
smallest-prime claim is inferred.

## Proof-level checks

The inverse-reciprocal equality is inherited from
[THM-4443](../../01-canon/theorems/THM-4443-arbitrary-jet-precision-and-dyadic-unit-boundary.md),
already independently audited in the sixth-round normalized-jet work.
It concerns the actual largest coefficient-recovery exponent, including
attainment, rather than only a local cancellation statistic. Translation,
p-adic unit changes of variable, and node permutation are unimodular for
the complete Hasse observer. Rational normalized units therefore present
no missing integrality restriction.

For a closest-pair endpoint, let its two differences have valuations
e and e+d. In the cubic numerator
`20(u+v)(u^2+uv+v^2)`, the smallest term after removing 20 has valuation
3e and is unique. Its top-order inverse cost is exactly
`11e+7d-[p=5]`. The lower-order upper bounds are
`8e+4d`, `9e+5d`, `10e+6d`; none exceeds that top cost because
e+d>=1. The third node contributes at most 11e. This includes the
shallow quinary tie `p5,e0,d1` and has no unexamined residue exception.

In the equilateral case the three C cubics are correct up to irrelevant
signs and unit denominators. Their two sum/remainder identities force a
common odd-prime zero to occur only at p3 or p7. At p3 the two relevant
linear factors add to three, while their quadratic factors are units;
the simultaneous valuation is exactly one at every lift. At p7 their
common roots are 2,4,6, and the sum is seven times a unit there; again
the common valuation is exactly one. At other odd primes there is a
unit cubic, and the external factor20 adds exactly one at p5. Lower
orders cost at most10e, so they cannot defeat `11e-beta` for e>=1.
At e0, the unit determinant correctly overrides the positive-depth form.

For the complete p7 partition, elimination of the four rows at zero
leaves the specified 8x8 residual matrix and four unit factors. Every
minor has a single factored scale `7^(We)` times an integral polynomial
in the normalized unit a. The weight minima are exactly
`(1,3,7,12,19,27)`, and the report's nine polynomials are the entire
equality-weight set. At ranks4–6 those nine witnesses have the asserted
seven-divisibility; every other minor pays at least one extra weight.
The step `e>=1` is essential and is stated. The companion witnesses
attain the lower bounds jointly for every a,a-1 unit: the two rank3
linear factors differ by the unit3 after the specified combination,
and both rank5 quadratics are units modulo7. The determinant48e and
the separately proved largest exponent determine D11. Successive
differences give the complete list, sorted for every e>=1. Its last
two entries tie exactly at e1,kappa1.

For the general packet, the Frobenius identity on
`Q_i(t_i+T)^(-m)` is valid modulo `T^(k+1)` because k<p and
m+k=p. It identifies the normalized top coefficient with a unit times
`[T^(p-1)]F(t_i+T)^k`. The binomial coefficient in this translation
vanishes modulo p unless its exponent is `sp-1`, when it equals one.
The declared universe satisfies `nk<=p(p-1)/2<p^2`; indeed the elementary
binomial argument supplied in the report works even without that upper
bound. After Frobenius evaluation, the resulting polynomial in t_i has
degree q-1<n, so vanishing at all n distinct residues is equivalent to
zero packet. At e>=1 only the top order can attain the ceiling. Thus
both directions of the claimed iff hold. Arbitrary subsets of F_p are
allowed; they need not occupy every residue. The complete-residue
polynomial `X^p-X` is a valid simultaneous-cancellation hostile.

For three nodes, `deg F^k=3(p-1)/2<2p-1`, so the packet has one
coefficient and its translation invariance is valid at the needed
degree. Direct expansion gives `(-1)^k H_p(a)`. The reciprocal and
reflection identities preserve its zero set under all six parameter
transforms with the stated nonzero unit factors. The values a=0,1
are excluded throughout. When p>=7 and p=3 modulo4, H_p(-1)=0,
while its degree is less than p-2; therefore both precision behaviors
occur among admissible residue classes. This is an unbounded analytic
family, without any assertion about higher p-adic cancellation depth.

## Primary-source check

The referee independently opened and read Roland Auer and Jaap Top,
[*Legendre elliptic curves over finite fields*, arXiv:math/0106273v1](https://arxiv.org/pdf/math/0106273),
Section3, printed pages5–6, particularly the paragraph immediately
before Proposition3.2. It identifies the roots of their sign-normalized
Deuring polynomial with supersingular Legendre parameters. The producer
uses the same root set. Only this three-node interpretation is cited;
the paper supplies none of the new Smith formulas or a vector-packet
interpretation for arbitrary n. No external priority is asserted.

## Independent exact controls

The [referee engine](../../04-computation/overnight8_20260906_jets_independent_audit.py)
uses only the Python standard library and imports no producer. It has
three separate computational paths:

- Every one of the12,804 rank1–6 weights is enumerated; all nine
  minimal-weight determinant polynomials are reconstructed by integer
  Bareiss evaluation and finite-difference interpolation with proved
  degree bounds. An extra value at a=-3 checks each reconstruction.
  This is separate from the producer's symbolic determinant path.
- Ninety-one full four-jet Hasse matrices are reduced by least-valuation
  p-local Gaussian pivots. Every row multiplier is checked p-integral,
  so the pivots give Smith exponents over the local ring. The bank
  contains close configurations, complete equilateral residue fibers
  at p3,5,7,11,13, fresh higher lifts, and e0 controls. Every matrix
  also satisfies the exact confluent determinant valuation. The p7
  hostile pair and its kernel exponents48 versus47 are recovered.
- One hundred fifty-five coefficient-packet controls compare ordinary
  integer polynomial translation against reciprocal coefficients solved
  recursively from the defining polynomial inverse, without the
  producer's negative-binomial formula. Thirty-three full Hasse matrices
  independently check the actual ceiling iff. The Deuring identities,
  full residue root sets, and both S3 generators are checked at eight
  odd primes through31.

These banks challenge the analytic proof; they are not a full matrix
height census or an extrapolation to unsampled depths. The normal and
optimized referee outputs are byte-identical LF and pass **26,520
explicit gates**.

```text
python -B 04-computation/overnight8_20260906_jets_independent_audit.py
python -B -O 04-computation/overnight8_20260906_jets_independent_audit.py

source 4e3d99d9ef917ca5026a4348e71058471e54ee20ff8af91a19a5ffb026197ca1
output a32847a7def96f074895538bbc4cffc6af908eb975481871c9a41955f9e60efa
literal bank 412964b708a5a3a37967900c5290dd6d9ca451ea2804ac3140fd136e59fa8d6a
```

All referee files remain outside the repository for parent-managed
integration. The shared worktree and all incoming proof surfaces were
left untouched by this audit.

**Filing:** root integrated these audited artifacts in the eighth checkpoint;
reproduction paths are relative to the repository root. Earlier outside-worktree
notes describe author provenance, not the present file location.
