---
id: THM-4440
title: "Signed duplication SOS and real-rooted Laurent return"
status: PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED
source: overnight-hexagon-sep05 fifth research wave
depends_on:
  - THM-4436-complete-factorial-row-simple-negative-roots-and-trinomial-phase-wall
proof: 05-knowledge/results/nc2_signed_duplication_overnight_hexagon_sep05.md
script: 04-computation/nc2_signed_duplication_overnight_hexagon_sep05.py
output: 05-knowledge/results/nc2_signed_duplication_overnight_hexagon_sep05.out
script_sha256: 654cabc326769d3852f35524ce143444f4c03aad7edf52fffa5b9bb4049a6c96
output_sha256: 60ee438bc8dd80355b4ad16a17e210d7a64f634b32951015b3bf4be0135e519a
hash_basis: raw LF bytes
---

# THM-4440 -- Signed duplication SOS and real-rooted Laurent return

**PROVED ANALYTICALLY + FINITE-EXACT + INDEPENDENTLY AUDITED.** The
[complete proof, consumers, hostile controls and classical Hermite comparison](../../05-knowledge/results/nc2_signed_duplication_overnight_hexagon_sep05.md)
are part of this theorem. The SOS has an elementary self-contained proof;
THM-4436 is needed for the full complex AP-phase consumer. No external
priority is claimed, and the Hermite identity is explicitly recovered from
classical Mehler, with an additional finite-subset limit derivation.

For arbitrary real r_1,...,r_n, signs and zeros allowed, and1<=k<=n, put
D=binom(n,k), C=binom(2n,2k)/D^2, beta_s=binom(n-s,k-s)/D and

```text
c_(k,s)=2^(2k-2s+1)(2s-2)!(k-s)!k!/[(s-1)!(2k)!].
```

Writing r^dup for the list with every entry repeated twice gives the exact
positive rational sum of squares

```text
C e_k(r)^2-e_(2k)(r^dup)
 =sum_(s=1)^k c_(k,s) sum_(|S|=s)
   [r_S e_(k-s)(r outside S)-beta_s e_k(r)]^2.
```

Equality means all k-subset products agree. For k<n this is exactly that
fewer than k entries are nonzero or all entries are equal; for k=n equality
always holds. The final square level has coefficient1/(2k-1).

The mechanism is the exact subset-distance matrix K_IJ=4^ell/binom(2ell,ell),
ell=|I\J|. Its quadratic form is the duplicated elementary coefficient,
not an arbitrary norm. Centering its constant mode and expanding
q^|I\J| over all common subsets gives the positive integral and every
square level. Individual opposite-coefficient products can be positive.

For a nonzero real-rooted real polynomial H with
ord_0(H)<k<deg(H),

```text
[s^k]H=0 implies [s^(2k)]H^2<0.
```

If H(0)!=0, the stronger bound is
[s^(2k)]H^2<=-H(0)^2 sum_(|I|=k)r_I^2/(2k-1).
Real-rootedness and the interior condition are essential:1+s^2 and s^2
give the respective hostile boundaries. Repeated real roots are allowed.

Consequently, for f(u)=u^(-a)R(u^d), R real-rooted with nonzero constant
term and exact degreeN, d>=1 and0<a<dN, set m0=d/gcd(a,d). The first
nonzero Laurent moment is exactly m0 or2m0, and is at most the endpoint
widthdN. For N=1 the first alternative always holds. At every congruent
positive mass, a zero moment has a strictly negative doubled moment in the
declared real gauge. Complex scalar/variable gauges preserve vanishing,
not an unnormalized real sign.

For arbitrary nonzero complex AP-trinomial coefficients on support
(-a,d-a,2d-a),0<a<d, the complete row at every admissible m is coprime
to the row at2m. All its roots lie at tau=alpha*gamma/beta^2<0 by
THM-4436; there the ordinary core1+s+tau*s^2 is real-rooted, so the
strict coefficient implication applies. The sharp family
(-2h,1,2h+2),h>=1, has h+1 first channels and h distinct cancellation
choices, each with first detection4h+2 equal to width. Both endpoint
sizes and first-fibre multiplicity are unbounded.

The real-rooted-core method does not cover general primitive trinomials:
their compressed degree must be<=3. It additionally covers the cubic
phase sector -4/27<=tau<0. The scope hostile1+s-s^3 is not real-rooted,
although its first cancellation has a nonzero doubled moment by THM-4432.
General doubled-row separation and the general sharp Laurent bound remain
OPEN outside the proved sectors.

The Gaussian limit recovers, with He denoting probabilists' Hermite,

```text
He_(2k)(sqrt(2)x)=2^k He_k(x)^2
 -sum_(s=1)^k binom(k,s)2^(k-s)(2s-3)!! He_(k-s)(x)^2.
```

Thus at every root of He_k the left side is at most -(2k-3)!!<0.
The scale sqrt(2) is essential, and the identity is a classical Mehler
specialization, not a new literature identity. The finite-subset limit
keeps the exact square coefficient and explains the retained strict margin.

Independent full written audits passed; standard-library normal/optimized
outputs match29,613 explicit gates:4,353 Gram entries,3,829 signed rows,
2,913 tuned zero carriers,360 full AP gcd controls,558 literal real-core
Laurent cases and32 Hermite identities/gcds. The finite tests challenge the
analytic proof; they are not its unbounded justification.
