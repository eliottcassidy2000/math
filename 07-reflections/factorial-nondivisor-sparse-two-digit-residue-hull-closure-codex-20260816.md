---
title: "Nondivisor residue-hull closure: automatic 3/5 lanes and a sparse two-digit family"
status: >
  PROVED structural lemmas + FINITE-EXACT census mining.  This is a
  theorem-ready no-ID proof packet, not canon.  It proves automatic rho
  admissibility for the p=3 and p=5 congruence lanes and an all-height exact
  Newton-polygon necessary degree barcode for d=p^K+p^s+2 at every odd prime
  p.  Coprimality is proved when the inherited all-divisor packet misses the
  resulting three positive addresses.  It does not prove that this packet
  condition holds for infinitely many parameters and does not extend the
  d<=10000 boundary.
source: /root/jc_level3_divisor_audit/factorial-followup/2026-08-16
depends_on:
  - THM-3483-factorial-nondivisor-residue-digit-pair-compiler
  - THM-3475-factorial-all-divisor-digit-polygon-and-pair-ledger-compiler
  - THM-3161-factorial-newton-euclidean-closure-through-r1998
  - THM-3142-prime-power-predecessor-newton-separation-and-composite-window-census
  - THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census
related:
  - THM-3474-factorial-binary-submask-polygon-and-prime-power-reset-families
script: 04-computation/factorial_nondivisor_sparse_two_digit_hull_followup.py
output: 05-knowledge/results/factorial_nondivisor_sparse_two_digit_hull_followup.out
script_sha256: 2c6cdecbba4ad4d4cd2fbd51d6d4ae1d03bde3a7a136de55b8a9736afb5aaf78
output_sha256: 6c7dcff54ca37f6b7605e303ca64ca2c3635daee4c8f9ea62fea131c005e2c33
hash_basis: LF-normalized bytes
---

# Nondivisor residue-hull closure

## 1. Outcome and exact scope

Use the exact factorial functional and resonant pair

```text
L(x^r)=r!,
A_n^(d)(v)=L((d-x+v x^2)^n),
F=A_(d-2)^(d),                    G=A_(d-1)^(d).             (1)
```

This packet proves two structural results.

1. If `p` is `3` or `5` and `d==2 (mod p)`, both actual Newton polygons in
   (1) equal their THM-3483 raw polygons.  Thus these two congruence lanes
   are automatically rho-admissible at every height.
2. More sharply, for every odd prime `p` and integers `K>s>=1`, put

   ```text
   P=p^K,              Q=p^s,              d=P+Q+2.         (2)
   ```

   The complete Newton-polygon necessary local degree barcode at `p` is
   exactly

   ```text
   D_p(F,G)={0,Q,P,P+Q}.                                    (3)
   ```

The second statement gives an infinite exact local hull/barcode family.  It
becomes a coprimality family under the explicit inherited sidecar

```text
P_div(d) intersect {Q,P,P+Q}=empty,                         (4)
```

where `P_div(d)` is the exact all-divisor necessary degree packet compiled
by THM-3475.  In particular, the cheaper condition

```text
P_div(d) subset (Q,P)                                      (5)
```

closes the row.

The statement is deliberately exact-support quadratic.  It says nothing
about a missing coefficient, translated or arbitrary support, `SFC(3)`, or
multivariable `FC(3)`.  A local allowed degree remains only necessary; it
does not construct a factor.

## 2. Inheritance pass and map

The closest proved mechanism is
`THM-3483-factorial-nondivisor-residue-digit-pair-compiler.md`: for a prime
`p` not dividing `d`, it writes the coefficient height as a raw height plus
the residue correction

```text
rho_p(n,j,d)
 =sum_(ell=0)^m0 C(m0,ell)(-d0^(-1))^ell(2j0+1)_ell mod p,  (6)
```

and proves that nonvanishing at every extreme raw-hull vertex makes the
actual and raw Newton polygons equal.  `THM-3161-factorial-newton-euclidean-
closure-through-r1998.md` supplies the least-used sidecar: when `p|n`, every
raw-hull vertex is congruent to `0` or `h=(p-1)/2` modulo `p`.

The source-to-target map is

```text
(d,p) -> raw weights -> extreme vertices -> rho residues
      -> exact F/G hulls -> common slopes/capacities
      -> local necessary degree set.                       (7)
```

The map forgets factor existence and cross-prime compatibility.  The
retained sidecar is

```text
(p,K,s,vertex coordinates,vertex rho values,
 complete common blocks,P_div(d)).                         (8)
```

The canonical hostile is `d=6518`: it lies in the automatic `p=3` lane but
its raw ledger has a denominator-one block and preserves the complete
divisor packet.  The corrected near miss is therefore "automatic
admissibility implies automatic closure"; it is false.

## 3. Automatic rho admissibility at p=3 and p=5

### Proposition 1

Let `p` be `3` or `5`, and suppose `d==2 (mod p)`.  Then the actual lower
Newton polygons of the two polynomials in (1) equal their raw THM-3483
polygons.

### Proof

Put `n=d-2`, so `n==0 (mod p)`.  By the proved plateau recursion in
THM-3161, every extreme raw-hull vertex of `F=A_n^(d)` has residue

```text
j0 in {0,h},                       h=(p-1)/2.               (9)
```

At `j0=0`, one has `m0=n-j==0`, so (6) is `1`.  At `j0=h`, every
positive-length rising product in (6) begins with

```text
2h+1=p,                                                       (10)
```

and vanishes modulo `p`; again (6) is `1`.  Thus every extreme raw vertex
of `F` is exact.

For `G=A_(n+1)^(d)`, the residue data are `n+1==1` and `d==2`.  Directly
evaluating the finite expression (6) for `j0=0,...,p-1` gives

```text
p=3:       (2,1,2),
p=5:       (3,1,1,2,2).                                  (11)
```

Every entry is nonzero.  THM-3483 now identifies both actual polygons with
their raw polygons.  QED.

This is a theorem at all heights, not an inference from the census.

## 4. The sparse two-unit-digit hull

### Proposition 2

Let `p` be any odd prime, let `K>s>=1`, and define `P,Q,d` by (2).  Put

```text
E_e=nu_p((2p^e)!)=2(p^e-1)/(p-1),
sigma_e=E_e/p^e.                                           (12)
```

Then the complete raw and actual hulls are

```text
F: (0,0), (Q,E_s), (P+Q,E_s+E_K),

G: (0,0), (1,0), (Q+1,E_s), (P+Q+1,E_s+E_K).              (13)
```

Their complete common block ledger is

```text
(slope,capacity,reduced denominator)
 = (sigma_s,Q,Q), (sigma_K,P,P),                           (14)
```

and hence their local necessary degree set is exactly (3).

### Proof of the raw F hull

For a degree `m`, write

```text
w_m(j)=nu_p C(m,j)+nu_p((2j)!).                            (15)
```

For odd `p`, consecutive raw weights satisfy the exact event identity

```text
w_m(j+1)-w_m(j)=nu_p(m-j)+nu_p(2j+1).                     (16)
```

Indeed, the `nu_p(j+1)` from the binomial ratio cancels the same valuation
from `2j+2`.

For `0<=j<=Q`, Kummer's carry formula gives

```text
nu_p C(P+Q,j)=nu_p C(Q,j),                                (17)
```

because the isolated digit at position `K` never participates in
subtracting `j<=Q`.  Hence `w_(P+Q)(j)=w_Q(j)`.  The prime-power estimate in
THM-3142 proves that every internal raw point of `w_Q` lies strictly above
the chord

```text
(0,0)--(Q,E_s).                                           (18)
```

At the final endpoint, the two nonzero digits of `2(P+Q)` do not carry
because `p` is odd, so Legendre gives
`nu_p((2(P+Q))!)=E_K+E_s`.

For the second interval, write `j=Q+t`, `0<=t<=P`.  Summing (16) gives

```text
w_(P+Q)(Q+t)-E_s-w_P(t)
 =sum_(u=0)^(t-1) [nu_p(2Q+2u+1)-nu_p(2u+1)].             (19)
```

This difference is nonnegative.  To see it, expand each valuation as a sum
of divisibility indicators.  At scale `p^e`, the events `p^e|(2u+1)` occur
at

```text
a_e=(p^e-1)/2 mod p^e.                                   (20)
```

The shifted events occur at `a_e-Q`.  If `e<=s`, this is the same residue.
If `e>s`, then

```text
0<=a_e-Q<a_e,                                             (21)
```

so in every period the shifted event occurs earlier.  Every prefix
`0<=u<t` therefore contains at least as many shifted events as unshifted
events, proving (19).

The prime-power estimate for `w_P(t)` is strict above its endpoint chord
for `0<t<P`.  Equation (19) therefore puts every internal second-interval
point strictly above

```text
(Q,E_s)--(P+Q,E_s+E_K).                                   (22)
```

Finally,

```text
sigma_s=2(1-p^(-s))/(p-1)
       <2(1-p^(-K))/(p-1)=sigma_K,                        (23)
```

so the two chords form the complete lower hull of `F`.

### Proof of the raw G hull

Since `p|(P+Q)` and `p` is odd, the binomial and factorial ratios give

```text
w_(P+Q+1)(j+1)=w_(P+Q)(j)+nu_p(2j+1).                    (24)
```

The right side is at least the shifted `F` height.  Equality holds at the
three `F` vertices `j=0,Q,P+Q`, since all three are divisible by `p` and
`2j+1` is a unit.  Thus the positive-index part of the `G` hull is the
one-step right shift of the `F` hull.  Its constant point adds the horizontal
edge `(0,0)--(1,0)`, proving the raw formulas (13).

### Residue exactness and degrees

At every displayed `F` vertex, `n0=j0=0`, so (6) equals `1`.  For `G`, the
vertex `j0=0` has residue `1-1/2=1/2`, while every other displayed vertex
has `j0=1` and `m0=0`, so its residue is `1`.  These values are nonzero for
every odd `p`.  THM-3483 therefore promotes the raw hulls to the actual
hulls.

The numerator of `sigma_e` is

```text
2(1+p+...+p^(e-1))==2 (mod p),                            (25)
```

so its reduced denominator is exactly `p^e`, equal to its block capacity.
A common-factor contribution must use either zero or the whole of each common
block.
Taking subset sums proves (3).  The unit constant vertices show that the
common coordinate-root capacity is zero.  QED.

## 5. Packet closure corollary

For every prime `q|d-1`, let `D_q(F,G)` be the complete THM-3475 local
necessary degree set and define

```text
P_div(d)=intersection_(q|d-1) D_q(F,G) intersect {1,...,d-2}. (26)
```

Any nonconstant rational common factor of `F,G` has one global degree, and
that degree must lie in every divisor-place set in (26) and in the
nondivisor set (3).  Therefore (4), or the stronger easy-to-read condition
(5), implies

```text
gcd_Q(F,G)=1.                                               (27)
```

By the exact resonance reduction in
`THM-3124-quadratic-factorial-moment-recurrence-and-shifted-window-census.md`,
(27) excludes three consecutive zero moments for the corresponding
exact-support quadratic window.

This corollary is a genuine all-parameter theorem with an explicit packet
hypothesis.  What is not yet proved is that (4) holds for infinitely many
pairs `(K,s)` for any fixed `p`.

## 6. Exact census readout and positive control

The three frozen THM-3483 blocks contain exactly `281` post-divisor rows.
Their final killer histogram is

```text
p=3:31, p=5:23, p=7:81, p=11:83, p=13:31,
p=17:18, p=19:9, p=23:4, p=29:1.                          (28)
```

Mining those rows gives

```text
lane                    rows    final killer p    other final label
d==2 (mod 3)              66          31                 35
d==2 (mod 5)              38          23                 15.          (29)
```

For `p=5`, three of the last column were already killed at `p=3`; the table
records final labels, not counterfactual scan order.  Proposition 1 proves
admissibility independently of this bookkeeping.

Exactly one of the `281` rows has the two-unit-digit form (2) for `p=3` or
`5`:

```text
d=6590,                 p=3, K=8, s=3,
d-2=3^8+3^3=6588.                                         (30)
```

Proposition 2 gives

```text
D_3(F,G)={0,27,6561,6588}.                                (31)
```

The frozen all-divisor packet is

```text
P_div(6590)={599,1198,1797,2396,2995,
             3594,4193,4792,5391,5990},                   (32)
```

which lies wholly in `(27,6561)`.  Thus (30) is an exact positive control
for the packet closure theorem, and the new proof explains its `p=3`
killer without coefficient construction.

The statement that (30) is the only sparse row is `FINITE-EXACT` over the
declared `281` rows.  Propositions 1 and 2 are `PROVED` at all heights.

## 7. Hostiles and sharp boundaries

### Congruence alone does not close: d=6518

The unique `p=29` census row is

```text
d=6518==2 (mod 3),
P_div(6518)={3087,3430,4802,5145,5488,5831}.               (33)
```

Its `p=3` residue test is automatically admissible, exactly as Proposition
1 predicts.  But its common raw ledger contains

```text
(8/9,9), (26/27,27), (80/81,81), (242/243,243),
(728/729,729), (2186/2187,2187), (1,3240),                 (34)
```

where each pair denotes `(slope,capacity)`.  The final denominator-one
block fills the would-be central region, and every degree in (33) survives
`p=3`.  The later primes `13,17,19` retain `5831`, `p=23` is inadmissible,
and `p=29` finally maps `{5831}` to the empty set.  This is the canonical
hostile to deleting the digit-shape/packet sidecar.

### Digit collision is load-bearing

The proof of Proposition 2 assumes `K>s`.  At `p=3,K=s=2`, one has
`d-2=18=2*3^2`, and the common blocks are

```text
(8/9,9), (1,9).                                            (35)
```

The resulting local degree set is every integer from `0` through `18`, not
just `{0,9,18}`.  The two isolated unit digits and the absence of a carry in
`2(P+Q)` are therefore load-bearing.

### The tempting 5q lane is false without an interval/digit sidecar

Taking `d-1=5q`, `q==2 (mod 3)`, does put the row in the automatic `p=3`
lane, and the large divisor can restrict degrees to multiples of `q` in
favorable one-digit cases.  It is not a theorem by congruence alone.  The
minimal cheap hostile

```text
d=56, q=11, d-2=54=2*3^3                                (36)
```

has common blocks `(26/27,27)` and `(1,27)` and retains all four candidates
`11,22,33,44`.  Likewise, the first `p=5` lane hostile

```text
d=22, q=7, d-1=3q
```

has common blocks `(2/5,10)` and `(1/2,10)` and retains both divisor
candidates `7,14`.

These examples locate the first failed implication:

```text
automatic rho admissibility
  -/-> a central raw-hull gap
  -/-> empty intersection with P_div(d).                   (37)
```

## 8. Rigorous stopping reason and next test

The present work does not claim an unconditional infinite coprimality
family from (2).  The missing assertion would be

```text
P_div(p^K+p^s+2) intersect {p^s,p^K,p^K+p^s}=empty        (38)
```

for infinitely many `(K,s)`.  Neither the base-`p` two-digit identity nor
the congruence `d==2 (mod p)` controls the prime factorization of
`d-1=p^K+p^s+1`.  The hostiles (33) and (36) show that suppressing this
sidecar would be mathematically false, not merely unaudited.

The cheapest next structural test is therefore not a larger `d` scan.  It
is a divisor-shape theorem: find a proved infinite subfamily of
`p^K+p^s+1` with one THM-3475 divisor barcode contained in `(p^s,p^K)`, or
prove that a specified algebraic factor of this number forces such a
barcode.  Until that factorization-to-barcode map is proved, (4) is the
sharp theorem-safe closure hypothesis.

## 9. Reproduction

Run

```text
python3 04-computation/factorial_nondivisor_sparse_two_digit_hull_followup.py
python3 -O 04-computation/factorial_nondivisor_sparse_two_digit_hull_followup.py
```

Both modes are byte-identical to
`05-knowledge/results/factorial_nondivisor_sparse_two_digit_hull_followup.out`.
The companion hash-pins the THM-3483 engine and all three stored census
outputs.  It checks the two rho tables, all `281` killer labels and packets,
`59` sparse hull cells over

```text
(p,K_max)=(3,8),(5,6),(7,5),(11,4),                       (39)
```

the `d=6590` positive control, the `d=6518` full hostile trace through
`p=29`, the digit-collision boundary, the `d=56` five-times-prime hostile,
and the `d=22` three-times-prime hostile.  The semantic digest is

```text
97c1967457cb9744b3bc04e04f0a9880d9b000041efc17161171fd67a9d3dfea. (40)
```

The finite grid is a hostile control for the proof; it is not the basis of
the all-height propositions.
