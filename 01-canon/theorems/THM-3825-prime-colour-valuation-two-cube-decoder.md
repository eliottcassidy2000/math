---
id: THM-3825
title: "Prime-colour valuations decode inert-scaled two-cube states"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  On the carrier
  where the primitive pair sum is inert and cube-free and the common scale
  is inert, quotient and remainder of every inert valuation modulo three
  recover the scale and primitive shell; the split valuations recover the
  Eisenstein cofactor.  One square test then recovers the coprime cube pair.
  Multiplication by an inert cube is the native scale law.  A 3-adic tag
  gives unordered and orientation-refined natural-number addresses for the
  5,855 support-two ratios.  On those tagged atlases, squarehood reads an
  exact label-parity/orientation bit on the unique square primitive address,
  while triangularity cuts out 31 states in either gauge.  Those states are
  fundamental points on two Pell sheets and never return to the pure-3
  ordinate under a later positive Pell power.  These are fixed label-order
  gauges, not LRC(14) certificates.  The separately proved
  THM-3818 recovers arbitrary scale inside the finite atlas and adds the
  selected facet, word, and residue-lattice packet.
source: root + lrc_graver_address + cube_decoder_audit / incoming-signal extension session, 2026-08-23
audit: >
  INDEPENDENT HOSTILE AUDIT PASS (cube_decoder_audit, 2026-08-23).  The audit
  rederived the valuation decoder, found the orientation collision in the
  proposed 78-channel covector address, verified the repaired 156-channel
  tag, and froze ramified-prime, split-scale, exponent-three, relabelling,
  and singleton-failure boundaries.  The assertion-free companion has
  2,066,834 active gates; normal and optimized raw LF streams match its
  frozen transcript.  A separate assertion-free 1,382,012-gate classical-
  number sidecar exhausts all 456,690 unordered and 913,380 oriented tags,
  proves the square parity rules, finds the two exact 31-state triangular
  censuses, and rewrites them as fundamental points on two Pell discriminant
  sheets whose later positive powers exit the pure-3 ordinate.  It also finds
  the unique triangular label-successor edge and checks that no tagged address
  is square-triangular.  Normal and optimized runs byte-match its frozen
  output.  An independent tag-sidecar audit (tag_sidecar_audit, 2026-08-23)
  rebuilt the atlas with `sympy.factorint`, rescanned every tag with `gmpy2`,
  verified all 31 fundamental solutions by an earlier-ordinate search, audited
  the Lucas no-return proof and the negative-Pell successor, and found no scope
  or transcript discrepancy.
depends_on:
  - THM-463-two-cube-representations-are-a-divisor-property-on-the-split-axis
  - THM-3793-inert-prime-sum-all-scale-two-cube-singleton
related:
  - THM-3743-lonely-runner-polyhedron-khinchin-flatness-relation-reduction
  - THM-3818-scaled-inert-cubeclass-support-two-pair-packet
script: 04-computation/two_cube_prime_colour_decoder_thm3825.py
output: 05-knowledge/results/two_cube_prime_colour_decoder_thm3825.out
script_sha256: a3e2a1950a2c684e55cbb2c15ff3f8ea90d2b7dad042a589c077dc6ef27e971e
output_sha256: e002a15b000ad8bd00011da54ac5b1b5b4f18ffdca20f68af66f5e44970b9e8b
tag_sidecar_script: 04-computation/lrc14_prime_colour_square_triangular_tags_thm3825.py
tag_sidecar_output: 05-knowledge/results/lrc14_prime_colour_square_triangular_tags_thm3825.out
tag_sidecar_script_sha256: e6e4f82e90a727b73ba9621f602eea03ab79bbe204303a372605b6c7525ecdc6
tag_sidecar_output_sha256: 8386c3c30dd106cb7b7a70c6bfccf906e025dd01f422d5e3cb2cd37f40c25d49
tag_sidecar_semantic_sha256: 57d484bf6f5e9c84ac174206bf48564c222e4646965c060bbbcdd540ade621db
hash_basis: raw LF bytes
---

# THM-3825 -- inert valuations are scale digits plus shell digits

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**  This is an
arithmetic decoder theorem and a typed natural-address theorem.  Its finite
LRC corollary does not assert loneliness, arrival, or a bad row.

## 1. Restricted carrier and exact decoder

Let `a<b` be coprime positive integers, put

```text
s=a+b,                       q=a^2-ab+b^2,             (1)
```

and assume

```text
every p|s satisfies p=2 mod 3 and v_p(s)<=2.           (2)
```

Let `g>=1` have only prime divisors congruent to `2 mod 3`, and define

```text
m=g^3(a^3+b^3)=g^3 s q.                               (3)
```

Then the prime factorization of `m` recovers `(g,a,b)` without divisor
enumeration.  For a positive integer `M` with `3` not dividing `M`, set

```text
G(M)=product_(p=2 mod 3) p^floor(v_p(M)/3),
S(M)=product_(p=2 mod 3) p^(v_p(M) mod 3),
Q(M)=product_(p=1 mod 3) p^v_p(M).                    (4)
```

On `(3)`, exactly

```text
(G(m),S(m),Q(m))=(g,s,q).                             (5)
```

Moreover `(4)` gives an exact image test for this restricted carrier.  Put

```text
D=(4Q(M)-S(M)^2)/3.                                  (6)
```

Then `M` has the form `(3)` if and only if `D` is the square `delta^2` of an
integer satisfying

```text
0<delta<S(M),             delta=S(M) mod 2.           (7)
```

In that case the unique decoded state is

```text
g=G(M),
a=(S(M)-delta)/2,          b=(S(M)+delta)/2.           (8)
```

The parity clause in `(7)` is automatic once `(6)` is an integer square, but
it is retained in the executable decoder as a hostile gate.

### Proof

THM-463 applied to the primitive pair `(a,b)` says that `q` is odd, every
prime divisor of `q` is `1 mod 3` except possibly `3`, and

```text
v_3(q)=1 iff 3|s.                                     (9)
```

Hypothesis `(2)` makes `3` divide neither `s` nor `q`; the inert scale `g`
is also 3-free.  Hence the inert and split factors of `(3)` are disjoint and

```text
v_p(m)=3v_p(g)+v_p(s)        for p=2 mod 3.           (10)
```

Because `0<=v_p(s)<=2`, quotient and remainder in `(10)` prove `(5)`.
Finally

```text
4q-s^2=3(b-a)^2,                                      (11)
```

which proves the forward square test.

Conversely, `(4)` gives `M=G^3SQ` by unique factorization.  Equations
`(6)--(8)` give positive distinct integers with sum `S`, and

```text
a^2-ab+b^2=(S^2+3delta^2)/4=Q.                       (12)
```

If a prime divided both `a` and `b`, it would divide `S` and its square
would divide `Q`; this is impossible because the prime supports of `S` and
`Q` are disjoint.  Thus `(a,b)` is primitive, and `(12)` proves `(3)`.
This proves the image test and uniqueness inside the restricted carrier.

THM-3793 supplies the stronger global consequence

```text
r_+(m)=1,                                             (13)
```

so no positive distinct two-cube representation outside the restricted
carrier competes with the decoded one.

## 2. Natural-number tags and their operation law

The bare value `m` is already a sparse injective natural address for
`(g,a,b)`.  To attach a fixed unordered coordinate pair `1<=i<j<=13`, put

```text
kappa(i,j)=binom(j-1,2)+(i-1) in {0,...,77},
A_un=3^kappa(i,j) m.                                  (14)
```

Since `3` does not divide `m`, `v_3(A_un)` recovers `kappa`, hence `{i,j}`;
stripping that power and applying `(4)--(8)` recovers `(g,a,b)`.  At
primitive scale this gives exactly

```text
78*5,855=456,690                                      (15)
```

unordered labelled states.

An unordered coordinate pair does **not** say which coordinate receives the
smaller speed.  Retain an orientation bit `epsilon`, where `epsilon=0`
assigns `(ga,gb)` to `(i,j)` and `epsilon=1` assigns `(gb,ga)`.  The repaired
tag is

```text
lambda=2kappa(i,j)+epsilon in {0,...,155},
A_or=3^lambda m.                                      (16)
```

This recovers the assigned pair and therefore the primitive support-two
covector, up to global sign:

```text
epsilon=0:  b e_i-a e_j,
epsilon=1:  a e_i-b e_j.                              (17)
```

There are `913,380` primitive oriented placements.  Both `(14)` and `(16)`
use the chosen coordinate order as a gauge; neither is invariant under
`S_13` relabelling.

The native scale operation is closed.  If `c` has only inert prime divisors,

```text
(g,a,b)->(cg,a,b)       sends m->c^3m,
A_un->c^3A_un,          A_or->c^3A_or.                (18)
```

Thus the arithmetic address carries inert dilation exactly.  Dense ordinal
successor and Stern--Brocot branching are different operations and are not
asserted to descend through `(14)` or `(16)`.

### 2.1 Square and triangular behaviour of the label operation

The incoming label construction itself has a classical-number signature.
Every primitive base `m=a^3+b^3` in this atlas is 3-free.  THM-3818's
independent classical-number sidecar proves that exactly one primitive base is
a square:

```text
56^3+65^3=671^2.                                     (18a)
```

Consequently coprimality of `3` and `m` turns squarehood into an exact parity
test.  Across all unordered tags `(14)`,

```text
A_un is square iff (a,b)=(56,65) and kappa is even,   (18b)
```

giving `39` square states.  Across the oriented tags `(16)`,

```text
A_or is square iff (a,b)=(56,65) and epsilon=0,       (18c)
```

giving `78` square states.  Thus, only on the unique square primitive fibre,
ordinary squarehood reads the orientation bit from the oriented address.  It
does not decode orientation on any other fibre.

Triangularity behaves differently.  Exhaust the declared tag universes by the
integer discriminant test

```text
3^ell m=T_N  iff  1+8*3^ell*m=(2N+1)^2.              (18d)
```

Split `ell=2r+eta`, with `eta in {0,1}`, and put
`X=2N+1,Y=3^r`.  Then `(18d)` is exactly the two-sheet Pell equation

```text
X^2-8*3^eta*m Y^2=1,             Y=3^r.              (18e)
```

For the oriented tag, `(r,eta)=(kappa,epsilon)`: the orientation bit selects
the discriminant `8m` or `24m`, while the coordinate label constrains the Pell
`Y`-coordinate to a pure power of three.  For the unordered tag, the same
split applies to the parity of its single label exponent.  This is an exact
reframing of the finite scan, not a claim that the pure-3 solutions are closed
under Pell multiplication.

In fact, continued-fraction reduction proves that all `31` displayed points
are the **fundamental** positive solutions on their respective Pell sheets.
There is a general no-return reason.  Write their positive Pell powers as

```text
X_n+Y_n sqrt(D)=(X+Y sqrt(D))^n,
Y_n=Y U_(n-1)(X),                                     (18f)
```

where `U_0=1,U_1=2X,U_(j+1)=2XU_j-U_(j-1)`.  For every `n>=2`, the factor
`U_(n-1)(X)` is not a power of three:

- if `n` is even, it is even;
- if `n` is odd and `3|X`, it is `(-1)^((n-1)/2) mod 3`;
- if `n` is odd, `3` divides neither `X` nor `n`, then
  `X=+1` or `-1 mod 3` gives `U_(n-1)(X)=n mod 3`; and
- in the remaining case `3|n`, the Lucas divisibility law makes
  `U_2(X)=4X^2-1` divide `U_(n-1)(X)` (immediate from the factorization of
  `(t^n-t^(-n))/(t-t^(-1))` when `3|n`).  If `U_2(X)` were a power of three,
  its coprime factors `2X-1,2X+1` would be powers of three differing by two,
  forcing `X=1`, whereas every displayed point has `X>1`.

Because `Y` itself is a power of three, `(18f)` proves that no later positive
Pell point has a pure-3 ordinate.  In particular the next one has
`Y_2=2XY`, the immediate exit checked by the companion.  Thus the tag-selected
triangular points are isolated in their fixed Pell cyclic groups, not a return
suborbit.

There are exactly `31` triangular unordered tags and exactly `31` triangular
oriented tags.  In both scans the numerical exponent/base hits coincide and
have exponent histogram

```text
ell:       0  1  2  3  4  5  6
count:     8  9  7  3  1  2  1.                     (18g)
```

Exactly one numerical hit has its immediate label successor in the same
triangular census:

```text
9^3+13^3=T_76,             3(9^3+13^3)=T_132.        (18h)
```

This is the fourth nontrivial point on the classical triangular-tripling
torsor.  Indeed

```text
T_M=3T_N  iff  (2M+1)^2-3(2N+1)^2=-2,                (18i)
```

and `(M,N)=(132,76)` gives `265^2-3*153^2=-2`.  Multiplication by
`2+sqrt(3)` generates the positive chain from `(5,3)`, corresponding to
`(M,N)=(2,1)`.  The uniqueness in `(18h)` is only inside the finite tag
atlas.  Its combinatorial meaning depends on gauge: `ell:0->1` moves from
unordered label `kappa=0` to `kappa=1`, but on the oriented address it keeps
`kappa=0` and flips `epsilon:0->1`.

This equality of counts is not equality of labelled states: in `(14)` the
exponent is `ell=kappa`, whereas in `(16)` it decodes as
`kappa=floor(ell/2), epsilon=ell mod 2`.  The oriented triangular set splits
`17/14` between `epsilon=0/1`.  The eight `ell=0` hits are precisely the
primitive triangular cube addresses of THM-3818; the remaining 23 are created
by the 3-adic label operation.  The full 31-row list is frozen in the sidecar
output.  None of these tagged addresses is square.  These are
**FINITE-EXACT** arithmetic subatlases, with no loneliness, arrival, or bad-row
consequence.

## 3. Sharp losses and relation to the larger packet

The first orientation hostile is already

```text
A_un=28 on labels {1,2}:
(1,3) gives covector (3,-1),
(3,1) gives covector (1,-3).                           (19)
```

The two covectors are not equal up to sign.  This is why `(16)` needs its
extra bit.  The fixed gauge also sends address `28` on `{1,2}` to `252` on
`{2,3}` under relabelling.

The arithmetic scope gates have small witnesses:

```text
ramified shell:       9=1^3+2^3,
ramified scale/tag:   756=3^3(1^3+3^3),
split scale:          9604=7^3(1^3+3^3),
exponent-three loss:  152=3^3+5^3, shell 8=2^3.       (20)
```

The first high-inert-exponent singleton failure is

```text
65728=12^3+40^3=31^3+33^3,       shell 64=2^6,       (21)
```

and the first failure with primitive exponent exactly three is

```text
515375=15^3+80^3=54^3+71^3,      shell 125=5^3.      (22)
```

THM-3818 is separately **PROVED**.  It uses finite cube-class separation to
handle arbitrary positive common scale inside the atlas and, after a selected
labelled placement, recovers the support-two covector, exposed facets,
THM-778 word, and an exact ambient pair-sum-lattice covering test.  THM-3825
instead gives a table-free factorization decoder under the stricter inert-scale
hypothesis.  Even its oriented arithmetic tag does not reconstruct the
off-lattice residue schedule, endpoint owner, root, phase, arrival, other
eleven speeds, or a loneliness predicate.  Neither theorem excludes an
LRC(14) row.

## 4. Exact verification

The assertion-free companion directly checks:

- all `94` admissible shells and `5,855` primitive pairs with `a+b<=356`;
- all inert scales through `100` and both `456,690`/`913,380` tag counts;
- the decoder and orientation round trips; and
- the hostile minima in `(19)--(22)`, including ramified- and split-scale
  singleton failures.

Run

```bash
python3 -B 04-computation/two_cube_prime_colour_decoder_thm3825.py
python3 -B -O 04-computation/two_cube_prime_colour_decoder_thm3825.py
```

Both raw LF streams equal
`05-knowledge/results/two_cube_prime_colour_decoder_thm3825.out`.

The square/triangular tag sidecar is reproduced by

```bash
python3 -B 04-computation/lrc14_prime_colour_square_triangular_tags_thm3825.py
python3 -B -O 04-computation/lrc14_prime_colour_square_triangular_tags_thm3825.py
```

Both streams equal
`05-knowledge/results/lrc14_prime_colour_square_triangular_tags_thm3825.out`.
