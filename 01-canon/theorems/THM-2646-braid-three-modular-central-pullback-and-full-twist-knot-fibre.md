---
id: THM-2646
title: "Braid-three modular central pullback and full-twist knot fibre"
status: >
  PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.  The
  exponent/modular map identifies B3 with the fibre product
  PSL2(Z) x_(C6) Z and classifies its conjugacy classes.  The missing integer
  central height is sharp: two shortest knot-closing lifts over one exact
  modular point close to the unknot and trefoil, and one modular fibre
  contains torus knots of unbounded genus.  Every affine V4 lift with the
  standard surjective S3 linear shadow kills the full twist.  The reduced
  Burau numerators along a central fibre obey an exact order-three
  recurrence.  A monic cubic discriminant winding supplies the same integer
  height, but finite S3 inertia plus height is insufficient without the full
  modular braid class.  This classifies fixed three-braids, not knots under
  Markov stabilization, and proves no LRC, G1, JC, DC, or knot-additivity
  conjecture.
source: root-2026-07-28-braid-modular-central-pullback
depends_on:
  - THM-2598-quartic-v4-resolvent-torsor-and-universal-cusp-boundary
  - THM-2632-farey-v4-theta-channel-and-hurwitz-crt-parity-sidecar
related:
  - THM-2455-quartic-swallowtail-scaffold-and-endpoint-corrections
  - THM-2641-modular-abelianization-theta-blindness-and-637-residue-no-go
  - HYP-2033-lrc-annular-braid-center-is-shift
external: "Birman--Menasco, arXiv:0802.1072, Theorem 1 (closed three-braid classification)."
script: 04-computation/braid3_modular_central_pullback_thm2646.py
output: 05-knowledge/results/braid3_modular_central_pullback_thm2646.out
script_sha256: 0466bc89f1be31feb3442de166273baa18c8883e3bda8a6d3f8295d40b8f3a5b
output_sha256: f7f2ab75c5693638ecc80b378a5f415bff01e1d192e1a4306e9a4492128d13e2
hash_basis: LF-normalized bytes
---

# THM-2646 -- the modular hexagon is a quotient of a braid helix

**PROVED + VERIFIED-EXACT + TWO INDEPENDENT HOSTILE AUDITS.**  This theorem
makes the user's binary/ternary co-occurrence literal:
the order-two and order-three modular generators are roots of one central
full twist before that twist is killed.

## 1. The exact central pullback

Use the Artin presentation

```text
B3=<sigma1,sigma2 | sigma1 sigma2 sigma1=sigma2 sigma1 sigma2>.
```

Put

```text
x=Delta=sigma1 sigma2 sigma1,
y=sigma1 sigma2,
z=Delta^2=(sigma1 sigma2)^3.                            (1)
```

Then

```text
B3=<x,y | x^2=y^3>,              z=x^2=y^3,             (2)
```

and `z` generates the centre.  Killing `z` gives

```text
pi:B3 -> <xbar,ybar | xbar^2=ybar^3=1>
       = C2*C3 = PSL2(Z).                                (3)
```

Let `e:B3->Z` be exponent sum, `e(sigma_i)=1`, and let
`chi:PSL2(Z)->C6` be the abelianization normalized by
`chi(pi(sigma_i))=1`.  Thus

```text
e(x)=3, e(y)=2, e(z)=6,       chi(pi(beta))=e(beta) mod 6. (4)
```

The homomorphism

```text
Phi=(pi,e):B3 -> PSL2(Z) x_(C6) Z
             ={(g,n):chi(g)=n mod 6}                     (5)
```

is an isomorphism.  Indeed, if `Phi(beta)=(1,0)`, then
`beta=z^k` and `0=e(beta)=6k`, proving injectivity.  Given a compatible
`(g,n)`, choose any lift `beta_0` of `g`; compatibility makes
`k=(n-e(beta_0))/6` integral and `z^k beta_0` is the unique required lift.

There is also an exact conjugacy version:

```text
beta ~_B3 gamma
 iff pi(beta) ~_PSL2(Z) pi(gamma) and e(beta)=e(gamma).    (6)
```

For the reverse implication, lift a modular conjugator.  The conjugated
braid and `gamma` differ by `z^k`, and equality of exponent sums forces
`k=0`.  Equivalently, the conjugacy classes of `B3` are the compatible pairs

```text
([g],n),             n=chi(g) mod 6.                      (7)
```

## 2. A sharp unknot/trefoil hostile in one modular fibre

Let `a=sigma1 sigma2`.  The braids

```text
beta_-=a^(-1),                  beta_+=a^2=z beta_-        (8)
```

have the same **exact** `PSL2(Z)` image.  Their closures are

```text
hat(beta_-)=T(3,-1)=unknot,     hat(beta_+)=T(3,2)=trefoil. (9)
```

Their symmetric Artin word lengths are `2` and `4`, and the maximum `4` is
sharp among distinct same-`pi` pairs whose closures are knots.  Any such
pair differs by `z^k`, so its exponent sums differ by a nonzero multiple of
six.  A knot closure has three-cycle strand permutation and hence even
exponent sum.  If both word lengths were at most three, both exponents would
lie in `{-2,0,2}`, which is impossible.  Equation (8) attains the bound.  No
uniqueness of the minimizing pair is claimed.

More strongly,

```text
beta_k=a^(3k-1)=z^k a^(-1),       hat(beta_k)=T(3,3k-1).   (10)
```

Every member has the same exact modular image and is a knot because
`gcd(3,3k-1)=1`.  Its genus is

```text
g(beta_k)=|3k-1|-1,                                      (11)
```

so one modular/Farey fibre contains unbounded knot complexity.  This is a
fixed-axis closed-braid statement, not connected sum.  Exponent sum changes
under Markov stabilization and is not by itself a knot invariant.

## 3. The helix over `C6` and the graceful boundary

Equation (5) says that modular abelianization remembers only `e mod 6`.
The integral braid state lies on the `Z`-cover of that hexagon: a helix over
`C6`.  The binary root `x`, ternary root `y`, and common full twist satisfy

```text
2 e(x)=3 e(y)=e(z)=6.                                    (12)
```

After forgetting group action and THM-2632's `L/R`--theta colouring, the
Cayley graph of modular abelianization is an abstract `C6`.  THM-2632 notes
that this uncoloured graph is nongraceful and every edge deletion is
graceful.  The pullback lifts only that uncoloured cycle to a `Z`-line:
cutting it chooses one interval in the cover.  It neither identifies the
abelianization with THM-2632's `S3` permutahedral action nor selects one of
its edges.  Every finite interval is a path; for `n` edges the labels

```text
0,n,1,n-1,2,n-2,...                                      (13)
```

have edge differences `n,n-1,...,1`.  The repair is therefore a choice of
central origin, not an invariant property of the quotient.

This also separates two centres already present in the repository.  The
annular LRC shift in HYP-2033 is conjecturally/experimentally invisible to
loneliness, whereas the ordinary `B3` centre changes (9).  Quotienting a
centre is observable-dependent and requires a literal descent check.

## 4. Why the affine `V4` torsor still cannot see the full twist

Let `V=F2^2`, so `AGL(V)=V4:S3=S4`.  Consider an affine representation of
`B3` whose linear part is the standard **surjective** mod-two shadow

```text
B3 -> PSL2(F2)=GL2(F2)=S3.                                (14)
```

Since `z` is central and its linear part is trivial, its affine image is a
translation by a vector fixed by every element of `S3`.  But

```text
V^(S3)=0,                                                 (15)
```

so every such affine lift kills `z`.  THM-2632's `V4` theta cuts, strand
permutation, and every such equivariant affine lift are therefore
centre-blind.

Surjectivity in (14) is load-bearing.  With

```text
A=[[1,1],[0,1]],       u=(0,1),       g(v)=Av+u,           (16)
```

the affine map `g` has order four.  Sending both braid generators to `g`
obeys the braid relation, has only a `C2` linear image, and sends `z` to
`g^2!=1`.

After marking strands, the three theta channels can label the three
unordered strand pairs.  The six signed band half-twists form
`{+,-} x {three pairs}`, but this is not canonically THM-2632's
`{L,R} x {theta}`: reduction mod two identifies a half-twist and its inverse.
A physical ancestor/writhe orientation is the missing sidecar.

## 5. The exact Burau central recurrence

Use the reduced Burau convention

```text
R(sigma1)=[[-t,1],[0,1]],        R(sigma2)=[[1,0],[t,-t]]. (17)
```

Direct multiplication gives

```text
R(z)=t^3 I.                                                (18)
```

For a fixed braid `beta`, define the unambiguously normalized numerator

```text
N_k(t)=det(I-t^(3k) R(beta))
      =1-t^(3k) tr(R(beta))+t^(6k) det(R(beta)).           (19)
```

Writing `X=t^3`, the entire central fibre obeys

```text
N_(k+3)-(1+X+X^2)N_(k+2)
       +(X+X^2+X^3)N_(k+1)-X^3 N_k=0.                    (20)
```

The three spectral channels are `1,t^3,t^6`.  For a closed three-braid the
fixed Burau normalization

```text
A_k(t)=((1-t)/(1-t^3)) N_k(t)                             (21)
```

represents its Alexander polynomial up to the usual convention.  Equation
(20) is asserted for `N_k` (and hence for this one fixed normalization), not
for arbitrarily unit-renormalized Alexander polynomials.  Split-link zeros
do not affect (20), but prohibit ratio arguments.  At `t=-1`, (18) is `-I`,
so this `SL2(Z)` shadow remembers only the parity of `k`.

## 6. The exact residual closure ambiguity is a flype edge

The following is a **CITED** consequence, not reproved here.  Birman and
Menasco's closed-three-braid classification says that a link of braid index
exactly three has either one `B3` conjugacy class of representatives or
exactly two, the latter precisely for a nondegenerate flype.  By (6), one
compatible pair `([g],n)` already determines one braid conjugacy class and
therefore its closure.  In the reverse direction, the fibre from compatible
pairs to a fixed braid-index-three link has size one or two; the two-point
case is exactly the involution

```text
sigma1^u sigma2^v sigma1^w sigma2^eps
 <-> sigma1^w sigma2^v sigma1^u sigma2^eps.                (22)
```

This is an undirected two-vertex matching, not a tournament.  Lower-braid-
index links form the separate `sigma1^k sigma2^(+/-1)` boundary.  See
J. Birman and W. Menasco, *A note on closed 3-braids*, Theorem 1,
<https://arxiv.org/abs/0802.1072>.

## 7. Discriminant winding is the central-height sidecar

Let a monic squarefree cubic vary around a loop, and follow its three roots.
The resulting root braid `beta` satisfies

```text
e(beta)=wind(prod_(i<j)(r_i-r_j)^2)=wind(Disc).             (23)
```

There is no extra factor two: one positive Artin half-twist rotates one root
difference through `pi`, so its square winds once.  A positive generic
meridian around a divisor `D` therefore has

```text
e(beta_D)=ord_D(Disc);                                     (24)
```

at crossings one uses the full signed intersection-number sum.  For
`f=a prod_i(T-r_i)`, the exact nonmonic formula is

```text
wind(Disc(f))=e(beta)+4 wind(a),
e(beta)=wind(Disc(f))-4 wind(a).                            (25)
```

THM-2598 proves exact equality between the general quartic discriminant and
its **monic integral resolvent** cubic discriminant.  Hence the quartic gives
the integer central height omitted by the finite quotient `S4/V4=S3`.
Once the full modular resolvent-braid conjugacy class is known, (6) and (24)
recover the complete local `B3` conjugacy class.

Two boundaries are essential.

1. Finite `S3` inertia plus exponent is insufficient.  The identity and
   `[sigma1^2,sigma2^2]` both have trivial strand permutation and exponent
   zero but distinct modular images.
2. The resolvent is still a root correspondence, not an affine polynomial
   Keller source.  Nothing here supplies the missing `V4` origin, a Jelonek
   inverse branch, or an `A4/S4` exclusion.

## 8. Transfer contract and exact stopping boundary

The source/target ledger is now explicit.

| source | quotient | preserved | destroyed | required sidecar |
|---|---|---|---|---|
| `B3` | `PSL2(Z)` | modular conjugacy, strand permutation | full-twist height | integer exponent/discriminant winding |
| signed band alphabet | mod-two theta channel | unordered strand pair | crossing sign | integral braid ancestor |
| quartic roots | resolvent matching triple | `S3` matching action, discriminant | `V4` origin | affine/Keller realization |
| closed 3-braid | unmarked link | closure type | braid axis/conjugacy | at most one flype edge at index three |

Thus the user's binary and ternary trees do meet on one object, but only
before quotient: they are two roots and three roots of one central full
twist.  The finite modular hexagon is the quotient after controlled
forgetting.  The exact next tests are to retain the exponent sidecar in the
quartic-resolvent local braid and to ask whether any physical LRC ancestor
identifies Farey `L/R` with signed half-twist orientation.  No such
identification, LRC row exclusion, Markov-invariant knot metric, `G1`,
`JC(2)`, or `DC(2)` conclusion is asserted.

## 9. Exact companion

The companion checks the Burau braid relation and full twist, the radius-
three hostile boundary, the torus-knot fibre, the `V4` fixed subspace, and
the symbolic recurrence (20) with explicit `require` gates:

```text
python 04-computation/braid3_modular_central_pullback_thm2646.py
python -O 04-computation/braid3_modular_central_pullback_thm2646.py
```

Both modes must byte-match
`05-knowledge/results/braid3_modular_central_pullback_thm2646.out`.
