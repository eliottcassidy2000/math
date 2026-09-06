# Width-two constant-term returns are parabolic petal cycles

**Status: CONDITIONAL during independent review.** The argument below is a
complete paper proof candidate, with its classical dynamical input checked in
the primary lecture notes. It is not a scarce-ID reservation or canon entry.
The exact algebraic audit is **FINITE-EXACT**. External priority is not claimed.

## Inheritance and concept board

- Closest proved mechanism: [THM-2111, effective compound-root bound](../../01-canon/theorems/THM-2111-effective-compound-root-bound-for-one-variable-constant-terms.md),
  especially its logarithmic small-root-product identity and exact first-return
  contact. Its current general estimate is `binom(M+N,min(M,N))`.
- Hostile: [THM-2070, horizontal Wick embedding and dihedral cancellation](../../01-canon/theorems/THM-2070-horizontal-wick-embedding-and-dihedral-return-cancellation.md)
  has reachable lengths whose coefficient sums vanish. Here the concrete
  polynomial `z^-2+z^-1+z-z^2` has first nonzero return `4`, with value `-12`.
- Corrected near miss: [THM-1650, Newton polygon](../../01-canon/theorems/THM-1650-newton-polygon-of-the-effective-dvdk-bound.md)
  refutes the claim that sparse endpoints always maximize the first return.
  MISTAKE-211 forbids splitting one scalar equation into separate atom equations.
- Forgotten sidecar: [THM-1595, Dickson ladder](../../01-canon/theorems/THM-1595-tnc-dickson-ladder-closures.md)
  retains the two-small-root factor; the present operation uses its root-swap
  involution instead of successively eliminating coefficients. Its old sample
  claims are not used as proved inputs.

Board: small-root product/contact; local root-swap involution; rational-map
iteration; parabolic critical budget; trinomial coefficient-phase gcd.
The first four now form an exact connection. The fifth produced a bounded
positive signal only, subsequently subsumed by incoming work as recorded below.

## Statement

Let `R in C[z]`, `r=R(0)!=0`, and `d=deg R>=3`. Define

```text
f(z)=z^-2 R(z),                  T(z)=-r z/R(z),
m*=min{m>=1: CT(f^m)!=0},
D_R=(R-zR')/gcd(R,R').
```

The gcd can have any nonzero constant normalization. Let `c_R` be the
number of distinct zeros of `D_R`, and `s_R` the number of distinct zeros
of `R`. Then the claimed bound is

```text
m* <= c_R <= s_R <= d.                                  (1)
```

In particular, the width bound `m*<=M+N` holds whenever `min(M,N)=2`,
by replacing `z` with `z^-1` if needed. No coefficient sign, reality,
squarefreeness, sparsity, or genericity assumption is present.

Together with THM-2111's width-one case, this resolves the entire
`min(M,N)<=2` region, subject to the independent review status above.
It leaves `min(M,N)>=3` open in this argument.

## 1. Local involution from exactly two small roots

Put `h(z)=z^2/R(z)`. Since `r!=0`, choose a holomorphic square root near
zero and set `w(z)=z/sqrt(R(z))`. This is a local coordinate and `h=w^2`.
Consequently the holomorphic germ

```text
iota(z)=w^(-1)(-w(z))
```

is an involution, with `iota'(0)=-1` and `h(iota(z))=h(z)`.
For sufficiently small nonzero `t`, the two small roots of

```text
Phi(X,t)=X^2-tR(X)
```

are precisely `z,iota(z)` when `t=h(z)`. If `Pi(t)` denotes their product,
the identity in THM-2111 §2 is

```text
log(Pi(t)/(-r t)) = sum_(m>=1) CT(f^m)t^m/m.             (2)
```

Existence of `m*` follows from THM-2111; only its existence and (2), not its
effective bound, are used here. Write `c=CT(f^m*)!=0`. Then

```text
Pi(t)=-r t(1+(c/m*)t^m*+O(t^(m*+1))).                  (3)
```

Substitute `t=h(z)` and divide `z*iota(z)=Pi(h(z))` by `z`. Since
`-r h(z)/z=T(z)=-z+O(z^2)`, it follows that

```text
iota(z)=T(z)-b z^(2m*+1)+O(z^(2m*+2)),
b=c/(m* r^m*) !=0.                                    (4)
```

Thus `T=iota+b z^k+O(z^(k+1))` with odd `k=2m*+1`. Composing twice gives

```text
T(T(z))
 =iota(iota(z))+iota'(iota(z))*b z^k+b*iota(z)^k
   +O(z^(k+1))
 =z-2b z^k+O(z^(k+1)).                                (5)
```

Here `iota^2=id`, `iota'(iota(z))=-1+O(z)`, and
`iota(z)^k=-z^k+O(z^(k+1))`. In particular the leading terms add; they do
not cancel. The exact relation is

```text
ord_0(T^2-id)=2m*+1,
[z^(2m*+1)](T^2-z)=-2 CT(f^m*)/(m* r^m*).              (6)
```

This is the missing transfer from a Laurent moment to an iterated map.

## 2. Cited dynamical input and its precise scope

**CITED:** John Milnor, *Dynamics in One Complex Variable: Introductory
Lectures*, Stony Brook IMS Preprint 1990/5, partially revised 1991-09-05,
[primary lecture-note PDF](https://abel.math.harvard.edu/archive/118r_spring_05/docs/milnor.pdf),
§7.2 (flower theorem), §7.10 (a critical point in each immediate petal
basin), and §10.3 plus the sentence immediately following it (printed
page `10-1`, PDF page 66). These statements were read on 2026-09-05.

The imported conclusion is: for a rational map of degree at least two,
a parabolic fixed point with primitive multiplier of order `q`, whose
`q`th iterate has fixed-point multiplicity `n+1`, has `n/q` petal cycles
and at least `n/q` distinct critical points in its immediate basin.
Points whose orbits land exactly on the parabolic point are excluded from
these basins, as stated in §7.4 and §7.5.

The same petal-cycle critical assertion is independently stated in the
introduction, page 2, of Adam Epstein,
[*Infinitesimal Thurston Rigidity and the Fatou-Shishikura Inequality*](https://arxiv.org/pdf/math/9902158),
1999. The more general refined Fatou-Shishikura theorem is unnecessary.

For the present `T`, degree is exactly `d`: its numerator and denominator
are coprime, since `R(0)!=0`. The fixed point `0` has multiplier `-1`,
so `q=2`; (6) gives `n=2m*`. Therefore its basin contains at least `m*`
distinct critical points of `T`.

## 3. The available critical-point budget

At a finite nonpole,

```text
T'(z)=-r (R(z)-zR'(z))/R(z)^2.                         (7)
```

Neither infinity nor any pole belongs to a parabolic attracting basin:

```text
{poles of T} -> infinity -> 0 -> 0.
```

These orbits land exactly on the parabolic fixed point. All critical
points supplied by the cited theorem must therefore be nonpole zeros of
`R-zR'`.

If `alpha!=0` is a root of `R` of multiplicity `e`, then `R-zR'` has
exact multiplicity `e-1` there: writing `R=(z-alpha)^e U` leaves leading
coefficient `-alpha e U(alpha)!=0`. Thus dividing by `gcd(R,R')` removes
exactly all pole roots, with their multiplicities. Its remaining zeros
are the critical points permitted in the parabolic basin.

Since `deg(R-zR')=d` and `deg gcd(R,R')=d-s_R`,

```text
deg D_R=s_R,
number of distinct eligible critical points = c_R <= s_R. (8)
```

Combining (8) with the lower bound `m*` from §2 proves (1).

The coarser Riemann-Hurwitz count gives the same `d` bound: infinity has
ramification multiplicity `d-2`, leaving `d` of the total `2d-2`.
Repeated poles remove another `d-s_R`. Equation (7) makes this counting
entirely explicit and avoids relying on that alternate computation.

## Boundaries, equality, and transfer ledger

- **Nonvacuity/sharpness in unbounded degree:** for odd `d>=3`,
  `f=z^-2+z^(d-2)` has first return `d`, with constant term `binom(d,2)`.
  This attains the bound for every odd width. The canonical four-term
  hostile above attains `d=4`. No universal equality classification is claimed.
- If `m*=d`, both `R` and `D_R` must be squarefree. Every eligible critical
  point must lie in a different petal-cycle basin. This is a necessary
  extremal condition, not a characterization.
- The omitted one-sided family `R=1+z` gives an exact Mobius involution
  `T=-z/(1+z)` and all positive constant terms vanish. It tests why the
  two-sided/degree hypothesis cannot simply be erased.
- Characteristic zero and complex dynamics are used. No finite-field
  version or Lean build is asserted.
- Source -> target: width-two Laurent polynomial -> rational map `T`.
  Map: `R=z^2f`, `T=-R(0)z/R`. Preserved predicate: first nonzero moment
  equals half the parabolic iterate-defect order minus one half, by (6).
  Lost information: this integer alone forgets higher moments and phase;
  the full rational map and leading coefficient in (6) are the sidecars.
  Cheapest decisive test: compare exact Laurent convolution with exact
  rational composition, including cancellation and repeated-root cases.
- The two-small-root product singles out the other root through
  `iota=Pi(h(z))/z`. For three or more small roots that quotient is a
  product of several roots and is not a local root-swap map. This is the
  first missing implication for wider negative support.
- In [THM-2022](../../01-canon/theorems/THM-2022-gmc2-frobenius-lowest-balanced-face.md),
  this replaces the lowest-face seed bound by `M+N` whenever the exact
  Laurent face has width two on one side. The subsequent good prime is
  still coefficient dependent; no uniform final Gaussian moment bound follows.

## Exact reproduction and bounded side lane

```bash
python3 04-computation/nc2_width_two_parabolic_seed_synthesis_sep05.py
python3 04-computation/nc2_trinomial_seed_probe_synthesis_sep05.py --height 35
```

The first script checks (6) by direct Laurent convolution and a separate
exact rational-composition polynomial. All 480 choices with degrees
`3..6`, endpoint coefficients `+-1`, zero constant Laurent coefficient,
and all other coefficients in `{-1,0,1}` pass. Maximum first returns by
degree are `3,4,5,6`. Additional controls include nonunit endpoints,
nonzero first moment, odd-binomial equality through degree 13, repeated
roots through degree 10, and an excluded exact one-sided involution.
Both scripts use explicit failure gates rather than Python assertions;
optimized replays retain all checks and reproduce the saved output exactly.

The side-lane trinomial probe covers all 17,761 primitive triples
`{-a,b,c}` with `1<=a,c<=35` and `1<=b<c`. Exactly 553 have multiple
monomials in the least reachable constant term. In every one, the first
and doubled-return coefficient polynomials have gcd a monomial, excluding
a common nonzero parameter root. This is **FINITE-EXACT**, not a uniform
proof of the forgotten THM-1680/1705 doubling conjecture. It supplies no
additional all-pattern theorem here. The larger incoming exact census below
subsumes this universe; the script is an independent replay/provenance artifact,
not a distinct new census claim.

Matching output artifacts are `nc2_width_two_parabolic_seed_synthesis_sep05.out`
and `nc2_trinomial_seed_probe_synthesis_sep05.out` in this directory. SHA-256:

```text
7eaa6ae69467252271e224ff6e5081d81e8ef97b8ba81860252f32aa1284a783  width-two script
89b7efc2c582aa7420653e7f14680e7d12bbe6a999f119b9d1ed11836d4b4085  width-two output
61977162398f44ed748114e77680e6d0407f2a7dfc1f6e5896eb38932d280af2  trinomial script
947c443ee00f53a69ebcf72200bf3a68107a3952689218854ff828960d9f1e1d  trinomial output
```

Repository duplication searches included `Sturmfels`, `compound-root`,
`min(M,N)=2`, the Dickson ladder, parabolic/Fatou/iteration paired with
Laurent constant terms, and first-return corrections. No existing
parabolic-map bound was found. That is repository recovery evidence,
not a claim of external novelty.

## Incoming-work integration: commit `566677ae1d`

Before filing this candidate, the live remote supplied
[`synthesis_20260905_moments_trinomial.md`](synthesis_20260905_moments_trinomial.md).
That report proves an all-exponent support classification and an exact
all-level carry profile for primitive trinomials. For charges `(-a,b,c)`
with `0<b<c`, write

```text
g=gcd(a+b,a+c), A=(a+b)/g, B=(a+c)/g.
```

Its first support return is collided exactly when `a-AB` belongs to
the numerical semigroup generated by `A,B`; in that case its mass is `g`
and `2g<=a+c`. Its all-trinomial two-rung coprimality assertion is still
OPEN, with a complete exact census through `a,b,c<=60`. The present
height-35 trinomial calculation lies entirely within that larger universe
and supplies no additional supported range.

This clarifies the overlap with the new dynamical argument. In a primitive
trinomial with sole negative charge `-2`, collision forces `2>=AB`.
Because `1<=A<B` are coprime, this means `A=1,B=2`, hence

```text
(-a,b,c)=(-2,b,2b+2), b odd,
g=b+2, a=AB, first support return g, width 2g.
```

The incoming freeness criterion places this entire collided family in
[THM-2639, equal-mass two-rung certificate](../../01-canon/theorems/THM-2639-gmc-equal-mass-two-rung-persistent-collision-certificate.md).
Thus the width-two trinomial bound already follows from that existing
two-rung theorem plus the new incoming classifier. A negative middle charge
gives the reflected pattern `(-c,1,2)`, whose consecutive sums have `g=1`
and a singleton first return; a neutral middle charge returns at mass one.
We do not count those trinomial cases as a separate new result.

The present advance covers arbitrary many nonzero terms at width two and
also proves the stronger distinct-critical-point/distinct-root bounds (1).
The incoming classifier neither gives those bounds nor subsumes the
rational-involution argument. Conversely, our width-two theorem does not
classify arbitrary-width trinomial support, recover the incoming carry
coordinates, or prove its all-trinomial two-rung conjecture.

The exact connection is therefore a scope intersection with two proof
mechanisms: channel freeness supplies two explicit polynomial equations
on sparse faces; rational dynamics bounds delayed contact on all width-two
faces. It is not an equivalence between those mechanisms.
