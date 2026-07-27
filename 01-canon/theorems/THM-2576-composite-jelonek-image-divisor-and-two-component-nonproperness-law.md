---
id: THM-2576
title: "Composite Jelonek image divisor and two-component nonproperness law"
status: >
  PROVED over C, with the image-divisor computation exact over Q.  For
  dominant polynomial maps the nonproperness set of a composite is proved to
  be the outer nonproperness set union the outer image of the inner one.  For
  the fixed sporadic Keller map, the closure of F(V(L)) is the irreducible
  degree-25 hypersurface V(H), where H is defined by an exact saturated
  resultant.  Consequently S_(F o F)=V(LH) has exactly two irreducible
  components.  The raw resultant's a^8 c^18 S^8 factors are inverse-chart
  artifacts, not composite discriminant multiplicities.  No discriminant,
  divisor-scheme, higher-iterate component-count, JC(2), or GMC(2) conclusion
  is claimed.
source: codex-2026-07-28
depends_on:
  - THM-2473-sporadic-keller-branch-tower-depressed-trisection-anatomy
  - THM-2570-jelonek-cusp-cylinder-normalization-and-conductor
related:
  - HYP-9033-discriminant-tower-and-genus-axis-of-the-keller-monoid
  - MISTAKE-287 (raw pullbacks and resultants require coefficient sidecars)
script: 04-computation/keller_composite_jelonek_thm2576.py
output: 05-knowledge/results/keller_composite_jelonek_thm2576.out
script_sha256: 8de7ebee702d641eb25b15159c48ba312db3e6972bf23f5bb7cd307d175a9727
output_sha256: b12cd9a780bb3baad72dd08be4dfa019f28f6f6bf221ac218a3f29dba4158577
---

# THM-2576 -- the sporadic square has exactly two Jelonek components

Let the asymptotic or nonproperness set of a polynomial map
`h:C^n -> C^n` be

```text
S_h={y : some x_k has ||x_k|| -> infinity and h(x_k) -> y}.  (1)
```

For the fixed sporadic Keller map `F`, THM-2473 proves

```text
S_F=V(L),
L=27a^2c^2-18abc+16a+b^3c-b^2.                    (2)
```

This theorem first proves the correct general set law for composition and
then computes the previously unnamed second component for `F o F`.

## 1. The nonproperness set of a composite

Let `f,g:C^n -> C^n` be polynomial maps and suppose that `g` is dominant.
Then

```text
S_(f o g) = S_f union f(S_g).                       (3)
```

This equality does not require a closure on its right side.

For the forward inclusion, take `x_k -> infinity` with
`f(g(x_k)) -> y` and put `z_k=g(x_k)`.  If a subsequence of `z_k` is
unbounded, then `y in S_f`.  Otherwise a subsequence converges to some `z`;
then `z in S_g` and `y=f(z) in f(S_g)`.

The inclusion `f(S_g) subset S_(f o g)` is immediate from (1).  For
`S_f subset S_(f o g)`, choose `z_k -> infinity` with `f(z_k) -> y`.
Dominance of `g` implies that its image contains a nonempty Zariski-open set
`U`.  Choose `w_k in U` arbitrarily close to `z_k`, close enough that

```text
||f(w_k)-f(z_k)|| < 1/k,       ||w_k-z_k|| < 1.     (4)
```

Choose `x_k` with `g(x_k)=w_k`.  Since `w_k -> infinity`, polynomial
boundedness on compact sets forces `x_k -> infinity`; now
`f(g(x_k))=f(w_k) -> y`.  This proves (3).

The set (1) is closed: if `y_j -> y` with `y_j in S_h`, choose one witness
`x_j` with `||x_j||>j` and `||h(x_j)-y_j||<1/j`.  Thus the apparently
nonclosed right side of (3) is automatically closed.  Iterating (3) gives
the exact set-level tower law

```text
S_(F^m) = union_(j=0)^(m-1) F^j(S_F)               (5)
```

for every dominant polynomial self-map `F`.  Equation (5) does not assert
that the displayed images are distinct or irreducible.

## 2. Normalize the first Jelonek component before taking its image

THM-2570 gives a set-bijective normalization of `V(L)` by the affine plane.
With parameters `(tau,lambda)`, write

```text
nu(tau,lambda)=(X,Y,Z),

X=lambda^2(3-tau lambda)/27,
Y=lambda(4-tau lambda)/3,
Z=tau.                                              (6)
```

Thus

```text
F(V(L)) = Psi(C^2),             Psi=F o nu.         (7)
```

The first two rows of `dPsi` have determinant `4/3` at
`(tau,lambda)=(1,0)`.  Hence the irreducible closure of (7) has dimension
two and is a hypersurface in target `C^3`.  The remaining task is to identify
its prime equation without confusing inverse-coordinate artifacts with image
components.

## 3. A small inverse norm defines the 361-term equation

Use target coordinates `(a,b,c)`.  For a source point `(X,Y,Z)` in a fibre,
put `s=XY`.  THM-1310's conic pair is

```text
C1 = 3aX^2-(1+s)(bX-s),
C2 = aX^3+c(1+s)^3-X(1+s)(2+s).                    (8)
```

Its exact degree-one subresultant in `s` is `D(X)s+N(X)`, where

```text
D = -3X^2ac+X^2b^2c-X^2b+2Xbc-2X+c,

N = -3X^3abc+4X^3a-6X^2ac+X^2b^2c-X^2b
    +2Xbc-2X+c.                                     (9)
```

The terminal subresultant is

```text
Res_s(C1,C2)=aX^3 E(X),
E(X)=LX^3+(4-3bc)X-2c.                             (10)
```

On the honest inverse chart, `s=-N/D`, while the third component of `F`
gives

```text
Y=s/X,                 Z=[X(2-3s)-c]/X^3.          (11)
```

Substitution into the *source* copy of the Jelonek polynomial gives
`L(X,Y,Z)=P(X,s)/X^6`, with

```text
P = 16X^7+296X^4s^2-360X^4s+108X^4
    +180X^3cs-108X^3c+27X^2c^2
    -3Xs^4+2Xs^3-cs^3.                              (12)
```

Define the polynomial

```text
Q(X)=D(X)^4 P(X,-N(X)/D(X)).                        (13)
```

It has `X`-degree `15` and `398` terms.  Finally define `H in Q[a,b,c]`
by the exact factorization

```text
Res_X(E,Q) = -a^8 c^18 S^8 H,
S=27ac^2-9bc+8.                                    (14)
```

The quotient in (14) is polynomial.  Exact factorization over `Q` proves
that `H` is primitive and irreducible, with

```text
multidegree_(a,b,c)(H) = (14,21,12),
total degree(H)        = 25,
number of terms        = 361,
gcd(H,L)               = 1.                        (15)
```

To specify the large polynomial independently of display order, list the
coefficients of `Poly(H,a,b,c).terms()` in SymPy's lexicographic monomial
order, one `monomial:coefficient` pair per line.  Its SHA-256 is

```text
b146c11f33e895c08f72303d282e2668d955e0a58d9268a1b445d4d5202016c2. (16)
```

Equations (8)-(14), rather than the hash, are the mathematical definition of
`H`; (16) is a transport checksum.

## 4. The image closure is exactly `V(H)`

On the dense inverse chart of `V(L)`, a point `p=(X,Y,Z)` and its target
`q=F(p)` satisfy both `E(q;X)=0` and `Q(q;X)=0`.  Hence their resultant
vanishes.  The chart is genuinely nonempty: the normalized point
`(tau,lambda)=(1,1)` has

```text
X(1+XY) a(q) c(q) S(q) D(q;X) != 0.                 (17)
```

Therefore (14) gives `H(F(p))=0` on a dense open subset of `V(L)`, hence on
all of its image closure.  Section 2 shows that this irreducible closure has
dimension two, while (15) says that `V(H)` is an irreducible hypersurface.
Consequently

```text
closure(F(V(L))) = V(H).                            (18)
```

The companion supplies an orthogonal pointwise control: it reconstructs `H`
from (14) and checks `H(Psi(tau,lambda))=0` at `25` exact hostile parameter
pairs.  It also checks the rank-two minor.  At `(tau,lambda)=(1,0)`,

```text
Psi(1,0)=(1,0,0),       H(1,0,0)=0,       L(1,0,0)=16, (19)
```

which cleanly separates the new image component from the old one.

The factors `a^8`, `c^18`, and `S^8` in (14) are therefore artifacts of the
chosen inverse chart, not additional image components.  Algebraically,

```text
(Res_X(E,Q)) : (acS)^infinity = (H).                (20)
```

This is the same saturation lesson as MISTAKE-287 at the next compositional
level.

## 5. Exact nonproperness set of the square

Apply (3) with `f=g=F` and use THM-2473:

```text
S_(F o F)=V(L) union F(V(L)).                       (21)
```

The left side is closed by the diagonal argument after (4).  Taking the
closure of the right side and using (18) therefore gives

```text
S_(F o F) = V(L) union V(H) = V(LH).               (22)
```

Because `L` and `H` are distinct irreducibles, (22) has exactly two
irreducible components.  This proves the set/component part of HYP-9033 P1
for the first composite grade.

## 6. Boundary of the result

- Equation (22) is a reduced set equality.  It does not determine the
  scheme structure or multiplicity of either component in a composite
  eliminant or discriminant.
- The exponents `8,18,8` in (14) belong to an inverse-coordinate resultant.
  They are not evidence for the predicted odd discriminant part `L H`.
- Equation (5) is an all-depth set law, but distinctness and irreducibility of
  the higher images `F^j(V(L))` remain open.  No grade-`m` component-count
  theorem follows from the `m=2` calculation.
- No statement here concerns a general Keller divisor, proves `JC(2)`, or
  changes the remaining planar-Jacobian branches.

## Reproduction

```bash
python3 04-computation/keller_composite_jelonek_thm2576.py
python3 -O 04-computation/keller_composite_jelonek_thm2576.py
diff -u 05-knowledge/results/keller_composite_jelonek_thm2576.out \
  <(python3 04-computation/keller_composite_jelonek_thm2576.py)
sha256sum 04-computation/keller_composite_jelonek_thm2576.py \
  05-knowledge/results/keller_composite_jelonek_thm2576.out
```

The ordinary and optimized runs reconstruct the exact factorization over
`QQ`, refactor `H`, verify its coefficient ledger, and replay all controls.
