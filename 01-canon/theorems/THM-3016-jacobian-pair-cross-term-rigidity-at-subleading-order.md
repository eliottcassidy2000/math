---
id: THM-3016
title: "Cross-term rigidity for planar Jacobian pairs at subleading order"
status: >
  PROVED + VERIFIED-EXACT / AWAITING INDEPENDENT HOSTILE AUDIT. Let (P,Q) in
  k[x,y] (char 0) satisfy Jac(P,Q) in k^*, n = deg P, m = deg Q, g = gcd(n,m),
  and let P_n = c H^a, Q_m = c' H^b be the forced leading-form factorisation
  (a = n/g, b = m/g coprime, deg H = g). Put, for a >= b,
  kappa = ca/(c'b), J = Jac(H, Q_{m-1}), and
  W = P_{n-1} - kappa H^{a-b} Q_{m-1}, a form of degree n-1. Then the
  subleading Jacobian condition together with the planar Plucker identity
  forces the single scalar identity  J * Jac(W, H) = 0. Consequently, if
  J != 0 and W != 0 then H is a pure power of a LINEAR form. Equivalently:
  whenever H has at least two distinct roots (K >= 2), the cross term
  Jac(P_{n-1}, Q_{m-1}) is not free -- the pair must satisfy
  Jac(H, Q_{m-1}) = 0  or  P_{n-1} = kappa H^{a-b} Q_{m-1}. Sampled
  automorphisms satisfy the second alternative identically (W = 0) and have
  K = 1. The mirrored statement holds for b >= a with the roles exchanged.
  Strengthening: J = 0 with Q_{m-1} != 0 also forces K = 1, so for K >= 2 the
  dichotomy collapses to  Q_{m-1} = 0 or P_{n-1} = kappa H^{a-b} Q_{m-1}; and
  with W = 0 the cross term equals kappa(a-b)H^{a-b-1}Q_{m-1}J, so every term
  of the degree n+m-4 condition carries an explicit power of H (exponents
  a-1, a-b-1, b-1) and the tower iterates with a Euclidean-reduced pair.
  The Q_{m-1}=0 escape also closes: for K>=2 it forces P_{n-1}=0, so the
  subleading layer is either entirely empty or entirely carried by Q_{m-1}.
  This constrains, but does not decide, JC(2); no bridge or reduction is claimed.
source: death-star-2026-07-31-coinC2
depends_on: []
related:
  - HYP-9070
  - THM-3003
  - THM-3004
external:
  - "Jung (1942), van der Kulk (1953): automorphisms of the affine plane."
script: 04-computation/jc2_cross_term_rigidity_thm3016.py
output: 05-knowledge/results/jc2_cross_term_rigidity_thm3016.out
script2: 04-computation/jc2_cross_term_iteration_thm3016.py
output2: 05-knowledge/results/jc2_cross_term_iteration_thm3016.out
---

# THM-3016 -- cross-term rigidity at subleading order

Throughout `k` has characteristic zero, `(P,Q)` is a **Jacobian pair**
(`Jac(P,Q) = P_x Q_y - P_y Q_x` a nonzero constant), `n = deg P`,
`m = deg Q >= 1`, and `P_j`, `Q_j` denote homogeneous components.

## 1. Setup: the forced leading form

Since `Jac(P_j, Q_l)` is homogeneous of degree `j + l - 2`, the degree
`n+m-2` part of `Jac(P,Q)` is `Jac(P_n, Q_m)` and must vanish. Two binary
forms with vanishing Jacobian are algebraically dependent, so there is a
form `H` and constants with

```text
P_n = c H^a,   Q_m = c' H^b,   deg H = g = gcd(n,m),
a = n/g,  b = m/g,   gcd(a,b) = 1.                                   (L0)
```

The degree `n+m-3` part gives `Jac(P_n, Q_{m-1}) + Jac(P_{n-1}, Q_m) = 0`,
i.e. with `J := Jac(H, Q_{m-1})` and `Y := Jac(P_{n-1}, H)`,

```text
c a H^{a-1} J = - c' b H^{b-1} Y,    so for a >= b:
Y = - kappa H^{a-b} J,      kappa := ca/(c'b).                       (L1)
```

## 2. The planar Plucker identity

For any `A,B,C` in `k[x,y]` the three gradients lie in a rank-2 module, and
Cramer's rule gives the identity (referee check P)

```text
Jac(B,C) grad A + Jac(C,A) grad B + Jac(A,B) grad C = 0.             (P)
```

Apply (P) to `A = H`, `B = P_{n-1}`, `C = Q_{m-1}`, and write
`X := Jac(P_{n-1}, Q_{m-1})` for the **cross term**:

```text
X grad H = J grad P_{n-1} + Y grad Q_{m-1}.                          (2)
```

## 3. The rigidity identity (PROVED)

Substituting (L1) into (2) and using
`H^{a-b} grad Q_{m-1} = grad(H^{a-b}Q_{m-1}) - (a-b) H^{a-b-1} Q_{m-1} grad H`:

```text
[ X - kappa (a-b) J H^{a-b-1} Q_{m-1} ] grad H  =  J grad W,          (3)
W := P_{n-1} - kappa H^{a-b} Q_{m-1}.
```

`W` is homogeneous of degree `n-1`, because
`(a-b) g + (m-1) = (n-m) + m - 1 = n-1`. Taking the determinant of (3)
against `grad H` annihilates the left side, so

```text
J * Jac(W, H) = 0.                                                   (R)
```

**Theorem.** For every planar Jacobian pair, (R) holds. In particular, if
`J != 0` and `W != 0` then `Jac(W,H) = 0`, so the binary forms `W` and `H`
are powers of a common form `G`: `W = alpha G^u`, `H = beta G^v`. Then
`deg G` divides `deg W = n-1` and `deg G` divides `deg H = g`, and `g | n`,
so `deg G` divides `gcd(n-1, n) = 1`. Hence `G` is **linear** and `H` is a
pure power of a linear form.

**Corollary (the cross term is not free).** Let `K` be the number of
distinct roots of `H` (the distinct directions at infinity). If `K >= 2`
then

```text
Jac(H, Q_{m-1}) = 0     or     P_{n-1} = kappa H^{a-b} Q_{m-1}.       (X)
```

For `b >= a` the mirrored statement holds with
`W = Q_{m-1} - kappa' H^{b-a} P_{n-1}` of degree `m-1`, `kappa' = c'b/(ca)`,
and `J' = Jac(H, P_{n-1})`; the divisibility argument is identical since
`g | m`.

## 4. Verification (VERIFIED-EXACT)

The referee checks (P) symbolically, then, on genuine Jacobian pairs built
as composites of affine and triangular maps, verifies `Jac(P_n,Q_m) = 0`,
the (L1) residual `= 0`, and the identity (R):

```text
(x+y^2, y+(x+y^2)^3)     deg (2,6)  g=2 (a,b)=(1,3) H=y^2 kappa=3  W=0  R: 0
deeper composite         deg (12,6) g=6 (a,b)=(2,1) H=y^6 kappa=2  W=0  R: 0
(x+(y+x^2)^4, y+x^2)     deg (8,2)  g=2 (a,b)=(4,1) H=x^2 kappa=4  W=0  R: 0
```

Two observations from the sample, recorded as EVIDENCE and not asserted in
general: every sampled automorphism has `K = 1` (consistent with HYP-9070's
`K = 1` gate), and every one satisfies the **second** alternative of (X)
identically, `W = 0`, even though `K = 1` does not force it.

## 4b. The dichotomy collapses, and the tower iterates (PROVED)

**(A) `J = 0` also forces `K = 1`.** If `J = Jac(H, Q_{m-1}) = 0` with
`Q_{m-1} != 0`, then `H` and `Q_{m-1}` are powers of a common form `G`, so
`deg G` divides `g` and `deg G` divides `m-1`; since `g | m`, `deg G`
divides `gcd(m, m-1) = 1`, so again `H` is a pure power of a linear form.
Hence the first alternative of (X) is unavailable when `K >= 2`, and the
dichotomy **collapses**:

```text
K >= 2   =>   Q_{m-1} = 0    or    P_{n-1} = kappa H^{a-b} Q_{m-1}.   (X')
```

That is: for a counterexample the subleading component of `P` is *determined*
by that of `Q` (or the latter vanishes). The mirrored statement holds with
the roles of `P` and `Q` exchanged.

**(A2) The remaining branch closes too: a clean dichotomy.** Suppose
`K >= 2` and `Q_{m-1} = 0`. Then (L1) reads
`0 = - c' b H^{b-1} Jac(P_{n-1}, H)`, so `Jac(P_{n-1}, H) = 0`. If
`P_{n-1} != 0` this makes `P_{n-1}` and `H` powers of a common form `G` with
`deg G` dividing `n-1` and `g | n`, hence `deg G | gcd(n-1,n) = 1` and
`K = 1` -- a contradiction. So `Q_{m-1} = 0` forces `P_{n-1} = 0`, and
symmetrically. The full statement for a counterexample is therefore the
**dichotomy**

```text
K >= 2  =>   (i)  P_{n-1} = Q_{m-1} = 0,   or
             (ii) both nonzero and  P_{n-1} = kappa H^{a-b} Q_{m-1}.   (X'')
```

In (i) the whole subleading layer is empty; in (ii) it is one-dimensional,
carried entirely by `Q_{m-1}`. Either way the layer has **no free
parameters** beyond `Q_{m-1}` itself.

**(A3) Branch (i) is killed by a translation.** Let `tau(x,y) = (x+s,y+t)`.
Then `Jac(P o tau, Q o tau) = Jac(P,Q) o tau`, so the pair stays Jacobian;
the leading forms are unchanged, so `H, g, a, b, K` are preserved; and the
subleading components shift by the directional derivative of the leading
form (verified symbolically):

```text
P_{n-1} -> P_{n-1} + (s d_x + t d_y) P_n,
Q_{m-1} -> Q_{m-1} + (s d_x + t d_y) Q_m.
```

With `P_n = c H^a`, `Q_m = c' H^b` and `L := (s d_x + t d_y) H` (a form of
degree `g-1`), a branch-(i) pair becomes
`P_{n-1}' = c a H^{a-1} L`, `Q_{m-1}' = c' b H^{b-1} L`, and then

```text
kappa H^{a-b} Q_{m-1}' = (ca/(c'b)) H^{a-b} c' b H^{b-1} L
                       = c a H^{a-1} L = P_{n-1}',
```

so the translated pair satisfies **branch (ii) exactly**. In characteristic
zero `L != 0` for generic `(s,t)` (Euler: `x d_x H + y d_y H = g H != 0` for
`deg H = g >= 1`). Concretely, `(P,Q) = (x, y+x^3)` has
`P_{n-1} = Q_{m-1} = 0`, and after `tau = (x+1,y+1)` it has
`P_{n-1} = 1`, `Q_{m-1} = 3x^2`, still with `Jac = 1`.

**So branch (i) is not a separate case:** after a generic translation every
counterexample with `K >= 2` satisfies

```text
P_{n-1} = kappa H^{a-b} Q_{m-1},   both nonzero,                     (X''')
```

and a pair arriving from branch (i) carries the sharper factorisation
`Q_{m-1} = c' b H^{b-1} L`, `deg L = g-1`. Script:
`04-computation/jc2_kill_branch_i_thm3016.py`.

**(B) With `W = 0` the cross term is explicit and the tower iterates.**
Using `Jac(fg, g) = g Jac(f,g)` (referee-verified) with
`f = kappa H^{a-b}`, `g = Q_{m-1}`:

```text
Jac(P_{n-1}, Q_{m-1}) = kappa (a-b) H^{a-b-1} Q_{m-1} J.              (4)
```

(Verified symbolically at `(a,b) = (4,1), (3,2), (5,2)`.) Substituting (4)
into the degree `n+m-4` condition

```text
c a H^{a-1} Jac(H, Q_{m-2}) + Jac(P_{n-1},Q_{m-1})
                            + c' b H^{b-1} Jac(P_{n-2}, H) = 0
```

shows that **all three terms carry explicit powers of `H`**, namely
`a-1`, `a-b-1`, `b-1`. Dividing by `H^{min(a-b-1, b-1)}` therefore produces
the next relation with a reduced exponent pair -- so the subleading tower
does not stall at the cross term: it iterates, and the exponents that appear
are exactly those of a Euclidean step on `(a,b)`. This answers the open
point flagged in HYP-9070 section 1.

## 5. Scope

(R) is an identity valid for all planar Jacobian pairs and is proved above
from (L0), (L1) and (P) alone; it uses no classification and no
characteristic-`p` input. The corollary (X) is a *constraint on
counterexamples*, obtained by combining (R) with the (classical) fact that
`K = 1` is the one-place-at-infinity situation. **Nothing here decides
JC(2)**, and no bridge, reduction, or equivalence to any other repo lane is
claimed -- cf. MISTAKE-237, which retracted exactly such a promotion. The
natural next questions are whether `J = 0` can be excluded, whether `W = 0`
is forced in general (the sample suggests it is common), and whether the
same wedge argument iterates to the degree `n+m-4` condition.

Referee: `ALL CROSS-TERM CHECKS PASS`. QED.
