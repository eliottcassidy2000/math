---
id: THM-3702
title: "Universal equal-step three-by-four Danielewski nonentry"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For
  every squarefree polynomial Sigma over C of degree at least two, no Darboux
  pair on c^2e=Sigma(b) can have grading supports {-5,-2,1} and
  {-5,-2,1,4}.  Endpoint commutation peels fifth, second, and fourth powers;
  the remaining middle row forces both weight -2 coefficients to have arm
  order at least two, so every scalar address vanishes.  The other three
  3-subsets of the four-weight window shear to previously closed cells.  This
  is all-degree and strictly subsumes the finite nine-function THM-3700 gate;
  arbitrary 3x4 support words and JC(2) remain open.
source: jc-sparse-direct-search / 2026-08-22
audit: >
  PASS -- root independently checked the bracket convention and all six
  fibres, the d=0 shear landing in <=2 x 3, exact UFD power exponents, the
  high-row sign, the middle primitive and polynomial identity, squarefree-arm
  unit division and all scalar orders, the other three subset shears, and the
  THM-3700 finite-span corollary.  Normal and optimized replays, hashes, and
  documentation checks pass.
depends_on:
  - THM-3569-danielewski-two-by-three-weight-darboux-nonentry
  - THM-3592-universal-exponent-two-three-by-three-weight-darboux-nonentry
  - THM-3695-y0-collision-ring-danielewski-embedding-and-seven-piece-floor
related:
  - THM-3696-y0-collision-ring-three-branch-conductor-and-graded-modules
  - THM-3699-y0-consecutive-four-weight-three-by-four-nonentry
  - THM-3700-y0-equal-step-two-weight-adjunction-pluecker-no-go
script: 04-computation/jc2_universal_equal_step_three_by_four_nonentry_thm3702.py
output: 05-knowledge/results/jc2_universal_equal_step_three_by_four_nonentry_thm3702.out
script_sha256: 8d215f37d9762eefacce2b182acf307004a3deb25a1a6482b5c5899a56797898
output_sha256: a15b0cf2519e933b085be677d688bd0f217862fd2a5b9b8260c10720ecc80284
semantic_sha256: 6bf3f14b9e3617d015a55e1ae8c138fc858b71938a2dc142e9c55cb70a0b190b
hash_basis: raw LF bytes for files; bucket table, differential peel identities, middle-row primitive, arm jets, and source-chart bracket for semantic hash
---

# THM-3702 -- universal equal-step three-by-four Danielewski nonentry

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**
This is the first all-degree closure inside the step-three four-weight window
suggested by THM-3700.  The proof does not depend on the special collision
polynomial `1-b^2`; it works in every squarefree exponent-two Danielewski
ring.  Its mechanism is a two-sided root peel followed by one exact
middle-row integral.

## 1. Universal setting and statement

Work over `C`.  Let `Sigma in C[b]` be squarefree of degree at least two and
put

```text
D_Sigma=C[b,c,e]/(c^2e-Sigma(b)),
wt(b,c,e)=(0,1,-2).                                      (1)
```

On `c!=0`, a homogeneous weight-`u` element has the unique form `c^u f(b)`.
It is regular on `(1)` precisely when

```text
Sigma^ceil(-u/2) divides f                 for u<0.       (2)
```

The Poisson coefficient rule is

```text
{c^u f,c^v g}=c^(u+v+1) W_(u,v)(f,g),
W_(u,v)(f,g)=v f'g-u f g'.                              (3)
```

After scalar weight-zero pieces have been removed, suppose nonzero
homogeneous pieces give

```text
P=c^-5 f+c^-2 p+c r,
Q=c^-5 g+c^-2 q+c s+c^4 t.                              (4)
```

Thus

```text
supp(P)={-5,-2,1},                 supp(Q)={-5,-2,1,4},
Sigma^3|f,g,                       Sigma|p,q.             (5)
```

Then

```text
{P,Q} != 1.                                               (6)
```

The same conclusion holds after exchanging the two outputs.  More generally,
no `3 x 4` Darboux pair has both supports in the window

```text
W={-5,-2,1,4}.                                           (7)
```

## 2. The six-fibre grammar

For the hard support orientation in `(5)`, the bracket output weights and
addresses are

```text
-9: (-5,-5)
-6: (-5,-2), (-2,-5)
-3: (-5,1),  (-2,-2), (1,-5)
 0: (-5,4),  (-2,1),  (1,-2)
 3: (-2,4),  (1,1)
 6: (1,4).                                               (8)
```

Hence the five nonscalar rows vanish and the middle three-address row equals
one.  The fibre sizes are `1,2,3,3,2,1`; the proof peels inward from the two
singleton endpoints.

## 3. Low endpoint: a fifth/second-power root

The weight `-9` row is

```text
W_(-5,-5)(f,g)=5(fg'-f'g)=0.                            (9)
```

Thus `g=lambda f` for some `lambda in C*`.  The adjacent weight `-6` row
then compresses exactly to

```text
W_(-5,-2)(f,q)+W_(-2,-5)(p,lambda f)
      =W_(-5,-2)(f,q-lambda p)=0.                      (10)
```

Put `d=q-lambda p`.  The case `d=0` is not silently divided away.  If it
held, the Darboux shear

```text
(P,Q) -> (P,Q-lambda P)                                (11)
```

would preserve the bracket, cancel the weights `-5` and `-2` from the second
output, and leave at most two of its weights against the three weights of
`P`.  This is the `2 x 3`-or-smaller cell excluded by THM-3569.  Therefore
`d!=0`.

For nonzero polynomials, `(10)` reads

```text
5 f d'=2 f'd.                                           (12)
```

Taking valuations at every irreducible polynomial in the UFD `C[b]`, and
using `gcd(5,2)=1`, gives nonzero `alpha,kappa in C` and `k in C[b]\{0}`
such that

```text
f=alpha k^5,                    d=kappa k^2,
q=lambda p+kappa k^2.                                  (13)
```

This is the exact common-power conclusion, including its exponents; no
root or coefficient has been discarded.

## 4. High endpoint: a fourth-power peel

The weight `6` row is

```text
W_(1,4)(r,t)=4r't-rt'=0,
```

so for `mu in C*`,

```text
t=mu r^4.                                               (14)
```

Substitute `(14)` in the adjacent weight `3` row.  Directly from `(3)`,

```text
W_(-2,4)(p,mu r^4)+W_(1,1)(r,s)
             =W_(1,1)(r,s-4mu p r^3)=0.                (15)
```

Because the two remaining weights are both one, `(15)` says that its second
argument is a scalar multiple of `r`.  Thus for `nu in C`,

```text
s=4mu p r^3+nu r.                                      (16)
```

## 5. The middle nonscalar row integrates once

Only the weight `-3` equation remains before the scalar row.  Substitute
`(13)` and `(16)` into

```text
W_(-5,1)(f,s)+W_(-2,-2)(p,q)+W_(1,-5)(r,g)=0.          (17)
```

Exact differentiation gives the identity

```text
left side of (17)
 =k^4 d/db [ 5alpha(nu-lambda)kr
             +20alpha mu kpr^3
             -2kappa p/k^2 ].                          (18)
```

Since `C[b]` is a domain and `k!=0`, the rational function in brackets has
zero derivative and is a scalar `rho in C`.  Multiplying by `k^2` yields the
polynomial identity

```text
(20alpha mu k^3r^3-2kappa)p
       =rho k^2-5alpha(nu-lambda)k^3r.                  (19)
```

This is the decisive gain: the nominally simple weight-`-2` coefficient is
forced to inherit two copies of the endpoint root.

## 6. Every scalar address dies on every arm

Let `beta` be a root of `Sigma`.  Squarefreeness gives
`ord_beta(Sigma)=1`.  From `(2)`, `(5)`, and `(13)`,

```text
3 <= ord_beta(f)=5 ord_beta(k),
```

so

```text
ord_beta(k)>=1.                                         (20)
```

The coefficient multiplying `p` on the left of `(19)` has value
`-2kappa!=0` at `beta`; it is a unit in the local ring.  The right side is
divisible by `k^2`.  Hence

```text
ord_beta(p)>=2,
ord_beta(q)>=2,                       by q=lambda p+kappa k^2.            (21)
```

The scalar row is

```text
W_(-5,4)(f,t)+W_(-2,1)(p,s)+W_(1,-2)(r,q)=1.           (22)
```

At `beta`, the three negative coefficients `f,p,q` have orders at least
`5,2,2`, respectively.  Differentiation can lower an order by at most one,
while `r,s,t` are regular polynomials.  Every summand on the left of `(22)`
therefore vanishes at `beta`.  Its value is zero, contradicting the value one
on the right.  This proves `(6)` for every squarefree `Sigma`.

## 7. The whole four-weight `3 x 4` window

There are four three-element subsets of `W`.  The hard subset in `(5)` omits
the top weight `4` and is closed by Sections 2--6.  Every other subset
contains `4`.  If it is the three-weight support of `P` while `Q` uses all of
`W`, the unique highest bracket row is the address `(4,4)`.  Its vanishing
forces the two weight-four coefficients to be proportional.  A linear
Darboux shear cancels the weight-four piece of one output, leaving `3 x 3`
or smaller support.  THM-3592 and THM-3569 exclude those cases.  Exchanging
the outputs changes only the bracket sign.  This proves the window claim
after `(7)`.

As a concrete corollary, consider the nine-dimensional span in THM-3700:

```text
(C,BC^2,AC^3; A,B^2,A^2C^2,ABC; A^2B,C^4),             (23)
```

whose weight profile is `{-5,-2,1,4}` and whose weight-four space is
one-dimensional.  The universal seven-piece floor leaves only `3 x 4`,
`4 x 3`, or `4 x 4`.  The first two are now closed.  In the last, a shear in
the common one-dimensional weight-four space reduces to at most `4 x 3`,
also closed.  Thus the all-degree argument strictly subsumes the finite
Pluecker no-go in THM-3700, while THM-3700 remains an independent exact
certificate.

## 8. Exact reproduction and frontier

Run

```bash
python3 -B 04-computation/jc2_universal_equal_step_three_by_four_nonentry_thm3702.py
python3 -O -B 04-computation/jc2_universal_equal_step_three_by_four_nonentry_thm3702.py
```

Both streams must agree byte-for-byte with the stored output.  The companion
checks the bucket table, every differential identity in `(9)--(19)`, generic
arm jets for `(22)`, and the bracket formula independently on the literal
`Sigma=1-b^2` collision source chart.

This theorem closes one absolute support word and its four-weight window, not
an arbitrary `3 x 4` placement.  It constructs no Darboux pair and no planar
Jacobian counterexample.  **QED.**
