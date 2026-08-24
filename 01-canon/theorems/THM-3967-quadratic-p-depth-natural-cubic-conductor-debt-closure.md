---
id: THM-3967
title: "Quadratic-P-depth natural cubics have no remaining conductor-debt escape"
status: >
  RESERVED / PROVISIONAL PROOF CANDIDATE / VERIFIED-EXACT / INDEPENDENT
  HOSTILE AUDIT REQUIRED. For every irreducible natural cubic
  T^3-3PT-q(P,t) with deg_P(q)<=2, the adjusted hidden polynomial is either
  squarefree, the P^2 zero-section case, or has a repeated polynomial graph
  root. The first case is normal and excluded by the global monogenic
  different; the second is the moving-P^2 normalization of THM-3963; and
  the third is exactly the fully closed THM-3964 family. Hence coefficient
  depth two leaves no same-function-field planar Keller chart. Depth three
  already admits a nongraph repeated factor, although the minimal example
  is reducible; arbitrary higher depth and JC(2) remain open.
source: jc-degree6-one-place / post-THM-3964 coefficient-depth synthesis, 2026-08-24
depends_on:
  - THM-3961-arbitrary-q-hidden-repetition-normality-and-conductor-debt
  - THM-3963-moving-p2-cubic-normalization-principal-ramification
  - THM-3964-polynomial-graph-hidden-double-root-normalization
related:
  - THM-3960-natural-one-parameter-cubic-normal-monogenic-closure
script: 04-computation/jc2_quadratic_p_depth_conductor_closure_thm3967.py
output: 05-knowledge/results/jc2_quadratic_p_depth_conductor_closure_thm3967.out
script_sha256: f27b181754b6425bf8c62e31465102910bcb85a30cd7b2d738aa0146028b934c
output_sha256: 5ce0d901cf18f717263230fbdb4680fd5638d89c7b865531ea4ebf158fee0778
semantic_sha256: 50b543e97c0f3da42d4396d687bce3f55fa61524e99ad274f845102a8137cbd5
hash_basis: raw LF bytes
---

# THM-3967 -- coefficient depth two has no conductor escape

**RESERVED / PROVISIONAL PROOF CANDIDATE / VERIFIED-EXACT / INDEPENDENT
HOSTILE AUDIT REQUIRED.** Work over an algebraically closed field `k` of
characteristic zero. Let

```text
q=a(t)P^2+b(t)P+d(t),
F=T^3-3PT-q,                    a,b,d in k[t],           (1)
```

and assume that `F` is irreducible over `k(P,t)`. Then the finite cubic
field defined by `(1)` admits no same-function-field planar Keller `A2`
chart.

This closes the full coefficient-depth-two slice of the natural globally
monogenic grammar. It does not assert the planar Jacobian conjecture.

## 1. The adjusted hidden polynomial

Put `K0=k(t)`. On the ramification parametrization `P=h^2,T=-h`, the hidden
polynomial is

```text
K(h)=q(h^2,t)-2h^3=a h^4-2h^3+b h^2+d.                 (2)
```

The THM-3961 trichotomy specializes as follows.

1. If `d!=0`, the adjusted polynomial is `L=K`.
2. If `d=0,b!=0`, then

   ```text
   K=h^2L,                    L=a h^2-2h+b.              (3)
   ```

   The displayed `h^2` is the harmless forced factor.
3. If `d=b=0`, then `q=aP^2`. Irreducibility excludes `a=0`, because then
   `F=T(T^2-3P)`. For `a!=0`, THM-3963 computes the full normalization and
   excludes a Keller plane by its principal ramification arm.

In the first two cases, squarefree `L` makes the natural cubic surface
normal, and THM-3961 excludes it by the global monogenic different. It
remains only to classify a repeated irreducible factor `M!=h` of `L`.
Such an `M` also occurs at least twice in `K`.

## 2. A repeated factor cannot be genuinely quadratic

Work in `K0[h]`. Since `deg_h(K)<=4`, a repeated irreducible factor has
degree at most two. If it had degree two, then necessarily `a!=0`, and after
making it monic one would have

```text
M=h^2+s h+p,                    K=aM^2.                  (4)
```

Comparison of the `h^3` and `h` rows in `(2)` and `(4)` gives

```text
2as=-2,                         2asp=0.                  (5)
```

Thus `s=-1/a!=0` and `p=0`, so `M=h(h+s)` is reducible. This contradiction
shows that every repeated `M!=h` is linear over `K0`.

Let its root be `r in K0`, with `r!=0`. From `K(r)=K_h(r)=0` one obtains
exactly

```text
b=3r-2ar^2,                    d=ar^4-r^3.              (6)
```

The full factorization and residual discriminant are

```text
K=(h-r)^2[a h^2+(2ar-2)h+r(ar-1)],
disc(residual)=4(1-ar).                                 (6a)
```

Thus `4ar=3` is exactly the triple-root seam, while a second double root
forces `ar=1` and hence `d=0`. This records the equality boundary rather
than merely locating the repeated component.

If `a=0`, then `b=3r` makes `r` polynomial and

```text
q=3rP-r^3,
F=(T+r)(T^2-rT+r^2-3P),                                (7)
```

contrary to the domain hypothesis. Hence the surviving repeated-root case
has `a!=0`.

## 3. The rational root has no polynomial denominator

Although `(6)` was derived over `K0`, its root is automatically in `k[t]`.
Write it in lowest terms as

```text
r=S/R,                    R,S in k[t],       gcd(R,S)=1. (8)
```

Clearing `(6)` gives

```text
2aS^2-3SR+bR^2=0,
dR^4=S^3(aS-R).                                        (9)
```

The first row and coprimality imply `R` divides `a`; write `a=Ra_1`.
After cancelling one `R`, reduction of the first row modulo `R` gives

```text
2a_1S=3 mod R.                                         (10)
```

The second row becomes `dR^3=S^3(a_1S-1)`, so its reduction modulo `R`
gives `a_1S=1 mod R`. Comparing with `(10)` yields `2=3 mod R`.
Therefore `R` is a unit and `r in k[t]`.

Equations `(6)` now identify the coefficient polynomial exactly:

```text
q=3rP-r^3+a(P-r^2)^2.                                  (11)
```

This is the polynomial-graph repeated-hidden-root family of THM-3964, with
its parameter `c=a!=0`. That theorem computes the regular normalization,
conductor, class lattice, and ramification addresses and excludes every
member from a same-field Keller plane.

The forced-square row has an especially rigid equality case: when
`d=0,b!=0`, the adjusted quadratic in `(3)` repeats exactly when `ab=1`.
Then `a,b` are nonzero constants and `r=1/a`, consistent with `(11)`.

## 4. Exhaustion and conclusion

The cases above are mutually exclusive and exhaustive:

```text
adjusted L squarefree       -> THM-3961 global-different obstruction;
q=aP^2                     -> THM-3963 principal-arm obstruction;
adjusted L not squarefree  -> (11) and THM-3964 boundary obstruction.
                                                               (12)
```

Consequently every irreducible cubic `(1)` is closed in this natural
same-function-field atlas:

```text
deg_P(q)<=2  =>  no same-function-field planar Keller A2 chart. (13)
```

The irreducibility assumption is essential to the organization of the
proof: reducible cubics do not define the required degree-three function
field.

## 5. Sharp next seam: a depth-three nongraph factor

The factor classification changes immediately at the next coefficient
depth. Put

```text
s=P+t,                         q=3sP-s^3.                (14)
```

Then `deg_P(q)=3` and

```text
q(h^2,t)-2h^3
 =-(h-s(h^2,t))^2(2h+s(h^2,t)),                         (15)
```

where

```text
h-s(h^2,t)=-(h^2-h+t)                                  (16)
```

is irreducible over `k(t)` because its discriminant `1-4t` is not a square.
Thus a genuinely nongraph repeated factor already exists at depth three.
Here it is neutralized by the exact factorization

```text
F=(T+s)(T^2-sT+s^2-3P),                                (17)
```

so it fails the domain gate. At depth three a repeated linear graph root
can also leave an affine-in-`P` multiplier after removing `(P-r^2)^2`; its
normalization has a new finite-pole packet and is outside THM-3964. This is
the sharp next conductor-debt problem. General `deg_P(q)>=3`, nonmonogenic
orders, and JC(2) remain open.

## Reproduction

```bash
python3 04-computation/jc2_quadratic_p_depth_conductor_closure_thm3967.py
python3 -O 04-computation/jc2_quadratic_p_depth_conductor_closure_thm3967.py
sha256sum 04-computation/jc2_quadratic_p_depth_conductor_closure_thm3967.py \
  05-knowledge/results/jc2_quadratic_p_depth_conductor_closure_thm3967.out
python3 agents/check_docs.py
```
