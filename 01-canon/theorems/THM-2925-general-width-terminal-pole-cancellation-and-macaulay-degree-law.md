---
id: THM-2925
title: "General-width terminal-pole cancellation and Macaulay degree law"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.  For every
  width M>=1 and tensor order m>=2, each
  mean-difference normalized factorial coefficient is cleared by
  prod_(r=1)^(M-1)(n+r)^(m-1)*(n+M)^(m-2), and the scaled polynomial
  has degree at most (m-1)M-1.  Interior pole control is termwise; at
  the terminal pole the all-top principal part -mC cancels the m
  one-nontop parts +C.  Consequently the canonical 20Q+10C+6F
  degree-seven Macaulay minor has degree at most 58M-36 at every width.
  No uniform minor nonvanishing or arbitrary-width SFC(4) is claimed.
source: codex-gmc-holotopy-extension-2026-07-29
depends_on: []
related:
  - THM-2452-indicator-idempotent-aggregate-endpoint-restoration
  - THM-2697-filtered-affine-handoff-germ-category-and-base-signature-holotopy-boundary
  - THM-2849-four-slot-first-window-macaulay-box
  - THM-2921-diameter-four-nonconsecutive-macaulay-newton-closure
  - THM-2922-diameter-five-macaulay-newton-atlas
  - THM-2924-diameter-six-macaulay-newton-atlas
script: 04-computation/gmc_general_width_terminal_pole_cancellation_thm2925.py
output: 05-knowledge/results/gmc_general_width_terminal_pole_cancellation_thm2925.out
script_sha256: 83d70a95f0943992d0e4b7027eede431d4dc968b66655e37b43fd0acfc692e47
output_sha256: 4fcccbf3998ae1b7249ccae8015b927596b8386f9e6d603ea77246e61f2d3722
constructor_dependency_sha256: 42e9b5ceddd677d1f2601a9d5d668c9437281596b65999ddcb8549d4e0b9bf64
hash_basis: LF-normalized bytes
---

# THM-2925 -- general-width terminal-pole cancellation and Macaulay degree law

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

## 1. Statement

Let

```text
f_j=s^j/j!,                         L(s^j)=j!,          (1)
```

and fix integers

```text
M>=1,                              m>=2.               (2)
```

For an ordered `m`-tuple

```text
0<=b_1,...,b_m<M,                                      (3)
```

form the normalized mean-difference tensor coefficient

```text
A_(b_1,...,b_m)(n)
 = L(prod_(i=1)^m(f_(n+b_i)-f_(n+M)))/L(f_n^m).       (4)
```

Then

```text
D_(m,M)(n) A_(b_1,...,b_m)(n) in Z[n],                (5)

D_(m,M)(n)
 =[prod_(r=1)^(M-1)(n+r)^(m-1)](n+M)^(m-2),          (6)
```

and

```text
deg[D_(m,M) A_(b_1,...,b_m)] <=(m-1)M-1.             (7)
```

The same claims hold after multiplying `(4)` by any multinomial count,
so `(5)--(7)` apply directly to every ordinary-monomial coefficient of
a mean-eliminated factorial moment form.

For the order-two, order-three, and order-four forms in the four-slot
degree-seven Macaulay construction, the scaled row degrees are therefore

```text
M-1,                         2M-1,                    3M-1. (8)
```

The canonical fixed row chart with allocation `20Q+10C+6F` has
determinant degree at most

```text
20(M-1)+10(2M-1)+6(3M-1)=58M-36.                     (9)
```

The four-distinct-slot Macaulay application is substantive for `M>=3`;
the tensor denominator and degree statements themselves include the
boundary widths `M=1,2`.

## 2. The normalized tensor terms

For offsets `0<=e_i<=M`, put `S=sum_i e_i`.  Direct factorial
normalization gives

```text
R_e(n)
 :=L(prod_i f_(n+e_i))/L(f_n^m)
  =(mn+1)_S/prod_i(n+1)_(e_i).                        (10)
```

Expanding `(4)` gives

```text
A_b(n)
 =sum_(J subset {1,...,m})(-1)^|J| R_(e(J))(n),       (11)
```

where `e_i(J)=M` for `i in J` and `e_i(J)=b_i`
otherwise.  Every finite pole of a term in `(10)` lies among

```text
n=-1,-2,...,-M.                                        (12)
```

## 3. Interior poles are controlled termwise

Fix `1<=r<M`.  At `n=-r`, the denominator in `(10)` has order

```text
N_r(e)=#{i:e_i>=r}.                                    (13)
```

If `N_r(e)<=m-1`, the desired pole bound is immediate.  If
`N_r(e)=m`, then every `e_i>=r`, so

```text
S>=mr.                                                  (14)
```

The numerator `(mn+1)_S` therefore contains the unique factor

```text
mn+mr=m(n+r),                                          (15)
```

which cancels one denominator factor.  Thus every individual tensor
term has pole order at most `m-1` at every interior divisor.  No
inclusion-exclusion cancellation is needed here.

## 4. The terminal principal parts cancel

Write

```text
n=-M+x,                         C=(mM-1)!/((M-1)!)^m. (16)
```

Any term of `(11)` with at most `m-2` top offsets already has terminal
pole order at most `m-2`.  Only the all-top term and the `m` terms with
exactly one nontop position can contribute to order `m-1`.

For the all-top term,

```text
(mn+1)_(mM)
 =m(-1)^(mM-1)(mM-1)! x+O(x^2),

(n+1)_M
 =(-1)^(M-1)(M-1)! x+O(x^2).                          (17)
```

Hence its Laurent coefficient before the inclusion sign is

```text
(-1)^(m-1)mC x^(1-m).                                 (18)
```

The all-top inclusion sign is `(-1)^m`, so its contribution to the
principal part of `(11)` is

```text
-mC x^(1-m).                                           (19)
```

Now retain one nontop offset `b<M` and use `M` in the other `m-1`
positions.  Then

```text
S=(m-1)M+b,

(mn+1)_S|_(n=-M)
 =(-1)^S (mM-1)!/(M-b-1)!,

(n+1)_b|_(n=-M)
 =(-1)^b (M-1)!/(M-b-1)!.                            (20)
```

After division by the `m-1` top denominators, the Laurent coefficient
before the inclusion sign is

```text
(-1)^(m-1)C x^(1-m),                                  (21)
```

independent of the retained offset `b`.  Its inclusion sign is
`(-1)^(m-1)`, so each such position contributes

```text
+C x^(1-m).                                            (22)
```

There are exactly `m` positions.  Their total `+mC` cancels `(19)`.
Therefore `(11)` has terminal pole order at most `m-2`.

Before reduction, every term in `(10)` has an integer numerator and a
monic denominator which is a product of the linear factors `n+r`.
Consequently `(11)` can be written with an integer numerator and a monic
denominator having only those factors.  Sections 3--4 bound the exponent
of each factor in its reduced denominator by the corresponding exponent
in `(6)`.  That monic denominator therefore divides `D_(m,M)` already in
`Z[n]`, not merely in `Q[n]`.  This proves the integrality statement
`(5)`.

## 5. Infinity gives the degree bound

Each term in `(10)` is bounded at infinity:

```text
R_e(n) -> m^S,                         n -> infinity.  (23)
```

Thus the signed sum `(11)` is `O(1)`.  By `(5)`,

```text
P(n)=D_(m,M)(n)A_b(n)
```

is a polynomial.  Since `P/D_(m,M)=O(1)`,

```text
deg P<=deg D_(m,M)
       =(m-1)(M-1)+(m-2)
       =(m-1)M-1,                                      (24)
```

which is `(7)`.  Substituting `m=2,3,4` gives `(8)`, and the exact row
allocation gives `(9)`.

## 6. Why this is an aggregate/holotopy theorem

At the terminal divisor, no individual all-top or one-nontop chart is
regular enough: all `m+1` displayed terms have the forbidden
order-`m-1` pole.  Regularity appears only after the full
inclusion-exclusion cospan is recombined.  The cancellation is
address-independent because every one-nontop face has the same principal
coefficient `C`.

This is the factorial-tensor analogue of THM-2452's aggregate-first
restoration and THM-2697's filtered handoff rule: the invariant object is
the recombined section, while an individual transition chart carries a
spurious terminal singularity.  Here the preserved predicate is the
scaled polynomial coefficient; the lost datum is which position retained
which lower offset.

## 7. Scope and the remaining frontier

This theorem proves denominator clearing and degree control for every
width and every tensor order.  It does **not** prove that any selected
Macaulay minor is nonzero.  In particular, `(9)` supplies the exact
finite interpolation budget but not:

```text
Gregory--Newton positivity,
one-chart persistence,
arbitrary-width SFC(4),
a shifted moment window,
SFC(5), or GMC(2).                                     (25)
```

THM-2921 and THM-2922 supply nonvanishing at widths four and five;
THM-2924 separately owns the width-six atlas.  The next uniform target
is to explain the observed shifted-Newton sign law or to replace one
minor by a finite Pluecker atlas when it crosses a chart wall.

## 8. Exact verification

The exact companion hash-pins only the proved ordinary-multinomial
constructor of THM-2921.  It then:

1. computes the terminal Laurent principal parts for
   `1<=M<=12`, `2<=m<=7`, and every `0<=b<M`, checking `468`
   one-nontop controls against `(19),(22)`;
2. proves exact divisor containment and `(7)` for every one of `2912`
   ordinary coefficients on all `56` four-slot shapes with
   `3<=M<=8` and `2<=m<=5`;
3. performs `936` numerical specialization checks against the direct
   normalized moment constructor; and
4. binds the Laurent and polynomial universes by separate SHA-256
   digests.

These finite checks audit the general proof; they do not replace its
quantifiers.  Run

```text
python 04-computation/gmc_general_width_terminal_pole_cancellation_thm2925.py
python -O 04-computation/gmc_general_width_terminal_pole_cancellation_thm2925.py
```

Normal and optimized executions byte-match the stored output with the
declared LF-normalized hashes.

An independent hostile audit rederived the normalized tensor formula,
every interior and terminal Laurent valuation (including repeated lower
offsets), and the monic-division/Gauss-lemma passage from rational pole
bounds to divisibility in `Z[n]`.  A separate SymPy route checked `49`
integer-polynomial controls, `789` interior valuations, `54` terminal
position controls, and `180` multinomial controls, including the boundary
cases `M=1,m=2`, `M=1,m=3`, and `M=2,m=2`.  It also reproduced every
declared count, digest, LF-normalized hash, and the normal/optimized/stored
output equality.  No load-bearing defect remained.

**QED.**
