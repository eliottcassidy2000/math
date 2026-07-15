---
id: THM-774
title: Folded-diamond obstruction for the two-sheet equality packet
status: PROVED (elementary sharp cap; exact independent audits; finite-exact low-core closure)
source: codex-2026-07-14-S3
depends_on:
  - THM-769
related:
  - THM-772
  - THM-773
  - THM-775
  - HYP-6820
verification:
  - 04-computation/lrc13_folded_diamond_two_sheet_obstruction_codex_S3.py
  - 04-computation/lrc13_s2_unbounded_odd_core19_closure_codex_S3.py
---

# THM-774 — Folded-diamond obstruction

This theorem is a metric refinement of THM-769's two-exception equality
packet.  Write

```text
A = 2U union {x,y},       |U|=10,       x>y positive odd,
a = (x+y)/2,              b = (x-y)/2.
```

Put

```text
G_U = {tau in R/Z : phi_U(tau)>1/13}.
```

## 1. The exact folded diamond

For every `tau in R/Z`, the two eligibility and opposite-colour conditions of
THM-769,

```text
||x tau|| <= 2/13,       ||y tau|| <= 2/13,
epsilon_x(tau) != epsilon_y(tau),
```

are together equivalent to the single closed inequality

```text
||a tau|| + ||b tau|| >= 11/13.
```

There is also a clearance-level form of the same fold.  If

```text
B_(x,y)(tau)
  = max_(j in {0,1}) min(||x(tau+j)/2||,||y(tau+j)/2||),
Q_(a,b)(tau) = ||a tau||+||b tau||,
```

then, without any eligibility hypothesis,

```text
B_(x,y)(tau) = (1-Q_(a,b)(tau))/2.                         (1a)
```

Thus `Q>=11/13` says exactly that neither of the two lifts can have both odd
runners clearer than `1/13`.

Consequently `A` is tight if and only if

```text
G_U subset H_(a,b),
H_(a,b) = {tau : ||a tau||+||b tau|| >= 11/13}.             (1)
```

In particular every point of `G_U` also satisfies the two individual
anti-tooth inequalities

```text
||a tau|| >= 9/26,       ||b tau|| >= 9/26.                (2)
```

### Proof

Put `epsilon=2/13`.  Suppose first that the two odd runners are eligible and
their unique nearest integers `N_x,N_y` have opposite parity.  Then

```text
M = (N_x+N_y)/2,       N = (N_x-N_y)/2
```

are half-integers.  With

```text
e_a = a tau-M,       e_b = b tau-N,
```

the signed errors of the two odd runners are `e_a+e_b` and `e_a-e_b`.
The elementary identity

```text
max(|e_a+e_b|,|e_a-e_b|) = |e_a|+|e_b|                   (3)
```

therefore says that both odd errors have magnitude at most `epsilon` exactly
when `|e_a|+|e_b|<=epsilon`.  Since `epsilon<1/2`, `M,N` are the unique nearest
half-integers to `a tau,b tau`.  If

```text
h_m(tau)=dist(m tau,Z+1/2)=1/2-||m tau||,
```

this is `h_a+h_b<=2/13`, equivalently (1).

Conversely, choose the unique nearest half-integers to `a tau,b tau` under
the inequality in (1).  Their sum and difference are integers of opposite
parity, (3) makes them the nearest integers to `x tau,y tau`, and both errors
are at most `2/13`.  This recovers exactly THM-769's closed eligibility and
opposite-colour predicate, including its endpoints.  The tightness criterion
now follows from THM-769.  Finally each norm in (1) is at most `1/2`, so the
other is at least `11/13-1/2=9/26`, proving (2).

For completeness, (1a) is a global version of the same calculation.  Choose
signed representatives `r,s` of `x tau/2,y tau/2`, and put
`u=min(|r|,|s|)`, `v=max(|r|,|s|)`.  The two lift-clearance minima are `u`
and `1/2-v`, while

```text
||r+s||+||r-s|| = min(2v,1-2u).
```

Therefore

```text
1-2B = min(1-2u,2v) = Q,
```

as claimed.  This signed-representative identity is essential: replacing it
by an unsigned phase sum before taking circle norms is not valid.

## 2. Exact measure formula

Let `mu` be normalized Lebesgue measure on `R/Z`.  Put

```text
g = gcd(a,b) = gcd(x,y),       alpha=a/g,       beta=b/g,
r = alpha-beta = y/g,          s=alpha+beta = x/g,
n = floor(2 alpha/13 + 1/2).
```

Then `alpha>beta` are coprime and of opposite parity, and

```text
mu(H_(a,b))
  = 2/(r s) sum_(j=1)^n min(4r/13, 4 alpha/13-(2j-1)).     (4)
```

### Proof

Multiplication by `g` preserves normalized measure, so it is enough to use
the reduced pair `alpha,beta`.  The set

```text
h_alpha(tau)+h_beta(tau) <= 2/13
```

lies in the disjoint `2/13`-neighbourhoods of the `alpha` half-grid points

```text
tau_k=(2k+1)/(2 alpha).
```

At these centres the positive values of `2 alpha h_beta(tau_k)` below
`4 alpha/13` are precisely

```text
q=1,3,...,2n-1,
```

each occurring twice.  This is the permutation of the odd half-grid residues
induced by coprimality and opposite parity.

Fix such a value `q`.  On the side moving toward the adjacent `beta`
half-grid point, `h_beta` first falls with slope `beta` and then rises; on the
other side it rises immediately.  Adding the two affine pieces shows that
the component at either of the two corresponding centres has length

```text
1/(r s) min(4r/13, 4 alpha/13-q).                          (5)
```

The first term in (5) is the branch that crosses the `beta` half-grid point;
the second is the branch that stops before it.  Summing (5) over the two
centres for every qualifying odd `q` gives (4).

## 3. The sharp universal cap

Uniformly over all distinct positive odd `x,y`,

```text
mu(H_(a,b)) <= 8/117.                                     (6)
```

The bound is sharp: equality is attained by the reduced pair

```text
(alpha,beta)=(5,4),       equivalently x:y=9:1,
```

and hence by its common odd dilates.  Therefore every tight two-sheet packet
obeys the new necessary condition

```text
measure(G_U) <= 8/117.
```

### Proof

The summands in (4) give the two bounds

```text
mu(H_(a,b)) <= min(
    8n/(13s),
    2n(4 alpha/13-n)/(r s)
).                                                         (7)
```

If `s>=9n`, the first term in (7) is at most `8/117`.  Suppose `s<9n`.
Because `s>alpha` and `u(2 alpha-u)` decreases for `u>alpha`,

```text
4rs > 36n(2 alpha-9n).                                    (8)
```

When `4 alpha>=23n`, the right side of (8) is at least

```text
117n(4 alpha/13-n),
```

so the second term in (7) is at most `8/117`.

It remains only to treat `4 alpha<23n`.  Since

```text
n <= 2 alpha/13+1/2,
```

this inequality is impossible for `alpha>=25`.  Direct substitution for
`alpha<25` leaves exactly the five rows below; maximizing (4) over the
admissible coprime opposite-parity `beta` gives

| `alpha` | `n` | maximum measure | maximizing `beta` |
|---:|---:|---:|---:|
| 4 | 1 | `6/91` | 3 |
| 5 | 1 | `8/117` | 4 |
| 10 | 2 | `8/169` | 3 |
| 11 | 2 | `16/273` | 10 |
| 17 | 3 | `2/39` | 16 |

This proves (6) and sharpness.

## 4. Scope and audit

The cap is a scalar necessary condition, not a closure theorem.  For example,
the ten-speed cores obtained from `{1,...,11}` by deleting `2`, `4`, or `10`
have respectively

```text
measure(G_U) = 569/10010, 499/10010, 746/15015,
```

all below `8/117`.  The proof-facing predicate remains the pointwise
containment (1), with the positions and widths of all components retained.
An ordinary tournament ranking candidate diamonds by size destroys precisely
that containment data.

Divisor transfer does not repair this scalar loss.  Among the `3,400`
divisor-complete ten-cores in `[1,19]`, exactly `52` already satisfy
`measure(G_U)<=8/117`; the minimum is

```text
measure(G_U)=41/858
at U={1,2,3,5,7,8,9,10,11,12}.
```

Section 5 nevertheless rules out every one of these cores pointwise.  Thus
the divisor pins plus the sharp measure cap still do not imply noncoverage;
component incidence is genuinely stronger.

The companion exact audit is

```text
04-computation/lrc13_folded_diamond_two_sheet_obstruction_codex_S3.py
05-knowledge/results/lrc13_folded_diamond_two_sheet_obstruction_codex_S3.out
```

It checks the predicate identity on every affine atom and endpoint for `1,600`
opposite-parity pairs through `a=80`, matches (4) to an independent
piecewise-linear integration for `5,217` reduced pairs through `alpha=160`,
verifies the complete exceptional table, and records Tournament Analysis with
the loss of pointwise containment stated explicitly.

## 5. Finite-exact closure for every quotient core in `[1,19]`

There is no tight two-sheet packet

```text
A=2U union {x,y}
```

with `U subset {1,...,19}` and `|U|=10`, even when the distinct positive odd
speeds `x,y` are otherwise unbounded.

### Exact reduction and certificate

Let `ell_max(U)` be the length of the widest connected component of `G_U`.
If one odd speed `w` is eligible throughout `G_U`, every component must lie
in one connected closed `w`-tooth, whose length is `4/(13w)`.  Therefore

```text
w <= floor(4/(13 ell_max(U))).                             (9)
```

This is an intrinsic finite cap, not a sampled speed cutoff.  The companion
certificate

```text
04-computation/lrc13_s2_unbounded_odd_core19_closure_codex_S3.py
05-knowledge/results/lrc13_s2_unbounded_odd_core19_closure_codex_S3.out
```

constructs every strict component with rational endpoints for all

```text
binomial(19,10)=92,378
```

cores, then checks every positive odd `w` permitted by (9), accepting equality
in the closed eligibility tooth.  Across `767,700` exact `(U,w)` incidences it
finds zero runner eligible throughout `G_U`.  Hence it is already impossible
for two runners to satisfy the stronger opposite-colour condition.  The
largest intrinsic cap is `60`, for

```text
U={1,2,3,5,7,8,11,13,18,19},
ell_max(U)=6/1183,
```

whose loose set has eighteen components.  The canonical endpoint/decision-row
digest is

```text
ec206bf06eda11b5f8ee5318b2bdbc97d61ae63c78e508c83610ce3a8a2dcf83.
```

No primitivity or divisor-completeness filter is used; only `3,400` of the
cores happen to satisfy THM-772's divisor pins.  Thus this is a complete
low-core slice, but it is not a uniform bound on `max(U)`.

## 6. Exact max-peel recursion

The sporadic maximum-deletion condition has two sharply different forms in a
two-sheet packet.

If `max(A)` is one of the odd exceptions, deleting it leaves

```text
P=2U union {z}
```

for the other odd exception `z`.  At a maximizer of `U`, the two lifts give
`z`-clearances `c` and `1/2-c`, so the better is at least `1/4`.  Therefore

```text
M(P)>=min(M(U),1/4)>=1/11>1/12.                           (10)
```

The last inequality uses the settled lower-dimensional LRC bound for the
ten-speed core.  Thus the sporadic condition is automatic in the odd-maximum
branch and supplies no extra rigidity.

Suppose instead that `max(A)=2R`, where `R=max(U)`, and put

```text
U^-=U\{R},
E_L(U^-)={tau:phi_(U^-)(tau)>L},
D_R={tau:||R tau||<=1/13},
Q(tau)=||a tau||+||b tau||.
```

Then tightness of `A` and sporadicity of its maximum deletion imply exactly

```text
E_(1/13)(U^-) subset D_R union {Q>=11/13},                (11)
E_(1/12)(U^-) intersect D_R intersect {Q<5/6} nonempty.   (12)
```

### Proof

For every `tau`, the two lifts under doubling have the same even-core
clearance.  Applying (1a) gives the exact max formulas

```text
M(A)=max_tau min(phi_(U^-)(tau),||R tau||,B_(x,y)(tau)),
M(A\{2R})=max_tau min(phi_(U^-)(tau),B_(x,y)(tau)).        (13)
```

Since `B=(1-Q)/2`, the global upper bound `M(A)<=1/13` is precisely (11).
The strict inequality `M(A\{2R})>1/12` is equivalent to the existence of a
point in

```text
E_(1/12)(U^-) intersect {Q<5/6}.
```

Because `5/6<11/13`, this point is outside the diamond term in (11), so it
must belong to `D_R`.  This proves (12).

Equations (11)--(12) are a quantifier-exact nine-core collar target, not a
claim that the collar is empty.  Proving its impossibility, or showing that it
forces THM-775's dyadic seam and a well-founded descent, would close the
sporadic `s=2` max-peel residue.
