---
id: THM-2379
title: "Anchored guard/unit deletion factor-repair current"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED. An anchored
  ordinary-unit or guard failure which retains the deepest safe factor
  has a two-coordinate danger-probe corner of fixed negative sign.
  Replacing the translated failure probe by its safe complement
  reverses every nonzero failure colour and forces a positive
  off-diagonal deep-by-repair coefficient at least 11rho/24336 for an
  ordinary unit and rho/2704 for the guard, together with exact
  squared-energy floors. Both displayed probe multipliers are coprime
  to 91. If a THM-2370 drift-selected physical deletion layer has this
  typing, then rho>=delta/n and the same bounds hold with rho replaced
  by delta/n. Retaining one of four labelled nondeep-blocker status
  words costs a factor four in mass/coefficient and sixteen in squared
  energy. Exact open chambers can have every blocker safe, so this is a
  factor-repair current, not a canonical owner or lawful target
  current. No scalar-row exclusion, ledger decrement, or LRC(14)
  consequence follows.
source: codex-2026-07-25-guard-unit-factor-repair
depends_on:
  - THM-2370-deletion-martingale-drift-conservation-and-sharp-clone-hostile
related:
  - THM-2269-marked-expiration-root-spectrum-and-branch-state-no-go
  - THM-2362-thirteen-shift-successor-statistic-and-role-jet-floor
  - THM-2364-anchored-corner-forces-mixed-deep-blocker-colour
script: 04-computation/lrc14_guard_unit_factor_repair_thm2379.py
output: 05-knowledge/results/lrc14_guard_unit_factor_repair_thm2379.out
script_sha256: 7fa8f2ad0053e84a77a0286e75846fc3ca3d3e0d6adcc68c8b2a9d4194ea4e5b
output_sha256: b9c50b4b68b38d4c2d545b1f0a53b06b79095249663bc88920f5ca43fcdae3da
hash_basis: working-tree bytes (LF)
---

# THM-2379 -- anchored guard/unit deletion factor repair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY HOSTILE-AUDITED.**

THM-2370 shows that drift lost while safe factors are inserted must
enter a physical deletion layer, but a deletion of the guard or an
ordinary unit lies outside the canonical blocker core. Such a layer is
not colourless. Its failed factor has an exact complementary
off-diagonal repair current:

```text
anchored guard/unit danger
  + retained deepest safe factor
  -> negative deep-by-danger corner
  -> positive deep-by-complement repair colour.                (1)
```

The repair colour belongs to a translated duplicate factor. It is not
yet a target coordinate or a canonical owner label. Section 7 gives
exact open chambers showing that this distinction is necessary.

## 1. Two anchored roles

Work on the circle `T=R/Z`. Put

```text
p=13,

zeta=exp(2*pi*i/p),

d_L(y)=1_(||y||<L/14),

u_L=1-d_L,                         L in {1,2}.        (2)
```

Here `L=1` is an ordinary danger factor and `L=2` is the guard-danger
factor. All support statements and identities below are almost
everywhere; strict-open endpoints form a finite null set.

Let `c,v` be positive integers and let `alpha,beta in T`. Let `w>=0`
be integrable, with

```text
rho=int_T w(x)dx>0,

support(w) subset {x:d_1(cx-alpha)=0}
                    intersection
                  {x:d_L(vx-beta)=1}                (3)
```

almost everywhere. Thus `w` retains the deepest safe role at `c` and
records a failure of the safe role at `v`.

For `r,s in F_p`, define

```text
K^-_(r,s)
 =int_T w(x)
    d_1(cx-alpha-r/p)
    d_L(vx-beta-s/p) dx,

K^+_(r,s)
 =int_T w(x)
    d_1(cx-alpha-r/p)
    u_L(vx-beta-s/p) dx.                            (4)
```

For `epsilon in {-,+}`, the normalized finite transforms are

```text
Khat^epsilon(a,b)
 =1/p^2 sum_(r,s in F_p)
    K^epsilon(r,s)zeta^(a r+b s).                   (5)
```

The signs in (5) are fixed throughout. Reversing both signs merely
conjugates every displayed coefficient.

The two anchors give

```text
K^epsilon(0,s)=0                    for every s,

K^+(r,0)=0                          for every r.     (6)
```

The second line is the promised off-diagonal repair: the safe
complement cannot coexist with the original failure at shift zero.

## 2. Exact thirteen-shift counts

For `L=1,2`, exact interlacing of the strict interval endpoints with the
thirteen translated roots gives

```text
sum_(s in F_p)d_L(y-s/p)
 =2L-d_L(py)                                      (7)
```

almost everywhere. Equivalently,

```text
sum_s u_L(y-s/p)
 =p-2L+d_L(py).                                    (8)
```

For later use, finite character orthogonality says

```text
sum_(a!=0)zeta^(ar)
 =p-1,       r=0,

 =-1,        r!=0.                                 (9)
```

Equations (7)--(9), rather than an asymptotic Fourier estimate, drive
the result.

## 3. The signed mixed corner

Sum (5) over all `a,b!=0`. Interchanging finite sums and the integral,
then using (9), factors the result into one root sum for each probe.
On the support (3), the deep bracket is

```text
p d_1(cx-alpha)
 -sum_r d_1(cx-alpha-r/p)

=-(2-d_1(p(cx-alpha))),                            (10)
```

while the failure bracket is

```text
p d_L(vx-beta)
 -sum_s d_L(vx-beta-s/p)

=p-2L+d_L(p(vx-beta)).                             (11)
```

Therefore

```text
S^-
 :=sum_(a!=0,b!=0)Khat^-(a,b)

 =-1/p^2 int_T w(x)
      (2-d_1(p(cx-alpha)))
      (p-2L+d_L(p(vx-beta))) dx.                  (12)
```

Every factor in the last integrand is real and nonnegative, with

```text
2-d_1(p(cx-alpha))>=1,

p-2L+d_L(p(vx-beta))>=p-2L.
```

Hence

```text
S^-<=-(p-2L)rho/p^2.                               (13)
```

This is an exact sign statement, not merely a nonzero energy.

## 4. Complement transfer and quantitative floors

For every `r,s`,

```text
K^-_(r,s)+K^+_(r,s)
 =int_T w(x)d_1(cx-alpha-r/p)dx,
```

which is independent of `s`. Consequently,

```text
Khat^+(a,b)=-Khat^-(a,b)             whenever b!=0. (14)
```

Combining (12)--(14),

```text
S^+
 :=sum_(a!=0,b!=0)Khat^+(a,b)
 >=(p-2L)rho/p^2.                                  (15)
```

There are `(p-1)^2=144` mixed colours. Thus some `a,b!=0` obey

```text
Re Khat^+(a,b)
 >=(p-2L)rho/(p^2(p-1)^2).                         (16)
```

At `p=13`, this is

```text
L=1, ordinary unit:
  Re Khat^+(a,b)>=11rho/24336;

L=2, guard:
  Re Khat^+(a,b)>=9rho/24336=rho/2704.             (17)
```

Cauchy--Schwarz applied to either signed corner also gives

```text
sum_(a!=0,b!=0)|Khat^epsilon(a,b)|^2
 >=(p-2L)^2 rho^2/(p^4(p-1)^2).                    (18)
```

The coefficient floor in (17) is stronger than an unsigned energy
existence statement: the complement repair has a chosen positive real
direction in the common cyclotomic embedding.

## 5. Exact probe multipliers

Assume now that `w` is a rational-interval step weight, as in the LRC
packets. Poisson-smooth the two displayed probes at a common radius
before expanding. At every finite radius the double series is
absolutely convergent, and (5) selects probe multipliers

```text
m=a mod 13,

n=b mod 13.                                        (19)
```

For `k!=0`,

```text
dhat_L(k)=sin(pi L k/7)/(pi k),

uhat_L(k)=-dhat_L(k).                               (20)
```

Since `L in {1,2}`, a nonzero factor in (20) has `7` not dividing `k`.
Equation (19) has `13` not dividing `m n`. Thus every nonzero summand
displayed by the two probes satisfies

```text
gcd(m,91)=gcd(n,91)=1.                              (21)
```

The boundary coefficient is their joint Poisson--Abel limit. Statement
(21) is deliberately limited to the deep and repair probe multipliers.
Fourier indices inside the collapsed ancestry weight `w` are not
controlled by this theorem, and `(a,b)` is not a relation target.

## 6. A typed THM-2370 deletion consequence

Use THM-2370's notation. Suppose a drift-selected layer `L_j` satisfies

```text
D(L_j)>=delta^2/N^2,                               (22)
```

where `N` is the number of inserted masks. Assume additionally that:

1. `L_j` is a nonnegative physical Boolean deletion layer;
2. it deletes an ordinary-unit-safe factor (`L=1`) or the guard-safe
   factor (`L=2`); and
3. the canonical deepest safe factor is still retained.

These are extra typing hypotheses; (22) alone does not imply them.
With normalized counting norm,

```text
delta/N
 <=||Q L_j||_2
 <=||L_j||_2
 <=||L_j||_infinity.                               (23)
```

Choose a maximizing cell `(r,s,t)`. A positive cell has `r!=t`, because
the retained deepest safe factor at shift `t` is disjoint from the
displayed THM-2365 danger probe `Delta_r`. Gauge the deepest safe shift
to `alpha=t/p`, remove only `Delta_r` from the cell integrand, and call
the remaining nonnegative integrand `w`. The deleted safe factor still
appears there as its danger complement at its selected shift; absorb
that shift into `beta`. Removing an indicator can only increase mass,
so

```text
rho=int_T w
 >=L_j(maximizing cell)
 >=delta/N.                                        (24)
```

The deleted safe factor appears in `w` as its danger complement, and
the retained deepest factor gives exactly (3). Applying (17) yields

```text
ordinary-unit deletion:
  Re Khat^+(a,b)>=11delta/(24336 N);

guard deletion:
  Re Khat^+(a,b)>=delta/(2704 N)                    (25)
```

for some `a,b!=0`, with the probe typing (21).

No assertion is made for an abstract Hilbert deletion, for a layer
deleting a blocker factor, or after the deepest safe factor has itself
been removed.

## 6a. A blocker-status word costs four in amplitude, sixteen in energy

At the maximizing cell, partition `w` by the danger/safe bits of the
two nondeep blockers. This is an exact disjoint four-way partition

```text
w=sum_(sigma subset {1,2}) w_sigma,

rho=sum_sigma rho_sigma.                            (26)
```

Each `w_sigma` remains nonnegative, rational-interval typed, deepest
safe, and anchored in the same guard/unit failure. Some status word has

```text
rho_sigma>=rho/4.                                   (27)
```

Applying the abstract theorem to that one word gives

```text
ordinary-unit deletion:
  Re Khat^+_sigma(a,b)>=11delta/(97344 N);

guard deletion:
  Re Khat^+_sigma(a,b)>=delta/(10816 N).            (28)
```

The corresponding squared-energy floors are

```text
ordinary-unit deletion:
  sum_(a,b!=0)|Khat^+_sigma(a,b)|^2
    >=121 delta^2/(65804544 N^2);

guard deletion:
  sum_(a,b!=0)|Khat^+_sigma(a,b)|^2
    >=9 delta^2/(7311616 N^2).                      (28a)
```

A nonempty `sigma` retains literal blocker-danger labels. It is still
not a canonical THM-2305 word, because the guard/unit failure remains
present. The empty status is genuinely possible.

## 7. Exact empty-blocker chambers

The stronger local claim

```text
guard/unit failure -> some blocker danger
```

is false even on open rational chambers. Take

```text
H=1,

(q_1,...,q_5)=(2,339,677,1015,1353),

(c_1,c_2,c_3)=(13,169,2197).                       (29)
```

At

```text
x=1/8
```

only the guard is dangerous; all five ordinary units and all three
blockers are safe. The same status persists for

```text
|x-1/8|<1/25256.                                   (30)
```

The tight margin in (30) is the `q_5=1353` distance from `1/8` to the
ordinary danger boundary.

At

```text
x=1/2
```

only `q_1=2` is dangerous; the guard, the other four ordinary units,
and all three blockers are safe. This status persists for

```text
|x-1/2|<3/15379.                                   (31)
```

Here the tight displayed margin comes from the largest odd blocker
`c_3=2197`. These are local factor-status hostiles, not LRC
counterexamples. They prove that scalar-cover information on the
all-safe core cannot by itself promote (25) to a blocker owner.

## 8. Exact boundary and next service

The theorem supplies:

```text
typed guard/unit deletion drift
  -> positive failure mass
  -> positive deep-by-complement repair colour
  -> two displayed 91-unit probe multipliers.       (32)
```

It preserves the failure ancestry, the chosen deletion clock, and, after
(26), a complete nondeep-blocker status word. It loses the lawful
two-dimensional target quotient, the endpoint phase, and canonical
owner/current service.

In particular, `K^+` is a translated complement-overlap probe. Because
`K^+(r,0)=0`, it does **not** restore the deleted safe factor at its
original shift. The cheapest next bridge is an identification of the
repair shift with one lawful target dipole while keeping the blocker
status and same-clock endpoint. Without that identification, THM-2369's
charged target edge is still unavailable.

No scalar row is excluded. The ledger remains `165`, and LRC(14)
remains open.

## 9. Exact companion

Run

```bash
python 04-computation/lrc14_guard_unit_factor_repair_thm2379.py
python -O 04-computation/lrc14_guard_unit_factor_repair_thm2379.py
```

and compare both transcripts, after LF normalization, with
`05-knowledge/results/lrc14_guard_unit_factor_repair_thm2379.out`.
The companion verifies:

- both `L=1,2` thirteen-shift count laws on the full rational cell
  refinement: `26+26` disjoint boundaries and `52` cells;
- all `315` ordinary-unit and `585` guard anchored refinement-cell
  pairs, representing respectively `69` and `161` distinct pairs of
  thirteen-bit profiles, including the signed corner factors and
  constants in (12)--(18);
- a symmetric-lift-index `324`-multiplier representative bank with
  `277` live lifts, `47`
  septimal zeros, and all `76,729` live ordered pairs satisfying the
  `91`-unit probe gate in (19)--(21);
- the unsplit and four-way THM-2370 constants; and
- every factor status and open radius in (29)--(31), together with the
  two THM-2269 control chambers.

Ordinary and optimized runs are byte-identical to the stored
transcript. The LF SHA-256 hashes are recorded in the frontmatter.

An independent hostile audit reconstructed the shift count, both
corner signs, all coefficient and energy denominators, the
THM-2370 norm-to-cell extraction, the four-way loss, the `91`-unit
scope, and both open hostile chambers. It caught four reporting
defects: the first line of (6) had a sign-placeholder typo; squared
energy pays sixteen rather than four under the status split; the
`315/585` counts are refinement-cell pairs rather than distinct
bit profiles; and the multiplier bank is symmetric in lift index, not
as an integer set. All four repairs are incorporated above and in the
exact companion. QED.
