---
id: THM-780
title: The quantitative strict-margin-to-measure bridge — phase pigeonhole gives mu_alpha >= ceil(1/(beta-alpha))^(-d), with Haar compactness as the relation-lattice limit proof
status: PROVED (elementary phase partition; independent Haar/character proof included; both proofs independently audited)
source: codex-2026-07-14-S10
depends_on:
  - LRC(<=13)
related: [THM-777, HYP-4452, HYP-4472, MISTAKE-117]
---

# THM-780 — Quantitative strict-margin-to-measure bridge

## Statement

Write `T=R/Z`, let `||x||` denote distance to `0` in `T`, and for
`p=(p_1,...,p_d) in Z^d` put

```
mu_alpha(p) = Leb{t in T : ||p_i t|| >= alpha for every i}.
```

> **Strict-margin measure theorem.**  Fix an integer `d>=1`, a class
> `A subset Z^d`, and `0<=alpha<beta<=1/2`.  If every `p in A` has a
> `beta`-deep time,
> `max_t min_i ||p_i t|| >= beta`, then, with
>
> `M=ceil(1/(beta-alpha))`,
>
> every `p in A` satisfies the explicit bound
>
> `mu_alpha(p) >= M^(-d)`.

No height bound, primitivity assumption, divisor hypothesis, or classification
of limiting subtori is required.

## First proof: finite phase pigeonhole

Fix `p in A` and choose a `beta`-deep time `t_0`.  Partition each coordinate
circle into the `M` half-open intervals

```
[j/M,(j+1)/M),       j=0,...,M-1,
```

and hence partition `T^d` into `M^d` half-open cubes.  For the orbit map

```
phi:T -> T^d,       phi(s)=s p,
```

the `M^d` measurable preimages partition `T`.  One cube `Q` therefore has

```
Leb(phi^(-1)(Q)) >= M^(-d).                               (1)
```

Put `A_Q=phi^(-1)(Q)` and choose `s_0 in A_Q`.  If `s in A_Q`, the two points
`phi(s)` and `phi(s_0)` lie in the same coordinate intervals.  Thus, for
`u=s-s_0`,

```
||p_i u|| < 1/M <= beta-alpha       for every i.          (2)
```

Translation preserves Lebesgue measure, so the simultaneous Bohr-return set

```
B={s-s_0:s in A_Q}
```

has measure at least `M^(-d)`.  For every `u in B`, the circle triangle
inequality and (2) give

```
||p_i(t_0+u)|| >= ||p_i t_0||-||p_i u|| > alpha.          (3)
```

Hence `t_0+B` is contained even in the strict `alpha`-safe set.  Equation (1)
proves `mu_alpha(p)>=M^(-d)`.  ∎

This proof exposes the missing finite object particularly cleanly: not a
height cutoff or a witness denominator, but one heavy cell in the joint phase
orbit.  Subtracting an anchor in that cell turns phase concentration into a
simultaneous return; the already-known deeper witness turns the return into a
positive safe interval set.

## Bohr-packet strengthening

The proof retains more than its scalar measure conclusion. For a chosen
`beta`-deep time `t_0`, define the simultaneous return neighborhood

```text
B_p(delta)={u in T : ||p_i u||<delta for every i},
delta=beta-alpha.
```

Then

```text
t_0+B_p(delta) subset {t : ||p_i t||>alpha for every i},     (BP1)
Leb(B_p(delta))>=M^(-d).                                    (BP2)
```

Indeed, (BP1) is the circle triangle inequality. For (BP2), the translated heavy
cell difference set `B` constructed above has measure at least `M^(-d)` and is
contained in `B_p(delta)` because `||p_i u||<1/M<=delta` on `B`.

Thus a theorem that needs incidence rather than mass may preserve the packet
`(t_0,B_p(delta))` instead of compressing immediately to `mu_alpha(p)`. A
coarser structural proxy is a heavy phase-cell address plus the
character-relation lattice, but only if the parametrized orbit map (equivalently
its degree/gcd and residue action) is also retained: `p` and `gp` have the same
annihilator subgroup but different pullback Bohr sets and event incidence. A
concrete next interface is intersection of the exact packet with an
owner-labelled endpoint/event coset; that is an inhomogeneous
residue/character question which keeps owner identity while discarding raw
height.

**Exact event-coset interface.** Let `u` be a positive integer event frequency,
`E(u,gamma)={t in T:ut=gamma}`, and `eta=gamma-u t_0`. Define

```text
psi_(p,u)(v)=(p_1v,...,p_dv,uv),       H=im psi_(p,u).
```

Then

```text
E(u,gamma) intersect (t_0+B_p(delta)) != empty
iff H intersects (-delta,delta)^d x {eta}.                 (BP3)
```

Equivalently, after choosing a real lift of `eta`, test the `u` candidates
`v_k=(eta+k)/u`, `k in Z/uZ`; (BP3) asks whether the affine cyclic code
`k |-> (p_i(eta+k) mod u)_i` enters the centered box of radius `u*delta`.
For Boolean existence the augmented annihilator
`ker((m,n)|->m dot p+n u)` and target phase describe the subgroup slice. For
component count or transport the orbit covering degree, gcd multiplicities,
and owner label must also survive. This is an exact recursive interface, not a
claim that every such slice is nonempty.

## Second proof: Haar compactness and the limiting relation lattice

The following longer proof gives the structural limit statement that motivated
the theorem.  It proves positivity without the explicit constant and audits
what information survives an unbounded-height sequence.

Suppose instead that `p^(n) in A` and `mu_alpha(p^(n))->0`.  Let

```
phi_n : T -> T^d,       phi_n(t)=t p^(n),
H_n = phi_n(T),         nu_n=(phi_n)_* Leb.
```

Thus `nu_n` is normalized Haar measure on the closed one-parameter subgroup
`H_n` (a circle subgroup or a point). Define

```
C_alpha = {x in T^d : ||x_i|| >= alpha for all i},
O_alpha = {x in T^d : ||x_i|| >  alpha for all i},
K_beta  = {x in T^d : ||x_i|| >= beta  for all i}.
```

Then `mu_alpha(p^(n))=nu_n(C_alpha)`.  By the hypothesis choose
`z_n in H_n intersect K_beta`.  After taking a subsequence, compactness gives
`z_n -> z in K_beta`.

Enumerate the countable character group `Z^d`.  A diagonal subsequence makes,
for every `m in Z^d`, the truth value of

```
<m,p^(n)> = 0
```

eventually constant.  Let `L` be the set of characters for which it is
eventually true.  This is a subgroup of `Z^d`; it is saturated, because
`k m in L` for a nonzero integer `k` implies eventually
`k<m,p^(n)>=0`, hence eventually `<m,p^(n)>=0`.

Put

```
H = L^perp = {x in T^d : <m,x>=0 in T for every m in L}
```

and let `nu_H` be normalized Haar measure on `H`.  For every character `m`,

```
hat(nu_n)(m) = 1_{<m,p^(n)>=0} -> 1_{m in L} = hat(nu_H)(m).
```

The final equality uses saturation of `L`, so that the annihilator of `H` is
exactly `L`.  Trigonometric polynomials are uniformly dense in `C(T^d)`;
therefore this pointwise Fourier convergence implies `nu_n => nu_H` weakly.

The deep witness survives the limit.  Indeed, if `m in L`, then for all large
`n`,

```
<m,z_n> = t_n <m,p^(n)> = 0 in T,
```

so `<m,z>=0`.  Hence `z in H`.  Since `z in K_beta` and `beta>alpha`,
`z in O_alpha`.  Consequently `H intersect O_alpha` is a nonempty relatively
open subset of `H`.  Haar measure has full support on a compact group, so

```
nu_H(O_alpha)>0.
```

On the other hand, Portmanteau applied to the open set `O_alpha` gives

```
0 < nu_H(O_alpha)
    <= liminf_n nu_n(O_alpha)
    <= liminf_n nu_n(C_alpha)
     = 0,
```

a contradiction.  This proves the theorem.  (For completeness, full support is
elementary: if a nonempty open identity neighborhood had Haar mass zero,
finitely many translates covering the compact group would also have total mass
zero.)  ∎

## The twelve-core consequence

Settled `LRC(<=13)` says that every 12-speed core has a `1/13`-deep time.  Take
`d=12`, `beta=1/13`, and `alpha=1/14`.  Since

```
1/13-1/14=1/182,
```

the first proof has `M=182` and gives every 12-core `P` the explicit bound

```
|G'_P| = Leb{t : min_{p in P} ||pt|| >= 1/14}
       >= 182^(-12)
        = 1/1320859596446125189798629376.                 (4)
```

Together with THM-777(1), this gives the genuinely uniform effective bound

```
rho(P) <= 12*182^12/pi < 5.046*10^27.                     (5)
```

Thus the qualitative asymptotic floor and boundedness of the normalized
`rho`-coordinate are proved with an explicit (deliberately crude) constant.
This does not by itself produce a finite or terminal atlas: endpoint-owner,
residue, and deck-transport fibres may still vary inside that scalar domain.
What remains open in THM-777 is the sharp identity `inf |G'_P|=7/858`, its
uniqueness statement, and a practically sized structural cutoff.

## What changed relative to the older compactness route

HYP-4472 tried to prove continuity of the safe functional on a hand-built
one-/two-torus compactification and then classify its zero locus.  Neither step
is needed.  Quantitatively, the heavy phase cell and its anchored difference
set already give a Bohr neighborhood of mass `M^(-d)`.  Structurally, the
correct limit state is the **entire stabilized character relation lattice**,
whose annihilator automatically supplies the limiting compact subgroup of any
rank.  Only the open-set half of Portmanteau is needed, not continuity of the
safe functional.  In both proofs the decisive datum is the strictly deeper
`beta`-witness: it absorbs a simultaneous return of radius `beta-alpha`.  This
is why the theorem applies at `1/14` using settled `1/13`, but does not prove the
older `2/25` floor.

## Assumption challenge and Tournament Analysis

The useful vertices here are not runners.  In the finite proof they are joint
phase cells; in the limit proof they are character obligations `m in Z^d`.  A
runner-level pair tournament forgets exactly the simultaneous datum used by
both: common-cell membership, or which integer characters vanish together.
The finite switch is “two orbit points occupy the same cell”; the limit switch
is `m·p^(n)=0`; their quotient objects are respectively the anchored Bohr-return
set and the annihilator `L^perp`.  These preserve the safe predicate while
discarding height.

A forced tournament adds no theorem-facing information here: same-cell and
relation switches are equivalence/subgroup closure, not antisymmetric
dominance.  The exact fingerprints are a phase partition and a saturated
lattice rather than scores, cycles, or Hamiltonian paths.  This is a concrete
case where challenging “vertices must be runners” finds the missing state and
also explains why Tournament Analysis should record an honest no-clean-switch
result instead of manufacturing an orientation.

## Audit

An independent proof audit checked both routes.  For the finite proof it
verified the half-open boundary and circle-wrap cases, translation of the
preimage mass, and the strict triangle-inequality margin.  For the limit proof
it checked diagonal stabilization, saturation and double annihilation, Fourier
determination of weak Haar convergence, survival of the deep witness, and the
direction of the Portmanteau inequality.  No admissibility or boundary gap was
found.
