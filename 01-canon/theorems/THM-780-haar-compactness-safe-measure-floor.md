---
id: THM-780
title: Haar-compactness turns every strict lonely-runner margin into a uniform safe-measure floor
status: PROVED (elementary compactness, character stabilization, and Portmanteau; independent audit pending)
source: codex-2026-07-14-S10
depends_on:
  - LRC(<=13)
related: [THM-777, HYP-4452, HYP-4472, MISTAKE-117]
---

# THM-780 — Haar-compactness safe-measure floor

## Statement

Write `T=R/Z`, let `||x||` denote distance to `0` in `T`, and for
`p=(p_1,...,p_d) in Z^d` put

```
mu_alpha(p) = Leb{t in T : ||p_i t|| >= alpha for every i}.
```

> **Strict-margin compactness theorem.**  Fix an integer `d>=1`, a class
> `A subset Z^d`, and `0<=alpha<beta<=1/2`.  If every `p in A` has a
> `beta`-deep time,
> `max_t min_i ||p_i t|| >= beta`, then
> `inf_{p in A} mu_alpha(p)>0`.

No height bound, primitivity assumption, divisor hypothesis, or classification
of limiting subtori is required.  The constant supplied by this proof is
non-effective.

## Proof

Suppose instead that `p^(n) in A` and `mu_alpha(p^(n))->0`.  Let

```
phi_n : T -> T^d,       phi_n(t)=t p^(n),
H_n = phi_n(T),         nu_n=(phi_n)_* Leb.
```

Thus `nu_n` is normalized Haar measure on the cyclic subgroup `H_n`.  Define

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
`d=12`, `beta=1/13`, and `alpha=1/14`.  There is therefore an absolute,
possibly non-explicit constant `c_12>0` such that every 12-core `P` satisfies

```
|G'_P| = Leb{t : min_{p in P} ||pt|| >= 1/14} >= c_12.
```

Together with THM-777(1), this gives the genuinely uniform but non-effective
bound

```
rho(P) <= 12/(pi c_12).
```

Thus the qualitative asymptotic floor and the existence of a bounded normalized
regime-2 atlas are proved.  What remains open in THM-777 is the sharp effective
identity `c_12=7/858`, its uniqueness statement, and any usable numerical atlas
cutoff.

## What changed relative to the older compactness route

HYP-4472 tried to prove continuity of the safe functional on a hand-built
one-/two-torus compactification and then classify its zero locus.  Neither step
is needed here.  The correct state object is the **entire stabilized character
relation lattice**, whose annihilator automatically supplies the limiting
compact subgroup of any rank.  The correct analytic input is only the
lower-semicontinuous open-set half of Portmanteau.  Most importantly, the
strictly deeper `beta`-witness is transported into that subgroup, where it gives
an interior safe point.  This is why the proof applies at `1/14` using the
settled `1/13` theorem, but does not prove the older `2/25` floor.

## Assumption challenge and Tournament Analysis

The useful vertices here are not runners.  A runner-level pair tournament
forgets exactly the data used by the proof: which integer characters vanish
simultaneously and which limit subgroup they cut out.  The alternative vertex
set is the countable family of proof obligations/characters `m in Z^d`, with
the binary stabilized switch `m·p^(n)=0` and the annihilator `L^perp` as the
quotient object.  This preserves all torus equations and the safe predicate in
the limit; it destroys finite height and the original one-dimensional
parameterization, neither of which the proof needs.  A forced tournament on
these vertices would add no information—the relation switch is a subgroup
closure, not an antisymmetric dominance relation—so the exact fingerprint is a
saturated lattice rather than scores, cycles, or a Hamiltonian path.  This is a
concrete case where challenging “vertices must be runners” finds the missing
state, while also explaining why a tournament encoding is not cleanly
available.
