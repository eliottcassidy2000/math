---
id: THM-589
title: The labeled-tournament Hamiltonian-path-count variance has an EXACT even-run closed form -- A_n(2) := sum_{pi'} c(pi')2^{j(pi')} = sum over edge-subsets of the integer path [1..n-1] with ALL maximal runs of EVEN length of 2^{#runs}(n-|S|)! (odd overlap runs CANCEL by orientation parity) -- giving CV(H)^2 = A_n(2)/n! - 1 with the precise rate CV(H)^2 ~ 2/n (n*CV^2 -> 2), the dominant fluctuation being a single length-2 overlap run
status: PROVED (inclusion-exclusion on descending successions + the per-run factor (1+(-1)^m)=2[m even]); VERIFIED -- the even-run form equals the permutation form (n<=9) and direct labeled-tournament enumeration (n<=5); n*CV^2 -> 2 to n=320 (1.92654, 1.98548, 1.99676, 1.99923, 1.99981, 1.99995). A_n(2)=1,2,8,32,158,928,6350,49752,439670,4327904,46963358 is NOT in OEIS (new).
source: klein-2026-06-29-S5
depends_on:
  - HYP-3560   # mac-mini: CV(H)^2 -> 0 via Poisson(1) adjacencies (this gives the EXACT structure + rate)
related:
  - THM-587    # the signed cycle index: H lives on the arc-hypercube; orientation parity is its signed structure
  - THM-588    # no linear invariant, one quadratic => proof effort is purely 2nd-moment (the cyclicity)
  - THM-579    # the LRC floor gatekeeper CV(N_R)^2: the witness-side analog (S_n-collapse) of the floor-side (Z_14-collapse)
  - THM-580    # the LRC 2-adic descent (even survives, odd peels): the floor-side mirror of this even-run parity
  - HYP-3554   # CV(N_R)^2 is set-dependent (Z_14 collapse); CV(H)^2 is the clean set-independent S_n collapse
results:
  - 04-computation/cvH_even_run_closed_form_klein.py
  - 05-knowledge/results/cvH_even_run_closed_form_klein.out
---

# THM-589 — the even-run closed form for CV(H)^2, and the 2/n concentration rate

## Statement

Over labeled tournaments on `[n]` (each arc oriented uniformly), let `H` = number of directed
Hamiltonian paths. `E[H] = n!/2^{n-1}`. Fixing a reference path and using relabeling symmetry,

```
CV(H)^2 = E[H^2]/E[H]^2 - 1 = (1/n!) sum_{pi' in S_n} c(pi') 2^{j(pi')} - 1
```

where `c(pi')=1` iff `pi'` has no DESCENDING consecutive-integer adjacency and `j(pi')` = number of
ASCENDING ones (the owner's identity; mac-mini HYP-3560). Define `A_n(2) := sum_{pi'} c(pi')2^{j(pi')}`.

**Even-run closed form.** Identify the consecutive-integer adjacencies with edges of the path graph
`1-2-...-n` (`n-1` edges). Then

```
A_n(2) = sum_{ S subseteq {edges}, every maximal run of S has EVEN length } 2^{#runs(S)} (n - |S|)!.
```

**Rate.** `CV(H)^2 = sum_{t>=2} W(n-1,t)(n-t)!/n!` (sum over even-run subsets of size `t`, weighted
`2^{#runs}`), whose leading term (a single length-2 run, ~`n` positions, weight 2, factor `(n-2)!/n!`)
gives

```
CV(H)^2 = 2(n-2)/(n(n-1)) + O(1/n^2) ~ 2/n,    so   n*CV(H)^2 -> 2.
```

`H` concentrates: `std(H)/E[H] ~ sqrt(2/n)`.

## Proof of the even-run form

By binomial inversion (forcing a set `L` of consecutive-adjacency successions glues them and leaves
`(n-|L|)!` permutations), the bivariate succession sum is
`sum_pi x^{asc} y^{desc} = sum_{L} (x-1)^{a(L)}(y-1)^{d(L)}(n-|L|)!`, where `L` ranges over orientations
of edge-subsets of the path that form directed linear forests (each maximal run oriented monotonically:
ascending or descending). At `(x,y)=(2,0)`: `(x-1)=1`, `(y-1)=-1`, so
`A_n(2) = sum_L (-1)^{d(L)}(n-|L|)!`. Grouping by the underlying edge-subset `S` and summing over the two
orientations of each maximal run of length `m`, the per-run factor is `1 + (-1)^m`, which is `2` if `m`
is even and `0` if `m` is odd. Hence only edge-subsets whose every maximal run is even survive, each
contributing `2^{#runs}(n-|S|)!`. ∎

## The parity mechanism (why it matters)

The cancellation is the content: when two Hamiltonian paths share a run of `m` consecutive-integer
arcs, the ascending and descending orientations of that run contribute with opposite sign, so an
ODD-length shared run cancels and only EVEN-length overlap survives. This is a Redei-flavored
orientation parity, and it is the **`S_n`-side (witness) mirror of the LRC 2-adic descent** (THM-580,
"even survives, odd peels"): the same even/odd-2-adic split, here on the Hamiltonian-path second moment
rather than the covering floor.

## Integration with the metagraph (THM-587/588)

`H` is a degree-`(n-1)` function on the arc-hypercube `Q_{C(n,2)}` (THM-584): `E[H]` is its level-0
(constant) Walsh component, and `Var(H)` is its higher-level energy. The orientation-parity cancellation
above IS the signed (hyperoctahedral) structure of the signed cycle index (THM-587) restricted to `H`.
THM-588 (no level-1 invariant, a single level-2 = cyclicity) explains why the proof effort sits at the
second moment; THM-589 evaluates that second moment in closed form. So `CV(H)^2 ~ 2/n` is the exact
witness-side companion of THM-588's spectral picture.

## Relation to the LRC floor (HYP-3554, HYP-3560, THM-579)

`CV(H)^2` (witness side) and `CV(N_R)^2` (THM-579 floor gatekeeper) are BOTH reference-collapsed
pair-overlap second moments:
- `CV(H)^2`: collapse two Hamiltonian paths by `S_n` -> a single reference-path permutation sum ->
  CLEAN, SET-INDEPENDENT, `~2/n -> 0` (this theorem). The collapse to one reference path is exactly the
  kind of set-independent reduction the LRC floor wants (mac-mini's Gamma_0(N) goal, HYP-3553).
- `CV(N_R)^2`: collapse two 14-sheets by `Z_14` -> a per-set autocorrelation sum -> SET-DEPENDENT,
  unbounded as `m_R -> 0` (HYP-3554). The difference is the vanishing-fiber corner that the transitive
  `S_n` collapse does not have.

So `H`'s concentration is the finite, fully-solved rehearsal (mac-mini HYP-3560) AND its even-run parity
names the mechanism the floor side must reproduce: a set-independent collapse plus an even/odd-2-adic
cancellation.

Caveat (MISTAKE-033): "consecutive-integer adjacency" here is relative to the FIXED reference labeling;
the cancellation is the orientation (complement) parity, the antipodal `Z_2` of THM-584.
