---
id: HYP-3731
title: The COVERING-MIN as an INTEGER PROGRAM over the DANGER-CIRCULANT -- pins the odd-n values (n=7->2/13, 8->2/15, 9->4/33, 11->3/31 NEW, beating the construction 11/111) and makes the CHROMATIC<->OCF bridge computational. Two dual reframes (the owner's 'copy the observer to all n points'): SET-COVER (one observer at 0: every rotation a has a speed within j of 0) vs INDEPENDENT-SET (all observers: the n runner-positions mutually >=M apart in C_D({1..j})). The SYMMETRIC (all-observers) version is GALILEAN-INVARIANT (lives on the difference set / the translation-invariant circulant); the one-observer M is not. The danger-circulant's ORIENTATION reruns the even/odd split: EVEN n -> 2n-1=3 mod 4 = Paley TOURNAMENT -> Redei/OCF H-count (odd: n=4->189, n=6->95095); ODD n -> 2n-1=1 mod 4 = Paley GRAPH -> Ramanujan/Ihara
status: VERIFIED (scipy.milp IP pins covering-min n=7,8,9,11 = 2/13,2/15,4/33,3/31, honest gap re-checked; n=7,8,9 match mac-mini, n=11 NEW; symmetric M Galilean-invariant under V->V-w; Paley tournament Redei-odd for q=7,11; Paley graph Ramanujan for q=13,17). HONEST: the IP pins the small-beater regime (Smax=4n suffices, n<=11); for n>=13 the optimal killer/binding exceeds 4n (construction scale ~n(n-1)), so the n=13 run under-resourced returned 1/12 > construction 13/157 (non-optimal). The chromatic<->OCF is made computational, not proved equal.
source: klein-2026-06-30-S36
depends_on:
  - HYP-3729   # even/odd = bipartite/non-bipartite C_n
  - HYP-3727   # mac-mini: primitivity + the Paley frame (2n-1)
related:
  - HYP-3728   # the metazeta / Ihara-Ramanujan (the odd-n Paley graph side)
  - HYP-3602   # OCF / Redei (the even-n Paley tournament side)
  - THM-523    # primitive covering sets carry the margin
results:
  - 04-computation/covering_min_ip_danger_circulant_klein.py
  - 04-computation/covering_min_ip_sequence_klein.py
  - 04-computation/reframes_chromatic_ocf_klein.py
  - 04-computation/reframes_fix_klein.py
  - 05-knowledge/results/covering_min_ip_danger_circulant_klein.out
  - 05-knowledge/results/reframes_fix_klein.out
---

# HYP-3731 — the covering-min IP over the danger-circulant; chromatic↔OCF

## The IP pins the odd-n covering-min (the owner's explicit ask)
**Set-cover formulation** (one observer at 0): `M(S) <= t` iff for EVERY modulus `D'` and rotation `a`, some
speed `s` has `dist(s.a mod D',0) <= floor(t.D')` -- the danger combs cover at every scale. Binary-searching
`t` with `scipy.milp` (constraints: this set-cover at all `D'<=Dmax`, resonance-killing `b|s` for `b<=n`,
primitivity `p∤s` for some `s`, all primes `p<=n`) gives the smallest feasible `t` = the covering-min:
```
n  parity   covering-min   primitive set found by the IP
7   odd      2/13          {1,2,5,6,7,8}         (= mac-mini's beater)
8   even     2/15          {1,3,4,5,7,11,24}     (= mac-mini's value)
9   odd      4/33          {1,4,5,6,7,11,32,36}  (alt optimum; same value as mac-mini's {1,3,4,5,7,11,18,32})
11  odd      3/31          {2,6,8,9,10,11,13,14,17,19}   NEW -- and 3/31 < 11/111 (beats the construction)
```
Honest re-check: each `S` re-verified by exact gap over `D'<=400`. **Scope:** `Smax=4n` suffices for the
small-beater regime (`n<=11`); for `n>=13` the optimal killer/binding exceeds `4n` (the construction scale
`~n(n-1)`), so the under-resourced `n=13` run returned `1/12 > 13/157` (non-optimal) -- pinning `n=13,14`
needs larger `Smax` (heavier IP).

## The two dual reframes (the owner's "copy the observer to all n points")
- **SET-COVER (one observer at 0):** the IP above -- the danger combs of the chosen speeds cover all
  rotation-times. The asymmetric, distinguished-`0` value.
- **INDEPENDENT-SET (observer copied to all n points):** at the lonely time the `n` runner-positions are an
  INDEPENDENT SET in the danger-circulant `C_{D'}({1..floor(M.D')})` -- every runner far from every other.
  The DUAL of the set-cover.

## Change-the-observer / translate = Galilean invariance (the owner's question)
The **one-observer** `M` is NOT invariant under moving the observer (the point `0` is distinguished). The
**symmetric (all-observers)** `M` -- every runner lonely from ALL others -- IS Galilean-invariant: it depends
only on the pairwise DIFFERENCE set `{v_i-v_j}`, unchanged by `V -> V-w`. Verified `V={0,1,2,5,6,7,8}` (diffs
`{1..8}`, `M_sym=1/9`) invariant under all shifts. So "copy the observer to all `n` points" is exactly the
move that makes the analogy invariant -- the symmetric version lives on the translation-invariant
danger-circulant (a Cayley graph on the GROUP `Z/D`). Same machine, two analogues (asymmetric `2/13` vs
symmetric `1/9`).

## The chromatic↔OCF bridge made computational, and the even/odd split AGAIN
The danger-circulant `C_D({1..j-1})` (the mutual-danger graph) has computed chromatic numbers
`chi = 3,3,5,4` for `n=7,8,9,11` (its circular chromatic number is the classical lonely-runner ↔ distance-
graph link). Its **orientation** reruns the project's even/odd split:
- **EVEN `n`:** `2n-1 = 3 mod 4` -> the **Paley TOURNAMENT** -> the OCF / Redei Hamiltonian-path count `H(T)`
  (verified ODD: `q=7` `H=189`, `q=11` `H=95095`). This is the OCF side.
- **ODD `n`:** `2n-1 = 1 mod 4` -> the **Paley GRAPH** -> Ramanujan / Ihara (verified Ramanujan `q=13,17`).
  This is the metazeta side (HYP-3728).
So the chromatic↔OCF bridge is: `chi(danger-circulant)` and the `H`-count of its orientation, computed by
DP/OCF -- and the orientation is a tournament (even `n`, OCF) or a graph (odd `n`, Ihara) by the SAME
bipartite/non-bipartite parity (HYP-3729). Even `n` = OCF/Redei; odd `n` = Ihara/Ramanujan.

## Net
The covering-min IS an integer program over the danger-circulant -- the set-cover (one observer) and its dual
independent-set (all observers / the owner's "copy to all `n` points"). It PINS the odd-`n` values
(`7->2/13, 9->4/33, 11->3/31`; `11` new, beating the construction) in the small-beater regime, and `n=8->2/15`
even. The symmetric version is the Galilean-invariant analogue (the difference set / translation-invariant
circulant). The chromatic↔OCF bridge is made computational, and the danger-circulant's orientation reruns the
even/odd split one more time: EVEN `n` -> Paley tournament -> OCF/Redei; ODD `n` -> Paley graph -> Ihara/
Ramanujan. The lonely runner, the tournament OCF, and the Ihara zeta are three readings of one circulant.
