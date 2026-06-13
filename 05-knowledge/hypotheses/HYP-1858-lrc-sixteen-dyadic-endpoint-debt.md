---
id: HYP-1858
status: OPEN
source: codex-2026-05-31-S390
related:
  - THM-367
  - HYP-1839
  - HYP-1841
  - HYP-1842
  - HYP-1844
  - HYP-1854
  - HYP-1857
  - HYP-1859
---

# HYP-1858: n=16 LRC counterexamples are killed by dyadic endpoint debt

## Statement

At LRC denominator `n=16`, every primitive full-measure speed set containing
a `16`-gate has a dyadic endpoint-debt leaf.  Equivalently, after building the
finite endpoint/interval protection system, the protection core always peels
to empty before an arithmetic endpoint cycle can close.

Together with the proved unit-gate branch below, this would prove the
`n=16` case:

```text
no speed divisible by 16  -> odd unit points a/16 survive;
some speed divisible by 16 -> dyadic endpoint debt leaf survives.
```

## Proved Branch

If no speed is `0 mod 16`, then every odd unit point `a/16` is a lonely
boundary witness.  Indeed, for every residue `r != 0 mod 16` and every odd
`a`, the distance `||ra/16||` is at least `1/16`.  Residues with odd part can
hit equality; even nonzero residues are farther away.  Only residue `0 mod 16`
kills the unit skeleton.

Thus an open-cover counterexample at `n=16` must contain a `16`-gate.

## Evidence

`lrc_sixteen_sedenion_proof_search_s390.py` runs six exact attacks:

1. unit-skeleton gate lemma;
2. half-turn parity split;
3. structured candidate audits;
4. one-gate replacement scan around the initial tight segment;
5. dyadic endpoint-debt cascade;
6. deterministic forced-gate random stress.

Main exact findings:

```text
initial segment:
  boundary_only, unprotected=8, layer=1

best 8-ladder:
  speeds=(1,8,16,24,32,40,48,56,64,72,80,88,96,104,120)
  gap/th=0.007576
  unprotected=140
  first exposed layer=8

best 4-ladder:
  gap/th=0.015152
  first exposed layer=4

best 2-ladder:
  gap/th=0.030303
  first exposed layer=2
```

The largest-proper-divisor ladder at `n=16` behaves exactly like the pure
doubling prediction: the first exposed layer is the half-dimension `8`.

The one-gate replacement scan tested all `960` sets

```text
{1,...,15} - {r} + {16q}, 1<=r<=15, 1<=q<=64.
```

Every one remained positive-gap.  The closest rows had
`gap/th=0.013672`, still not boundary-only or open-cover.

The forced-gate random stress sampled `48` primitive speed sets with one to
four `16`-gates; all were positive-gap and the closest endpoint audits still
peeled to `coreE=0`.

The local endpoint-cover audit also found that if `v=16` is treated as a
maximum-speed branch and only lower protectors `p<v` are allowed, exactly nine
lower residues are needed to cover all `32` endpoints:

```text
(1, 3, 5, 7, 8, 9, 11, 13, 15).
```

For `v=2,4,8`, lower residues cannot cover all endpoints at all.  This is not
yet a proof for arbitrary largest speeds, but it is the first hard local
shape of the dyadic debt inequality.

## Dyadic Debt Mechanism

For a speed `v`, its endpoints have labels

```text
(16m +/- 1)/(16v).
```

A protector `p` covers such an endpoint exactly when

```text
|p(16m +/- 1) - a*16v| < v
```

for some integer `a`.  The S390 residue audit shows the shape:

```text
super-gate p = 0 mod 16v  protects all endpoints;
half-gate residues        protect the largest non-super chunks;
lower dyadic residues     protect descendant fragments.
```

If `v` is maximal in the speed set, it cannot use a super-gate.  Its endpoint
debt must descend to lower dyadic layers.  The proof target is to turn this
into a positive-divergence inequality: with only `15` speeds, the debt cannot
descend, close the unit layer, and also form a leafless arithmetic endpoint
cycle.

THM-367 now proves the pure dyadic local count law.  For owner `u=2^k` and
protector `p=2^j q mod 16u`, the protected endpoint count is:

```text
p = 0 mod 16u: all 2u endpoints;
j >= k, non-super: 0;
k-j >= 3: 2^(k-2);
k-j = 2: 2^(k-1) only for q = +/-1,+/-3 mod 16;
k-j = 1: 2^k only for q = +/-1 mod 16.
```

The theorem also sharpens the local cover picture: `u=2,4,8` cannot be closed
by lower protectors, while `u=16` has exact lower-cover number `9`, forced by
private endpoints.  Higher pure dyadic owners have self-similar nine-covers,
so the remaining proof target is the global debt-flow inequality recorded in
HYP-1859.

## Predictions

1. Any exact branch-and-bound for `n=16` should terminate fastest when ordered
   by dyadic endpoint layer, not by speed magnitude.
2. A full-measure set with a `16`-gate, if it exists, must show a nonzero
   endpoint core before it can be dangerous; all sampled forced-gate systems
   have `coreE=0`.
3. The useful certificate is a debt-flow inequality across layers
   `1,2,4,8,16`, with half-turn parity making odd-speed protection one-sided.
4. The `n=16` proof should be easier than `n=18` because the odd torsion
   payload is absent; there is no CRT side-channel for the debt to escape.

## Sources

- `04-computation/lrc_sixteen_sedenion_proof_search_s390.py`
- `05-knowledge/results/lrc_sixteen_sedenion_proof_search_s390.out`
- THM-367
- `04-computation/lrc_n16_dyadic_endpoint_formula_s391.py`
- `05-knowledge/results/lrc_n16_dyadic_endpoint_formula_s391.out`
- `07-reflections/lrc-sixteen-sedenion-proof-search-s390.md`
- `07-reflections/lrc-n16-dyadic-endpoint-formalization-s391.md`
- HYP-1841
- HYP-1842
- HYP-1844
- HYP-1854
- HYP-1857
- HYP-1859
