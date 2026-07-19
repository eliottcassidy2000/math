---
id: THM-1151
title: Exactly two thirteen-carriers reduce to thirty normalized switch obstructions
status: PROVED — elementary integer/rational classification.  For two normalized carriers u<v above a core maximum 1<=M<=12, the translated thirteen-grid switch of THM-1136 exists except for exactly 30 pairs: 8<=M<=12, u=M+1, and 13M+14<=v<=15M-1.  Thus every r=6 row with exactly two 13-divisible killers is lonely outside this explicit list.  A dependency-free exact referee replays the constructive certificate and brute classification on 2,081,164 integer rows
source: codex-2026-07-18-S67 continuation
depends_on: [THM-1136]
related: [THM-1135, THM-1137, THM-1145]
script: 04-computation/lrc14_r6_two_carrier_switch_referee_codex_S67.py
output: 05-knowledge/results/lrc14_r6_two_carrier_switch_referee_codex_S67.out
---

# THM-1151 — the exact two-carrier switch classification

Let `P` be a nonempty subset of `{1,...,12}` and put `M=max(P)`.  Suppose an
`r=6` clustered row has exactly two killers divisible by `13`.  Write them

```text
z=13u,  h=13v,       M<u<v.                           (1)
```

For a source `s in {u,v}` and an integer switch

```text
1<=c<=floor(s/M),                                    (2)
```

call `(s,c)` **carrier-safe** when

```text
dist(cu/s,14Z)>=1  and  dist(cv/s,14Z)>=1.            (3)
```

This is exactly THM-1136(B), after cancelling the common factor `13`.  If a
carrier-safe switch exists, the translated thirteen-grid argument absorbs the
four noncarrier killers and produces a `1/14`-lonely time.

## The classification

> **Two-carrier switch theorem.**  A carrier-safe switch exists unless
>
> ```text
> 8<=M<=12,
> u=M+1,
> 13M+14<=v<=15M-1.                                  (4)
> ```
>
> Conversely every triple in (4) has no carrier-safe switch.  There are
> exactly
>
> ```text
> sum_{M=8}^{12} 2(M-7)=2+4+6+8+10=30               (5)
> ```
>
> exceptional triples.

### Constructive proof

First try source `u` and `c=1`.  The source itself is safe at equality.  If
the other carrier is also safe, we are done.  Otherwise, with `r=v/u`,

```text
dist(r,14Z)<1.
```

Since `r>1`, this puts `r` in one of the open bands

```text
14m-1<r<14m+1,       m>=1.                            (6)
```

In particular `r>13`.  Put `q=ceil(r)` and change the source to `v`.

If `q` is not divisible by `14`, choose `c=q`.  The source is safe because
`dist(q,14Z)>=1`, while

```text
1<=q/r<2                                                (7)
```

makes the smaller carrier safe.  The cap is automatic: `u>=M+1` and
`r>13>=M+1` give

```text
v/M=r u/M >= r+r/M > r+1 >= q.                        (8)
```

Now suppose `14|q`.  Choose `c=q+1` whenever the cap permits it.  Again the
source is safe at distance one and `1<(q+1)/r<2` makes the smaller carrier
safe.  If `r>=2M`, then

```text
v/M >= r+r/M >= r+2 >= q+1,                           (9)
```

so the cap permits this switch.  Inequality `r>=2M` holds throughout every
band with `m>=2`, because then `r>27>24>=2M`; it also holds in the `m=1`
band when `M<=6`.

It remains to inspect `m=1` and `M>=7`.  Here `q=14`.  If `u>=M+2`, then
the integrality of `v>13u` gives

```text
v>=13u+1>=13M+27>=15M,                                (10)
```

where the last inequality uses `M<=12`.  Thus `c=15` is allowed.  The same
is immediate if `v>=15M`.  Consequently failure is possible only when

```text
u=M+1,       13u<v<15M.                               (11)
```

After using integrality this is exactly (4), and the interval is nonempty
exactly when `M>=8`.

For the converse, assume (4).  Then

```text
floor(u/M)=1,       floor(v/M)=14,       13<v/u<14.    (12)
```

The only small-source switch is `c=1`, and the larger carrier lies less than
one unit below `14Z`, so it obstructs.  At source `v`, every `1<=c<=13`
satisfies `0<cu/v<1`, so the smaller carrier obstructs.  The last possible
switch `c=14` is divisible by `14` and the source obstructs itself.  No
switch exists.  This proves both directions.

## Consequence for the r=6 frontier

> **Corollary.**  Every clustered `r=6` row with exactly two `13`-divisible
> killers is `1/14`-lonely unless its normalized carrier pair belongs to the
> thirty rows (4).

This follows immediately from THM-1136(B).  Notice what has disappeared:
the positions, magnitudes, and residue collisions of all four noncarriers.
Each excludes at most two of the twelve observer multipliers, so their union
cannot exhaust the chart.  The only remaining obstruction is the exact
two-carrier band-and-cap interaction above.

The theorem does **not** assert that the thirty carrier pairs are LRC
counterexamples.  It says precisely that this particular translated-grid
certificate is impossible there.  Those thirty pairs are now the finite
target for a second chart or a component argument.

## Structural viewpoint and Tournament Analysis

The two carrier vertices with the speed-order gauge form the unique
two-vertex tournament: score multiset `{0,1}`, no directed cycle, two
singleton SCCs, and one Hamiltonian path.  That quotient is completely
uninformative.  It cannot distinguish a certified pair from any of (4).

The faithful vertices are the proof obligations `(source,c)`, with labelled
edges to whichever carrier violates (3).  The source-dependent cap (2) is
part of the object.  On an exceptional row, source `u` has one obligation,
source `v` has fourteen, the first thirteen large-source obligations point
to `u`, and the last is a self-obstruction at `v`.  This labelled bipartite
or directed-incidence object preserves exactly the LRC certificate; a naked
tournament on runners, carriers, or ratios destroys the band positions and
the cap.

Candidate vertex sets explicitly challenged here are runners, carriers,
carrier ratios, switches, danger bands, and proof obligations.  The switch
obligations are the minimal faithful choice.

## Exact replay

The referee clears every denominator and checks safety by integer remainders
modulo `14s`.  It exhausts a box containing all thirty exceptions and several
higher `14m` bands, compares the closed formula against literal switch search,
and independently checks the constructive branch certificate on every row.
Ordinary and `PYTHONOPTIMIZE=1` executions are required to be byte-identical.

```text
04-computation/lrc14_r6_two_carrier_switch_referee_codex_S67.py
SHA-256 367e25402a52ccf7e3476b7cdd982641eac63de798360b6713c722b017cbf6a9

05-knowledge/results/lrc14_r6_two_carrier_switch_referee_codex_S67.out
SHA-256 b1feaa6ac2bbcb187cb89f14a3c6bf9ab2298879975bae97e114f0b0659ac6d3
```
