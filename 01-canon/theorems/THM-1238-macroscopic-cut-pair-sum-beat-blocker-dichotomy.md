---
id: THM-1238
title: THE MACROSCOPIC-CUT PAIR-SUM BEAT DICHOTOMY — a Kakeya path edge emits a nonzero mixed-clock blocker unless it is locked in one negative-curvature beat cell
status: PROVED (all-scale two-speed beat theorem; Kakeya cut composition; exact singleton residual and infinite parity guardrail; dependency-free referee; sorry-free Lean arithmetic core)
source: codex-2026-07-19-S78 continuation with tangent-stalk agent
depends_on: [THM-1156, THM-1217, THM-1219, THM-1241, THM-1236]
related: [THM-1240, THM-1242, HYP-7870]
script: 04-computation/lrc14_macroscopic_cut_pair_sum_beat_dichotomy_thm1238.py
output: 05-knowledge/results/lrc14_macroscopic_cut_pair_sum_beat_dichotomy_thm1238.out
formalization: 04-computation/lean/TournamentH7/TournamentH7/LRCMacroscopicCutPairBeat.lean
script_sha256: 28632f623004cc767499e9752a3f2552f075e7ddac3c060b340564d1e79b800b
output_sha256: 03b83b0c078d4ba033b705738a91c07fa93cafb60c7e886c2c46f9906a671cd9
formalization_sha256: 6c3d3029cf3b91ee18274cfd7e60d96d000c7f98d0e90aefca99acd42dcd360d
---

# THM-1238 — macroscopic cut to pair-sum blocker

## 1. Two-speed beat dichotomy

Let `c<x<y` be positive integers, let

```text
G=G_k(c)=[(14k+1)/(14c),(14k+13)/(14c)],
q=x+y,
P_q(k)={p in Z:p/q in G}.                              (1)
```

The pair sum satisfies `q>2c`, so `P_q(k)` is a nonempty consecutive integer
block.  For `p in P_q(k)`, put

```text
D(p)=min(xp mod q,q-(xp mod q)).                        (2)
```

Since `y=-x mod q`, the two defining speeds have equal circle distance at
`t=p/q`:

```text
||xp/q||=||yp/q||=D(p)/q.                              (3)
```

The signed active-pair curvature is exactly

```text
14D(p)-q.                                              (4)
```

If either

```text
q>=7c/3                    or                    y>6x, (5)
```

then there is `p in P_q(k)` with

```text
14D(p)>=q.                                             (6)
```

Thus both defining speeds are safe at `p/q`.  If `x,y` are two members of a
six-comb cover of `G`, one of the other four speeds must be dangerous there:
(6) emits a literal third-blocker obligation.

If (6) fails for every beat in the block, the residual is exact:

```text
P_q(k)={p},
2c<q<7c/3,
y<=6x,
14D(p)-q<0.                                            (7)
```

Moreover equality in (6) can occur only when

```text
14|q.                                                  (8)
```

This is precisely THM-1156's zero-curvature seam/third-support branch.  When
`14` does not divide `q`, every supplied blocker in (6) sits over a strictly
positive active edge and therefore covers an open neighborhood of the beat.

## 2. Proof

The scaled interval `qG` has length

```text
|qG|=6q/(7c)>12/7,                                     (9)
```

which proves nonemptiness of `P_q(k)`.  First suppose `y<=6x`.  The normalized
step of `xp mod q` obeys

```text
1/7<=x/q<1/2.                                          (10)
```

Two consecutive residues therefore cannot both lie in the strict symmetric
danger arc of length `1/7`.  If also `q>=7c/3`, (9) is at least two, so the
block contains two consecutive integers and one satisfies (6).

Now suppose `y>6x`.  Put

```text
R=q/(7x)>1.                                            (11)
```

Here `0<x/q<1/7`.  A consecutive run of defining-danger residues has length
at most `ceil(R)`: after `ell-1` positive steps inside one lifted danger arc,

```text
(ell-1)x/q<1/7.                                        (12)
```

But `x>c` gives

```text
|qG|=6q/(7c)>6R,                                       (13)
```

so `P_q(k)` contains at least `ceil(R)+1` integers.  It cannot be wholly
dangerous, proving (6).

If both suppliers fail, then `q<7c/3`, `y<=6x`, and (9) lies strictly between
one and two.  The beat block can contain only one or two consecutive
integers.  The latter case would trigger (10), so failure of (6) forces the
singleton (7).  Finally, `D(p)` is integral, so `14D(p)=q` immediately gives
(8).

## 3. Mixed-clock lift and the Kakeya cut

At a supplied beat, the reduced defining period is

```text
Q=q/gcd(x,y)>=3.                                       (14)
```

For a full six-speed packet `d1,...,d6`, let

```text
d0=gcd(q,d1,...,d6),       L=q/d0.                     (15)
```

Then `Q|L`.  The safe beat has nonzero residue modulo `L`: if `L|p`, every
speed divisible by `d0`, including `x`, would be integral at `p/q`, contrary
to (6).  Thus branch (6) is a genuinely nonzero mixed-clock blocker, not the
common-zero residue seen by a collapsed beat mask.

THM-1241 chooses an adjacent path edge

```text
x=dj<y=d_(j+1),
y-x>d6/210.                                            (16)
```

Applying the theorem to that pair produces a finite cut-index dichotomy:

1. a nonzero positive/zero-curvature beat with a third blocker; or
2. the one-cell negative-curvature phase lock (7).

If at least two speeds cross the full-lap threshold `7c/6`, the terminal pair
`(d5,d6)` has `d5+d6>7c/3`, so branch 1 is forced independently of which edge
THM-1241 selected.  Only the `r=1` suffix branch can retain the universal
singleton obstruction.

## 4. The singleton branch is real

For every `m>=1`, set

```text
c=420m+1,                 k=210m,
(d1,...,d6)=
(428m+2,440m+2,452m+2,464m+2,476m+2,500m+2).           (17)
```

The slow gap is centered at `1/2`; `c` is odd and every fast speed is even.
The packet satisfies the THM-1241 macroscopic requirements:

```text
d5<7c/6<d6,
sum_i(d6-di)=240m>d6/14,
max_j(d_(j+1)-d_j)=24m>c/180.                          (18)
```

For every one of the fifteen fast-fast pairs, `q=x+y` is even, `qG` is
centered at the integer `q/2`, and

```text
1<|qG|<2.                                              (19)
```

Hence `P_q(k)={q/2}`.  Since `x,y` are even,

```text
x(q/2)/q=x/2 in Z,
y(q/2)/q=y/2 in Z,                                    (20)
```

so `D=0` and the curvature is `-q`.  Every reduced period in (14) is still at
least three.  This is not a cover construction; it is an infinite exact
obstruction to deducing a positive active edge from macroscopic speed
separation plus nontrivial quotient periods alone.

The lesson is structural: the fast-fast Hamiltonian path is incomplete.
Positive curvature can be exported to the two carrier-boundary roles in
THM-1219's augmented path.  THM-1240 captures those carrier spokes directly
and avoids the parity lock.

## 5. Tournament and carrier audit

The speed-order tournament merely selects the path edge (16).  It destroys
the beat block, determinant, clock residue, curvature sign, and blocker
identity.  The faithful vertices are adjacent path edges lifted to obligations

```text
(q,Q,L,P_q(k),p mod L,14D-q,blocker label).             (21)
```

This quotient preserves exact defining-pair safety and third-support truth at
the beat.  It still destroys off-beat coverage.  The parity family challenges
the assumption `macroscopic cut = active edge`: the cut can coexist with a
common-zero singleton on every fast-fast sum clock.

For Tournament Analysis one may orient blocker incidence rather than speed
order once branch 1 is present.  That directed relation can have cycles and
is developed on the unconditional carrier spokes in THM-1240.  A transitive
runner tournament cannot see it.

## 6. Verification and scope

The exact referee checks `1,070,405` two-speed/gap rows, including both
supplier regimes, `7,179,634` positive safe-beat occurrences, `10,856`
zero-curvature occurrences, and every forced singleton.  It independently
checks all fifteen fast-fast pairs of (17) for `m=1,...,1000`.  Normal and
optimized Python outputs are byte-identical.

The Lean module kernel-checks the moderate/far step bounds, scaled-gap length
regimes, exact seam divisibility, third-blocker logic, even half-beat
factorization, and the symbolic parity-family Kakeya ledger.  The lifted-arc
run bound (12) remains the explicit paper lemma; there are no proof
placeholders.

Frozen hashes are

```text
source         28632f623004cc767499e9752a3f2552f075e7ddac3c060b340564d1e79b800b
output         03b83b0c078d4ba033b705738a91c07fa93cafb60c7e886c2c46f9906a671cd9
formalization  6c3d3029cf3b91ee18274cfd7e60d96d000c7f98d0e90aefca99acd42dcd360d
```

THM-1238 converts the Kakeya cut into a precise blocker/phase-lock alternative.
It does not make the blocker hypergraph inconsistent and does not prove
LRC(14).
