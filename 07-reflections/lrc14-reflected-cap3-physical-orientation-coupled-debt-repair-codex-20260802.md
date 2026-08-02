# LRC(14) reflected cap-3 physical-orientation and coupled-debt repair

**Status:** **PROVED + FINITE-EXACT + VERIFIED**, within the pinned reflected
`k=1` sufficient family of THM-2941.  This is not a physical-survivor census,
does not close the other LRC(14) sectors, and is not a proof of LRC(14).

## Outcome

The corrected exact referee proves

```text
D >= 6 and 2m >= D  ==>  every reflected THM-2941 k=1 packet closes.
```

The proof covers all `3,003` six-label bodies: the prior arbitrary-positive-
level bank covers `2,442`, and the new referee independently rebuilds and
closes the remaining `561`.  Consequently the reflected sufficient-
certificate failure wedge is now confined to

```text
561 bodies, D >= 6, 1 <= m < D/2.
```

The old cap-3 conclusion is not inherited.  Its structural data are pinned as
hostile controls, while its CSP partition, physical orientations, and every
tail obligation are recomputed.

## Inheritance and the exact missing coordinate

The closest proved mechanism was the unaffected cap-`7/4` full-interval
closure.  MISTAKE-347 supplied the canonical hostile example: reversing both
the ratio interval and ordered label pair can encode the same physical order
twice.  The least-used useful sidecar was the failure digraph of a family of
ordered ratio lanes.

For an ordered slot pair `(i,j)`, the script uses

```text
r = Q/P = q_j/q_i.
```

A low lane on `(1,0)` and a high lane on `(0,1)` therefore both cover
`q_0 <= q_1`; they do not form a physical disjunction.  Two exact all-distinct
corner witnesses survive the repeated-level theorem and expose the defect:

- `H=(1,2,3,4,6,12)`, `q=(3,4,5,8,6,9)`: all seven stored Farey lanes fail;
- `H2=(1,3,4,6,8,12)`, `q=(5,4,3,6,7,9)`: both stored side lanes fail.

After exhaustively reselecting every old base-debt lane, `23/26` tail bodies
admit a single whole `[1/3,3]` lane.  Exactly three do not:

```text
H  = (1,2,3,4,6,12),
H2 = (1,3,4,6,8,12),
H3 = (2,3,4,6,8,12).
```

The failure graph for `H` is the acyclic chain `2 -> 1 -> 0`; that for `H2`
is `3 -> 2 -> 1 -> 0`.  Thus no disjunction made only from the available old
full-side lanes can cover all distinct assignments.  In contrast, `H3` has
the inverse low lanes `(0,1)` and `(1,0)`, whose failure edges form a directed
2-cycle.  Equal selected levels are already discharged by the independently
rechecked repeated-level theorem.

## Why this is a selector obstruction, not a boundary runner obstruction

For both acyclic bodies, the referee checks every six-distinct assignment of
six values from `{3,...,9}` with minimum `3` and maximum `9`.  These are the
`3,600` assignments at the first integer cap-3 corner `(m,D)=(3,6)`.

Every one has a strictly positive direct reflected pair certificate.  The
weakest exact margins are

```text
H : 151480151406464/2641051743800701
    at q=(9,8,7,4,5,3), pair=(1,3),
H2: 285989563781330/4777236334436643
    at q=(9,7,8,5,4,3), pair=(2,4).
```

The first obstruction was therefore in the certificate selector, not in the
underlying corner assignments.

## The coupled-debt lemma

Let the selected physical levels be `q_i=P g`, `q_j=Q g`, with `(P,Q)`
reduced, and put `s=min(q_i,q_j)`.  On a closed ratio lane `r=Q/P in [u,v]`,
the selected-pair maximum satisfies

```text
pair_max >= alpha s,

alpha = 1/v,  if v <= 1,
        u,    if u >= 1,
        1,    if u <= 1 <= v.
```

Indeed, below one the maximum is `q_i=s/r >= s/v`; above one it is
`q_j=rs >= us`; a lane crossing one has the tautological lower bound `s`.
The cap-3 cone says `q_max <= 3m`, so

```text
m >= ceil(alpha s/3).                                      (1)
```

The singleton debt is nonincreasing in the available global minimum.  Hence
the tail may charge debt at

```text
ell(s) = max(3, ceil(alpha s/3)),
```

instead of uniformly charging the worst level-three debt.  This is the
missing coordinate: the ratio lane controls not only phase transport but also
how much singleton debt remains possible.

## Exact monotonicity

For a lane with endpoint numerator bound `N`, body ruler `L`, and first label
`a`, the lower envelope is

```text
E(s) = 1/49 - 12/(49s^2) - N/(L-a/s) - debt(ell(s)).
```

The exact increment from `s` to `s+1` is the sum of three nonnegative terms:

```text
12/49 (1/s^2 - 1/(s+1)^2),
N a / ((Ls-a)(L(s+1)-a)),
debt(ell(s)) - debt(ell(s+1)).
```

The first two are strictly positive and the third is nonnegative.  Thus one
positive threshold proves the entire lane suffix, including across ceiling
jumps in `(1)`.

For each reduced head channel the verifier separately fixes a primitive
skeleton, checks every earlier common scale directly, and proves the infinite
suffix using decreasing phase displacement, decreasing transport displacement,
and nonincreasing singleton debt.  The favourable inverse transport factor is
discarded only after the primitive bracket is positive.

## Assignment-complete tail atlas

The repaired atlas has `35` lanes and `1,600` exact head controls.

- `23` bodies use one whole `[1/3,3]` lane, for `736` heads.
- `H3` uses inverse low lanes `(0,1)` and `(1,0)`, for `56+20=76` heads.
- `H` uses six adjacent closed intervals, all on the same physical pair
  `(0,1)`:

  ```text
  [1/3,3/8], [3/8,3/7], [3/7,1/2],
  [1/2,2/3], [2/3,1], [1,3],
  ```

  with first tail scales `56,35,26,23,19,10`, for `686` heads.
- `H2` uses four adjacent closed intervals on `(0,1)`:

  ```text
  [1/3,2/5], [2/5,3/5], [3/5,1], [1,3],
  ```

  with first tail scales `10,10,10,9`, for `102` heads.

The fixed-pair interval chains cover the literal physical ratio, rather than
an orientation-renamed copy.  Closed endpoint overlap is harmless.  The
weakest transported tail margin is

```text
508374102558227956/234246433355178759099585
```

on `H`, pair `(0,1)`, channel `(81,29)`.  The largest directly checked prefix
has six common scales.

## Recomputed constrained layer

The cap-3 corner transport maximum is `35/164`, attained on
`H`, pair `(0,5)`, at `(D,m)=(6,3)`.  Rebuilding from the complete `8,796`-
channel primitive bank gives `5,666` constrained edges and the exact split

```text
561 residual bodies
  = 282 CSP-closed
  + 279 traps
  = 282 + 253 constrained policies + 26 unconstrained tails.
```

Forward and reverse CSP search orders agree.  A constrained policy stores an
unordered gain `rho=max(q_i,q_j)/min(q_i,q_j)`.  The repaired audit expands
every allowed `rho` to both ordered channels `(P,Q)` and `(Q,P)`.  All `3,062`
physical-orientation controls are strictly positive; the weakest margin is

```text
4389035249640369833/461756458190054127555
```

at body `(3,4,6,8,9,12)`, pair `(0,1)`, channel `(9,26)`.  This explains why
the old orientation error was localized to the unconstrained tails.

## Reproduction and scope

Run

```bash
python3 04-computation/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.py
python3 -O 04-computation/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.py --output /tmp/cap3-opt.out
cmp 05-knowledge/results/lrc14_j7_reflected_half_cone_orientation_enrichment_scout_thm2941.out /tmp/cap3-opt.out
```

Ordinary, optimized, and stored transcripts are byte-identical.  LF-normalized
hashes are

```text
source   d53b9ebb76930fc6eb4cdc697c11f57413747e834e0ee4a11ea81ee47c884398
output   a5cfde4f3ba447c92d5486017bddd484a0cebfad87330cd19c25d307327a5d87
semantic 1a8a6cac81420388034c849ca0ba852a9b87416c85eeedcb12e1cbafea572f5a
```

This proof restores the cap-3 reflected sufficient-family cone and supersedes
the *conclusion* of the quarantined half-cone referee.  The old source remains
valuable as a hostile control and correction-lineage dependency.  No claim is
made about physical survivors outside the reflected THM-2941 sufficient
family, the six-body/seven-tail rung, the projected `k=2,3` sectors, or
LRC(14) itself.
