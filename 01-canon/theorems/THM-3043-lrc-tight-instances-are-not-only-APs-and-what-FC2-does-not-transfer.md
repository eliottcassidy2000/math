---
id: THM-3043
title: "LRC: measure-zero is tightness, the finite tight census is not only APs, and its witnesses have denominator n+1"
status: >
  FINITE-EXACT + VERIFIED-EXACT + CORRECTION (MISTAKE-343) + an honest
  NEGATIVE on transferring the FC(2) proof.
  (R1) MAX-NOT-MEAN, instantiated. mu(Safe) = 0 does NOT mean the danger arcs cover: the
  witness condition ||v t|| >= 1/(n+1) is CLOSED, so the safe set can be a nonempty finite
  set of POINTS. Every measure-zero instance found below is TIGHT -- LRC holds there with
  equality -- and none is a covering failure. Any census that reads "mu = 0" as "covering"
  is misclassifying tight instances, which is exactly the repo's recorded MAX-not-MEAN trap.
  (R2) THE TIGHT FAMILY IS STRICTLY LARGER THAN DILATED APs. Exact census over all speed
  sets from [1,B]: n=3, B=14 gives 4 measure-zero sets, all dilated APs; but n=4, B=14 gives
  5, of which (1,3,4,7) (and its dilate (2,6,8,14)) is NOT a dilated AP; and n=5, B=12 gives
  3, of which (1,3,4,5,9) is NOT. So an inverse theorem phrased as "only dilated APs are
  extremal" is FALSE as stated; the target class contains sporadic members.
  (R3) EVERY TIGHT WITNESS SITS AT t = a/(n+1). In all six tight instances checked, the safe
  points have reduced denominator n+1: (1,2,3) -> 1/4, 3/4; (1,2,3,4) and
  (1,3,4,7) -> a/5; (1,2,3,4,5,6) -> a/7; (1,3,4,5,9) -> 1/6, 5/6; and
  (1,2,3,4,5,6,7) -> a/8 for a in {1,3,5,7}.
  The reduced witness denominator is exactly n+1, not merely "some q".
  (R4) THE FC(2) TRANSCENDENCE MECHANISM DOES NOT TRANSFER. FC's engine is that a
  counterexample pins a period to a RATIONAL value and transcendence contradicts it
  (THM-3031/THM-3039). The LRC analogue is quantisation: mu(Safe) lies in (1/N)Z with
  N = (n+1) lcm(v), verified exactly. But the pinned value is 0, and 0 carries no arithmetic
  rigidity, so no contradiction is available; and the gap principle "mu = 0 or mu >= 1/N" is
  far too weak in practice (observed mu/(1/N) = 30, 58, 66, 11020). What DOES transfer is
  only the ARCHITECTURE: exclude non-rigid configurations, reduce to a rigid family, kill
  that family separately -- which is LRC Route B, and (R2) says its rigid family was
  mis-specified.
source: death-star-2026-08-01-coinC2
depends_on: []
related:
  - THM-1123
  - THM-1132
  - THM-1134
  - THM-3031
  - THM-3039
external:
  - "Lonely Runner Conjecture; the classical tight instances for small n."
script: 04-computation/lrc_tight_instance_census_thm3043.py
output: 05-knowledge/results/lrc_tight_instance_census_thm3043.out
script_sha256: f2bf9680be73400bc075268890ea90d2cc37f2c544400990fc81fac02c9c3b30
output_sha256: b1a47938475b504cf1ba8fe805f3c842db56649a55819d7997ebde7c42201c98
hash_basis: LF-normalized bytes
---

# THM-3043 -- tight is not covered, and the tight family is bigger than APs

## 1. The correction (R1)

For integer speeds `v_1..v_n` and threshold `1/(n+1)`, write

```text
Safe = { t in [0,1) : ||v_i t|| >= 1/(n+1) for every i }.
```

`Safe` is defined by **closed** conditions, so it may be a nonempty finite set of
points and still have `mu(Safe) = 0`. **Measure zero therefore means TIGHT, not
covered.** Concretely, every measure-zero instance in the census below has a
nonempty `Safe`, so LRC holds there -- with equality.

This is the repo's recorded **MAX-not-MEAN** failure mode ("existence is a
supremum") in its most direct form, and it is easy to walk into when computing
measures rather than witnesses. Any covering census must certify `Safe = empty`,
not `mu(Safe) = 0`.

## 2. The census (R2)

All `n`-subsets of `[1,B]`, exact rational interval arithmetic, threshold `1/(n+1)`:

```text
n = 3, B = 14:   364 sets,  4 with mu = 0   -- all dilated APs
n = 4, B = 14:  1001 sets,  5 with mu = 0   -- 3 dilated APs, and
                                               (1,3,4,7),  (2,6,8,14) = 2*(1,3,4,7)
n = 5, B = 12:   792 sets,  3 with mu = 0   -- 2 dilated APs, and (1,3,4,5,9)
```

So **the extremal/tight class is strictly larger than the dilated APs**. An
inverse theorem stated as "the only extremal configurations are dilated
arithmetic progressions" is **false as phrased**; `(1,3,4,7)` and `(1,3,4,5,9)`
are exact counterexamples to that phrasing. (These are, as expected, the
classically known sporadic tight instances at small `n`; the point here is that
the repo's frontier map should name them, since they are precisely what an
inverse theorem must accommodate.)

## 3. The witnesses (R3)

For each tight instance, enumerating `Safe` exactly on the common-denominator
grid:

```text
(1,2,3)         thr 1/4 :  t = 1/4, 3/4
(1,2,3,4)       thr 1/5 :  t = 1/5, 2/5, 3/5, 4/5
(1,3,4,7)       thr 1/5 :  t = 1/5, 2/5, 3/5, 4/5
(1,2,3,4,5,6)   thr 1/7 :  t = a/7
(1,3,4,5,9)     thr 1/6 :  t = 1/6, 5/6
(1,2,3,4,5,6,7) thr 1/8 :  t = 1/8, 3/8, 5/8, 7/8
```

**The reduced witness denominator is exactly `n+1` in every case** -- not merely "some
`q`". This sharpens the repo's `t = 1/q` reduction: on the tight locus the
denominator is pinned, and the sporadic instance `(1,3,4,7)` has *the same*
witness set as the AP `(1,2,3,4)`, which is a strong hint about what the inverse
theorem is really classifying (witness sets, not speed sets).

## 4. What does and does not transfer from FC(2) (R4)

**Does not: the transcendence mechanism.** FC(2)'s engine (THM-3031, THM-3039) is
that a counterexample pins an exponential period to a **rational** value -- the
simplex volume `1/(n-1)!` -- and transcendence of that period contradicts it. The
LRC analogue of "the value is rational" is genuine and was verified exactly:
every endpoint of a danger arc is a rational with denominator dividing
`(n+1)v_i`, so

```text
mu(Safe)  in  (1/N) Z,     N = (n+1) * lcm(v_1..v_n),
```

confirmed integral on every sample. But the pinned value is **0**, and `0` has no
arithmetic rigidity -- there is nothing for a transcendence statement to
contradict. The associated gap principle "either `mu = 0` or `mu >= 1/N`" is also
far too weak to use: measured `mu/(1/N)` was `30, 58, 66, 11020` on the
non-covering samples, i.e. the quantum is orders of magnitude below the actual
measure. **This route should not be pursued.**

**Does: the architecture only.** FC(2)'s proof is (i) exclude non-rigid
configurations by an inverse/monodromy argument, (ii) reduce to a rigid family,
(iii) kill that family by a separate arithmetic argument. That is exactly LRC
Route B's shape. What FC contributes is therefore not a tool but a warning about
step (ii): **the rigid family must be specified correctly**, and (R2) shows the
repo's working description ("dilated APs") is not it.

## 5. Scope

The census is finite and exact (`n <= 5`, `B <= 14`); it proves nothing outside
that box, and in particular does not classify tight instances for larger `n`.
(R1) and (R4)'s quantisation are proofs. (R3) is an exact observation on six
instances, not a theorem for all `n`. Nothing here bears on the `r`-ladder
(THM-1123/1132/1134) or on `LRC(14)` itself, which concerns `n = 13`.

## 6. Corrected exact companion

Run

```text
python 04-computation/lrc_tight_instance_census_thm3043.py
python -O 04-computation/lrc_tight_instance_census_thm3043.py
```

Both modes byte-match the stored transcript.  The companion rebuilds every
finite census row, requires a nonempty exact safe-point set for every
zero-measure row, checks the six displayed witness sets and eight
quantisation samples, and uses explicit runtime exceptions rather than
truth-bearing Python assertions.
