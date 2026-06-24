# LRC14 Multiplicity-Moment Dual

Created: 2026-06-24T16:41:40Z
Author: codex-2026-06-24-S156
Hypothesis: HYP-2971

This is a non-apex proof attempt, complementary to HYP-2970's endpoint-credit
winding-cycle dual.  For a row `S`, define

```text
X_S(t) = #{v in S : ||v t|| < 1/14}.
```

Then strict lonely mass is exactly `Pr[X_S=0]`.  If a polynomial satisfies

```text
P(0)<0,   P(k)>=0 for k=1,...,13,
```

then `E[P(X_S)]<0` is a rigorous certificate that `M(S)>1/14`.

Script and result:

```text
04-computation/lrc14_multiplicity_moment_dual_codex_s156.py
05-knowledge/results/lrc14_multiplicity_moment_dual_codex_s156.out
```

Default audit:

```text
rows = 17104
zero_safe_rows = 2
positive_safe_rows = 17102
certified_by_barrier_family = 17102
```

Only AP and GW are zero-safe.  Every other sampled qdiv-hard AP-neighborhood
row gets a negative admissible moment barrier from odd binomial barriers or
root barriers of degree at most 7, with a small `B13` tail.

The proof target is now:

```text
counterexample
  => moment-indistinguishable from AP/GW
  => satisfies all admissible multiplicity moment inequalities
  => must expose a genuine NORK/F7 packet or state-lift.
```

Tournament note: vertices are proof carriers, not runners.  The transitive
order is integer multiplicity histogram, dual moment barrier, AP/GW zero-safe
atom, NORK pinch template, qdiv witness, apex aperture, raw safe mass, raw
runner set.  This quotient intentionally forgets endpoint location so it can
serve as a global dual obstruction before labelled packet machinery starts.
