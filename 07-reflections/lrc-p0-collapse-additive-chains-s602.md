---
source: codex-2026-06-03-S602
status: exact bounded audit plus new classification subproblem
tags: [LRC, p0-collapse, additive-chain, unit-boundary, transversals, nonunit-shells]
---

# LRC `p0` Collapse And Additive Chains

The user's note changes the shape of the tight-family problem:

```text
the p_0=0 collapse family is larger than AP.
```

That is right in the endpoint-collapse sense already present in S357/S359.  If
`p_0=0` means "no open safe interval, but boundary witnesses survive and the
visible witness quotient collapses to `n`," then the family contains the AP
rows and also sparse additive chains.

Incoming integration during close-out: opus-S599 made `p_0` literal as the
zero-depth cell in the covering-depth distribution

```text
p_k(delta) = meas{t : exactly k runners are within delta of the origin}.
```

This S602 note is the endpoint-critical refinement of that same object.  At the
reduced threshold, `p_0=0` means the open zero-depth cell has vanished; the
extra S602 condition is that isolated boundary witnesses remain and collapse
to the unit quotient `Z/nZ`.

The two named rows are not random exceptions:

```text
(1,3,4,7)      = seeds 1,3; then 4=1+3; top 7=3+4
(1,3,4,5,9)    = seeds 1,3; then 4=1+3; 5=1+4; top 9=4+5
```

Both have the unit boundary skeleton, both are perfect transversals modulo
`C=2n-1`, and both have flip-set `{2}`.  In the S553 language, they are exactly
the primitive one-flip floor transversals for `n=5` and `n=6`.

## What S602 Checked

I wrote `04-computation/lrc_p0_collapse_additive_chains_s602.py` and stored
the exact output in `05-knowledge/results/lrc_p0_collapse_additive_chains_s602.out`.
The predicate is deliberately operational:

```text
max open safe gap = 0,
boundary witnesses exist,
witness denominator lcm = n,
boundary witnesses = {a/n : gcd(a,n)=1}.
```

The targeted boxes recover nine primitive collapsed rows:

```text
n=4: AP
n=5: AP, (1,3,4,7)
n=6: AP, (1,3,4,5,9)
n=7: AP
n=8: AP, (1,2,3,4,5,7,12), (1,4,5,6,7,11,13)
```

All nine are two-seed addition chains.  AP is therefore not wrong; it is just
the all-lower, most regular member of a larger additive-chain boundary family.

Read this together with S599, not as a competing definition.  S599 proves the
first moment `E[depth]=2m delta` is configuration-independent, checks the named
chains as `p_0=0`, gives the non-chain control `(1,2,4,7)` with `p_0=6/35`,
and fits the near-collapse exponent `p_0(delta_c-epsilon) ~ epsilon^1` for AP
and `(1,3,4,7)`.  S602 adds the unit-boundary quotient, flip-set, and nonunit
shell labels that the raw depth distribution forgets.

## The New Split

There are two visibly different non-AP branches.

**Prime/transversal branch.**  For `n=5,6`, the sporadics are perfect
antipodal transversals modulo `C=2n-1` with flip-set `{2}`.  They still occupy
one residue from each shell, so the unit-inverse witness machinery sees them
the same way it sees AP, except the additive chain changes which top cap closes
at the boundary.

**Composite/nonunit branch.**  For `n=8`, `C=15` is composite.  The two non-AP
rows are not perfect transversals:

```text
(1,2,3,4,5,7,12)       misses nonunit shell 6
(1,4,5,6,7,11,13)      misses nonunit shell 3
```

They still collapse to the unit boundary skeleton, but the `C`-support story
has passed through nonunit shell ramification.  This is the small model for
the n=14 problem, where `C=27=3^3`.

## Why The Additive Chain Matters

The tempting old lens was "small sumset" or "AP-like."  S572 already weakened
that: `(1,3,4,5,9)` has positive sumset excess but stays at the floor.  S602
sharpens the replacement lens:

```text
small raw sumset       is too coarse;
two-seed generation    captures the collapsed rows in these boxes;
top-sum pincer         identifies the sparse non-AP examples;
C-shell labels         decide transversal versus nonunit branch.
```

The top equation is especially suggestive.  In both user examples, the top is
the sum of the two previous chain terms.  In one n=8 row this remains true:
`12=5+7`.  In the other, the top has a sum pair `13=6+7` but not the immediate
previous-two form.  So "top=sum of two below" is a real marker, while
"previous two" is too narrow for the composite branch.

## S603: Relation To The Master-Object Baseline

The S599b/HYP-2154 master-object pass gives this subproblem a cleaner
probabilistic shape.  The covering-depth distribution `{p_k}` has a free
baseline: if the danger arcs behaved independently, depth would look roughly
Poisson with mean `2n delta`, and `p_0` would remain positive near the LRC
threshold.  The additive chains are the opposite edge: structured arithmetic
correlation drives the bulk `p_0` cell to zero.

So the new classification problem is not "find all AP variants."  It is:

```text
which additive/shell correlations can empty p_0 without killing the unit
boundary witness floor?
```

The two named rows are the first clean anti-Poisson examples.  They are
correlated enough to close every open lonely interval, but not so pathological
that the `1/n` unit skeleton disappears.  That is exactly the large-deviation
tail that HYP-2154 says the LRC proof has to control.

## Tournament Analysis / Assumption Challenge

The script makes the quotient explicit.  It does not put runners or arcs at
the tournament vertices.  It uses proof lenses:

```text
unit boundary, two-seed chain, top-sum, C-transversal, AP, flipset {2},
nonunit C-hole, C-prime.
```

The pairwise observable is:

```text
(collapse hits, -false positives, lens name).
```

This produces a transitive tournament with no directed 3-cycles and singleton
SCCs.  The ranking is unsurprising but useful:

```text
unit_boundary_skeleton >
two_seed_addition_chain >
top_has_sum_pair >
perfect_C_transversal >
AP_only > ...
```

The assumption challenge is the main point.  If we force the vertex set to be
"runners" or "arcs," the additive-chain regularity is hard to see.  If we use
proof obligations as vertices, AP-only immediately loses recall.

## Next Tasks

1. Formalize the user's `p_0=0` notation and compare it with the operational
   endpoint-collapse predicate.
2. Prove the `{2}` one-flip classification in the prime/transversal branch, or
   find the first larger-`n` exception.
3. Classify composite nonunit additive chains by the missing shell and gcd
   stratum of `C=2n-1`.
4. Carry this to `n=14`, where the relevant composite shell modulus is
   `C=27=3^3`.
