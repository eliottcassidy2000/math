---
id: HYP-1803
status: EXPLORATORY
source: codex-2026-05-30 S357
related:
  - HYP-1781
  - HYP-1783
  - HYP-1785
  - HYP-1794
---

# HYP-1803: Caccetta-Haggkvist Criticality Is Return Residue

## Statement

Near-extremal Caccetta-Haggkvist examples should be controlled by a delayed
return residue:

```text
first_return_time  = directed girth
growth_profile     = sizes of out-neighborhood layers or sumsets
return_residue     = the first boundary where expansion can no longer avoid
                     returning to the start
```

In a Cayley digraph `Cay(G,A)`, a directed cycle of length `ell` is a zero-sum
word of length `ell` in `A`.  The conjecture becomes:

```text
if |A| >= |G|/k, then 0 in ell*A for some ell <= k.
```

The hypothesis is that any genuinely hard non-Cayley instance must mimic the
same additive equality pattern locally: out-neighborhood layers grow at nearly
the minimum possible rate until a forced return boundary appears.

## Evidence

`04-computation/caccetta_haggkvist_residue_probe_s357.py` exhaustively checks
cyclic Cayley connection sets

```text
A subset {1,...,n-1}, |A| <= floor(n/3), n <= 24.
```

The probe found:

- `1,612,949` connection sets checked;
- `0` violations of `first_zero_sum_length <= ceil(n/|A|)`;
- `62,966` CH-tight sets with `first_zero_sum_length = ceil(n/|A|)`;
- `1,114` tight sets whose Kemperman slack is zero all the way until the
  final return step;
- the canonical interval family `A={1,...,d}` in `n=d(g-1)+1` has exact
  growth

```text
|j(A union {0})| = jd + 1,  j=1,...,g-1,
```

and first zero-sum length `g`.

The interval family is the clean model: no early zero-sum, minimum additive
growth, then return exactly at the CH boundary.  Non-interval tight sets exist
abundantly, but their profiles still show the same residue vocabulary:
pre-return coverage, Kemperman slack, and first return time.

## Predictions

1. Cayley-tight sets should decompose into low-dimensional progression or
   subgroup-quotient mechanisms after quotienting by obvious symmetries.
2. General minimal CH counterexamples, if they exist, should have very small
   variance in rooted reachability-layer growth.
3. Triangle-free high-outdegree digraphs should concentrate around minimal
   second-neighborhood surplus; flag-algebra certificates are detecting
   violated residue inequalities in that surplus ledger.
4. Rainbow versions of CH should be expressible as transversal return residues:
   each color class is an out-star substitute, and a rainbow cycle is a
   return using distinct source fibers.

## Test Plan

1. Extend the cyclic probe with canonical quotient classification of tight
   sets under units and subgroup factors.
2. Add a general digraph sampler that records, for each vertex `v`,

```text
L_j(v) = vertices reached from v in <= j steps,
R_j(v) = number of return walks from v to v of length j,
collision_j(v) = edges from layer j back into earlier layers.
```

3. Search small triangle-free oriented graphs with large minimum outdegree and
   rank them by second-neighborhood surplus.
4. Compare those profiles to Cayley interval profiles and to the harmonic
   nonuniform-degree bound of Aharoni-Berger-Chudnovsky-Guo-Zerbib.

## Sources

- `04-computation/caccetta_haggkvist_residue_probe_s357.py`
- `05-knowledge/results/caccetta_haggkvist_residue_probe_s357.out`
- `07-reflections/caccetta-haggkvist-return-residue.md`
- Nathanson, "The Caccetta-Haggkvist conjecture and additive number theory",
  arXiv:math/0603469.
- Aharoni, Berger, Chudnovsky, Guo, Zerbib, "Non-uniform degrees and rainbow
  versions of the Caccetta-Haggkvist conjecture", arXiv:2110.11183.
