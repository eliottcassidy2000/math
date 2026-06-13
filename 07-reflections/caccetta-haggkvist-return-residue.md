# Caccetta-Haggkvist as Return Residue

**Session:** codex-2026-05-30 S357
**Mode:** long creative pass with one cyclic Cayley probe
**Primary web anchors:** AIM workshop statement, Nathanson/Hamidoune additive
Cayley proof, nonuniform/rainbow CH, triangle flag-algebra progress.

The Caccetta-Haggkvist conjecture says:

```text
an n-vertex digraph with minimum outdegree at least n/k
has a directed cycle of length at most k.
```

The repo's residue calculus suggests a reframing:

```text
CH is not first about counting cycles.
It is about how long a directed expansion can avoid returning.
```

The residue is the first return boundary.

## Sources Read

The AIM workshop page gives the compact statement and emphasizes the network
of variants and additive-number-theory connections.  Nathanson's expository
paper records Hamidoune's proof for Cayley and vertex-transitive digraphs, the
part that matters most for this session.  Aharoni-Berger-Chudnovsky-Guo-Zerbib
prove the nonuniform bound

```text
directed_girth(D) <= 2 * sum_v 1/(outdeg(v)+1)
```

for digraphs without sinks, and formulate rainbow analogues.  Hladky-Kral-
Norin's triangle-free flag-algebra paper shows how far density/counted-flag
technology has gone in the `k=3` case.

These sources point in the same direction: CH has several faces, but the
strongest proven face is additive growth.

## Feedback Loop

### Iteration 1: Translate CH Into Return Time

For a digraph `D`, define:

```text
g(D) = directed girth
d    = minimum outdegree
```

CH predicts:

```text
g(D) <= ceil(n/d).
```

If `g(D) > k`, then for every vertex `v`, the first `k` layers of directed
walks avoid returning to `v`.  This makes every rooted expansion tree a
return-free transport system.  The conjecture says that with row budget `d`,
such transport cannot stay return-free longer than `n/d` steps.

That is already residue language:

```text
out-neighborhood growth = projected shadow
back edge / closed walk  = return residue
no short cycle           = vanishing return residue through time k
```

### Iteration 2: Use Cayley Graphs As The Clean Laboratory

For `Cay(G,A)`, every vertex has the same out-neighborhood.  A closed directed
walk of length `ell` is a word

```text
a_1 + ... + a_ell = 0,  a_i in A.
```

Minimal such `ell` gives a simple directed cycle, because any repeated
interior vertex would split off a shorter zero-sum subword.

So the Cayley version is:

```text
if |A| >= |G|/k, then 0 lies in ell*A for some ell <= k.
```

Hamidoune's proof, as presented by Nathanson, uses additive growth of
`B=A union {0}`.  The core obstruction to a counterexample is:

```text
no zero-sum return  ->  repeated sumsets of B must keep growing.
```

This is the same pattern as Lonely Runner S356:

```text
forbidden return fibers occupy a quotient;
witness/cycle appears when a boundary residue cannot stay empty.
```

### Iteration 3: Run The Cyclic Probe

`04-computation/caccetta_haggkvist_residue_probe_s357.py` enumerates all
cyclic connection sets

```text
A subset {1,...,n-1}, |A| <= floor(n/3), n <= 24.
```

For each `A`, it computes:

```text
first_zero_sum_length(A)
CH_bound = ceil(n/|A|)
growth_j = |j(A union {0})|
Kemperman_slack_j = growth_j - (j(|A|+1)-j+1)
pre_return_residue = n - |(CH_bound-1)(A union {0})|.
```

The sample census:

```text
total_sets = 1,612,949
CH violations = 0
CH-tight sets = 62,966
non-interval tight sets = 62,064
tight with zero pre-return residue = 18,900
tight with zero Kemperman slack until return = 1,114
```

The canonical tight family is the cyclic interval construction:

```text
n = d(g-1)+1
A = {1,...,d}
first_zero_sum_length = g
|j(A union {0})| = jd + 1 for j=1,...,g-1.
```

This is almost too clean: the pre-return layers are exactly arithmetic
progressions.  The first return occurs precisely when the progression must
wrap around the circle.

## What The Probe Says

The naive hope was:

```text
CH-tight Cayley sets are just intervals up to units.
```

False.  There are many non-interval tight sets even for small `n`.

The better hope is:

```text
Every hard case has a quotient where the growth profile is progression-like.
```

Non-interval tight sets can have positive Kemperman slack, zero pre-return
residue, or both.  That suggests three regimes:

1. **Pure progression:** zero slack until return.  The interval family lives
   here.
2. **Early expansion / late return:** slack is positive, but the zero-sum
   boundary still waits until the CH limit.
3. **Covered-before-return:** sumsets cover the whole group before the first
   zero-sum word using nonzero letters.  This is a boundary distinction between
   `B=A union {0}` and `A`.

The third regime is the interesting residue.  `B` can cover zero trivially
using all-zero padding, while `A` still has no nontrivial return.  This is the
same boundary/tight-witness distinction that appeared for Lonely Runner.

## General Digraph Translation

For arbitrary digraphs there is no group operation.  But a minimal CH
counterexample would still have local data:

```text
L_j(v) = vertices reachable from v in <= j steps
R_j(v) = return walks from v to v of length j
C_j(v) = edges from layer j back into earlier layers
```

No short directed cycle means:

```text
R_j(v) = 0 for j <= k.
```

The conjecture says this cannot persist with minimum outdegree `d >= n/k`.
So the task is to build a fake additive quotient out of the layer data.  The
candidate residue vector is:

```text
CH_res(v) = (
  |L_1(v)|, ..., |L_k(v)|,
  second-neighborhood surplus,
  backward-collision profile,
  return-walk profile,
  harmonic weight 1/(outdeg(v)+1)
).
```

This makes the Aharoni-Berger-Chudnovsky-Guo-Zerbib harmonic bound feel less
like a separate theorem and more like a weighted residue checksum.

## Connections Back To The Repo

### Deletion-Residue Rank

HYP-1785 says exact kills and dangerous near-kills should be separated by the
rank of what survives after deletion.  For CH, the analogous operation is not
deleting a vertex but deleting the first return:

```text
return-free layers = exact kill of closed walks
first return       = residue boundary
```

The triangle case should have an especially cheap residue:

```text
N+(v), N++(v), and edges N++(v) -> v.
```

This is close to the second-neighborhood conjecture and to the flag-algebra
triangle-free work.

### Quotient Gaps

HYP-1783 says empty quotient fibers are zero transport rows/columns.  A CH
counterexample would be a tower of empty return fibers:

```text
return length 1 empty,
return length 2 empty,
...
return length k empty.
```

The proof problem is to show the transport matrix cannot keep all those
columns empty while every row has budget `d`.

### Endpoint/Rainbow Transfer

Rainbow CH replaces vertex out-stars by color classes.  That is almost exactly
the endpoint-transfer viewpoint: each row/color/fiber wants to contribute one
edge to a cycle.  A rainbow cycle is a return residue using distinct source
fibers.

The recent endpoint collision hypergraph work in HYP-1790--HYP-1793 suggests
a possible engineering analogy:

```text
non-private collision columns in endpoint transfer
support-3 rainbow triangle obstructions
leaf-peelability as a rank proof
```

This is speculative, but it is the sort of speculation with a computation
attached.

## New Thesis

The Caccetta-Haggkvist conjecture is a return-residue theorem:

```text
High outdegree gives transport budget.
No short directed cycle kills return residues.
Growth must then be progression-like.
Progression-like growth cannot avoid wrapping/returning past n/d.
```

For Cayley digraphs this is additive number theory.  For general digraphs the
missing object is a quotient that makes rooted layer growth behave like a
sumset.  The creative task is not to count all short cycles; it is to find the
right layer quotient where the absence of short cycles becomes impossible
Kemperman equality.

## Next Work

1. Classify the `62,966` cyclic tight sets by unit action and subgroup
   quotient.
2. Add a fast profile extractor for arbitrary digraphs:

```text
rooted layer sizes,
return-walk counts,
second-neighborhood surplus,
collision edges into previous layers.
```

3. Generate small triangle-free oriented graphs with large minimum outdegree
   and rank them by CH residue.
4. Build a rainbow-cycle toy extractor and compare its collision hypergraph to
   the endpoint-transfer collision hypergraphs.
5. Try to formulate a local lemma:

```text
If every rooted layer profile has Kemperman-equality growth up to time k,
then the digraph admits a cyclic interval quotient.
```

That lemma would not solve CH.  But it would tell us what a minimal
counterexample is pretending to be.
