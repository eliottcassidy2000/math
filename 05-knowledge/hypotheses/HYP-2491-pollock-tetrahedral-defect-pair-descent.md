# HYP-2491 - Pollock tetrahedral five-sum reduces to four-defect triangular-spacing descent

**Status:** OPEN proof program with exact finite scout through `10^6` and shell-stencil audit through `k=1200`.
**Source:** codex-2026-06-13.
**Extends:** HYP-1962, HYP-1963, HYP-2453, HYP-2487..HYP-2490, THM-498, T797.
**Artifacts:** `04-computation/pollock_tetrahedral_defect_descent_codex.py`, `05-knowledge/results/pollock_tetrahedral_defect_descent_codex.out`, `07-reflections/pollock-tetrahedral-defect-pair-descent.md`.
**External anchors:** Pollock tetrahedral conjecture remains open; MathWorld/OEIS report `241` known four-tetrahedra defects, apparently ending at `343867`: <https://mathworld.wolfram.com/PollocksConjecture.html>. Brady et al. survey the generalized tetrahedral variant and note the same frontier: <https://www.tandfonline.com/doi/full/10.1080/00029890.2021.1982635>.

**Concurrent context:** kind-pasteur-2026-06-13-S3 claimed HYP-2487..HYP-2490 for the Pollock/Linnik/circle-method and cycle-spectrum hierarchy.  HYP-2491 is the complementary defect-pair route: it keeps the modern analytic template as the likely tail engine, but isolates a weaker no-long-triangular-self-correlation lemma that would already imply the five-tetrahedra result once paired with a finite stencil certificate.

## Statement

Let

```text
Te_k = k(k+1)(k+2)/6,  k >= 1,
R_4 = {sums of at most four Te_k},
D_4 = positive integers not in R_4.
```

Pollock's tetrahedral conjecture follows from either of two defect-set
statements.

### Strong Tail

```text
D_4 is finite and max(D_4)=343867.
```

This is the Salzer-Levine/OEIS-strength frontier.  The computation rediscovers
exactly `241` defects through `10^6`, with largest `343867`, and verifies that
every integer through `10^6` is a sum of at most five tetrahedral numbers.

### Weaker One-Back Descent

For all sufficiently large shell indices `k`, there is no pair

```text
r, r + tri(k) in D_4
with 0 <= r < Te_{k+1} - Te_k,
tri(k)=Te_k-Te_{k-1}=k(k+1)/2.
```

Indeed, if `n=Te_k+r` lies in the `k`th tetrahedral shell and `r in R_4`, then
`n=Te_k+r` is already a five-term representation.  If `r notin R_4` but
`r+tri(k) in R_4`, then

```text
n = Te_{k-1} + (r + tri(k))
```

is a five-term representation.  Thus a one-back failure is exactly a
triangular-spacing pair inside `D_4`.

## Computed Frontier

The scout verifies:

```text
missing from <=5 tetrahedra through 1000000: 0
min-term histogram: {1: 180, 2: 14152, 3: 444283, 4: 541144, 5: 241}
four-term defect count: 241
largest four-term defect: 343867
```

The last four-term defects through `10^6` are:

```text
44833, 47627, 48043, 56467, 56842, 58613, 59077, 62158,
64752, 65253, 65567, 71157, 74687, 78003, 78787, 83603,
84023, 85993, 91128, 106277, 113062, 134038, 148437, 343867.
```

Among these `241` defects, there are `601` pairs separated by a triangular
number.  The largest triangular separation is:

```text
3142 -> 343867
gap = 340725 = tri(825).
```

Thus, relative to the known defect frontier, the last one-back obstruction
occurs at shell `k=825`.

## Shell-Stencil Audit

For shells

```text
[Te_k, Te_{k+1}),  1 <= k <= 1200,
```

define an anchor-offset stencil: offset `j` covers `n` if
`n-Te_{k-j} in R_4`.

Exact bitset slicing gives:

```text
minimal needed offset width histogram: {0: 4, 1: 1035, 2: 156, 3: 5}
after offsets <=0: failures=1196, last=1200
after offsets <=1: failures=161, last=825
after offsets <=2: failures=5,   last=121
after offsets <=3: failures=0
```

The only shells needing the third fallback are:

```text
15, 24, 56, 89, 121.
```

This suggests a practical proof decomposition:

1. Certify finitely many shells up to `825` by a width-3 stencil.
2. Prove the tail one-back triangular-spacing lemma for `k>825`.

## Tournament Analysis

The offset tournament used vertices

```text
j = 0,1,2,3, meaning anchor Te_{k-j}.
```

Pairwise observable:

```text
private shell positions covered by offset i and not by offset j,
summed over shells k <= 1200.
```

Gauge:

```text
i -> j iff i has larger private coverage.
```

Tie Hamiltonian path:

```text
0 -> 1 -> 2 -> 3.
```

Result:

```text
covered counts: {0: 289186020, 1: 289432265, 2: 289437717, 3: 289439464}
edges: 1->0, 2->0, 3->0, 2->1, 3->1, 3->2
score histogram: {0:1, 1:1, 2:1, 3:1}
directed 3-cycles: 0
SCCs: [(3,), (2,), (1,), (0,)]
Hamiltonian paths: [(3,2,1,0)]
```

The older anchor dominates because it asks the four-sum set at a larger
remainder, where `R_4` is denser.  This mirrors the LRC blocking-height theme:
dominance grows in raw cumulative load, but the proof-relevant object is the
tail defect/correlation, not the scalar dominance ranking alone.

## Octahedral Sibling

The same script includes a small Pollock octahedral scout.  Through `100000`,
every integer is a sum of at most seven octahedral numbers, with exactly `12`
sampled cases needing seven:

```text
17, 35, 55, 61, 73, 79, 80, 200, 206, 213, 225, 309.
```

This is only a side diagnostic, but it suggests reusing the defect-spacing
method for octahedral numbers, where the known analytic result gives all
sufficiently large integers and leaves a finite computation.

## Proof Program

The high-leverage next lemma is:

```text
No long triangular self-correlation of D_4:
for k > 825, D_4 cap (D_4 - tri(k)) has no element r < tri(k+1).
```

Possible routes:

1. **Density tail:** prove that above a threshold every interval
   `[x, x+tri(k)]` of the relevant shape hits `R_4` at one endpoint.
2. **Modular lift:** show that a hypothetical pair of four-sum defects with
   triangular separation forces incompatible local representations modulo a
   selected family of moduli.
3. **Convolution lift:** encode `R_4` as a four-atom convolution/tiling image;
   a defect pair is then two adjacent missing boundary totals separated by the
   tetrahedral shell gap.
4. **Circle-method finite tail:** adapt the octahedral/cubic Waring strategy:
   prove all sufficiently large numbers are in `R_4`, then close the remaining
   finite interval with the bitset certificate.

The key shift is that Pollock need not be attacked as an arbitrary five-sum
problem.  It can be attacked as a sparse forbidden-spacing theorem for the
four-sum defect set.
