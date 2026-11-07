# Session222: L-shell versus two hypotenuse diagonals

## Useful result first

**PROVED, elementary geometry:** adding a one-cell L on both legs of a staircase is exactly equinumerous with adding two consecutive hypotenuse-parallel diagonals. Both turn the n-vertex tournament tile carrier into the (n+2)-vertex carrier and add **2n−1 cells**. An explicit cell permutation identifies the two constructions.

**PROVED distinction:** that permutation is not a tournament equivalence and does not preserve Hamiltonian path count H. Mode B's zero L-boundary has H unchanged; keeping the old cells anchored and adding two zero diagonals need not. The smallest nontrivial old carrier, n=3 with its sole bit 1, gives H=3 versus H=9 at order 5.

## Inheritance and provenance

Source revision: `aa977c095d7651e1f27c3b626fa61bed75cf4ae9`.

- Closest proved mechanism: `01-canon/theorems/THM-291-mode-b-multilinear-recursion.md`: source/sink extension, H_(n+2)(inner,0)=H_n(inner).
- Exact coordinate definition: `01-canon/theorems/THM-280-grid-reflection-is-complement.md`, Setup. Cells (x,y), 1≤y<x≤n, x−y≥2; bit zero means x→y, bit one y→x; consecutive descending arcs fixed.
- Corrected near miss: `01-canon/MISTAKES.md`, line 1848, Oct 4 tournament-recursion audit. Each boundary arm has n−1 cells, not n−2. It explicitly warns that Mode B is vertex order +2, not H+2.
- Least-used relevant sidecar: `01-canon/theorems/THM-793-mode-b-line-tower.md`, section (1): restriction to inner labels is induced-subtournament descent, commutes with bit complement and grid reflection, with free fiber 2^(2n−1) for n→n+2. Do not inherit its older unrelated OEIS labels or self-line conjecture: the same document records later corrections.
- Hostile mechanism: `01-canon/theorems/THM-284-linear-coefficient-formula.md` (title: Linear Coefficient of H in Tiling Bits = 2^(skip−1)). Single reversed arc has H=1+2^(x−y−1).
- Related genuine doubling: `01-canon/theorems/THM-452-doubling-tiling-coset-law.md`. Vertex doubling n→2n is a different operator; its canonical tiling is not generally blue, contrary to the original conjecture.
- Terminology hazard: THM-477 calls the reflection-fixed line x+y=n+1 the “hypotenuse”; in the present coordinate model it is the symmetry median a=b. The actual staircase hypotenuse runs parallel to x−y=constant. THM-291 uses hypotenuse for the fixed-base-path boundary. Always give equations rather than relying on that word.

### Historical recovery, pinned to the same commit

The local snapshot omitted two cited reflections, but both were retrieved directly from GitHub at `aa977c095d7651e1f27c3b626fa61bed75cf4ae9` and inspected in full:

- [Everything Is the Triangle](https://github.com/eliottcassidy2000/math/blob/aa977c095d7651e1f27c3b626fa61bed75cf4ae9/07-reflections/everything-is-the-triangle.md), kind-pasteur-2026-03-22-S20cm, blob `ca847958f3c63ccc8e35bc3ec180718aa87fbe82`. “Side 1” explicitly defines the staircase hypotenuse as cells `(i,j)` with `i+j=n−3`; “Side 2” and “Side 3” define the coordinate-axis legs. This is the historical geometric model recovered here. Our `(a,b)` are those nonnegative staircase coordinates, with a possible exchange of axes and reversed vertex-label convention. The displayed `H=1+2^d` is safely inherited only through THM-284's proved single-defect formula, with `d=skip−1`; its broad geometric rhetoric and asymptotic/graph claims are not imported. In particular, “the meta-graph is the quotient ... by the symmetric group” is not used as a group-action theorem on the fixed-path bit cube.
- [Unlocking G_n/Z_2 at All n](https://github.com/eliottcassidy2000/math/blob/aa977c095d7651e1f27c3b626fa61bed75cf4ae9/07-reflections/unlocking-gn-at-all-n.md), kind-pasteur-2026-03-23-S20cs, blob `6820ee57c09f570e448438e726b5326c1279c860`. This contains a text/ASCII dependency flowchart, not an L-versus-diagonal cell diagram. Its master-synthesis claims are historical, not mathematical dependencies of this report.

Neither fetched reflection contains a linked image artifact. The accessible geometric ASCII diagram in THM-291 was inspected alongside its corrected cell ranges. `shells.svg` is a new, exact coordinate illustration, not a reproduction of an older image. The Oct 4 results file cited by MISTAKES was not needed for the proof; the correction entry and repaired THM-291 supply the inherited statement.

Concept board: staircase geometry / induced subtournament / Hamiltonian path H / framed bit cube / dyadic scaling. Geometry preserves cells and bit freedom, not graph incidence; induced-subtournament bookkeeping is the missing coordinate.

## Exact shell theorem

Put Q_m={(a,b)∈Z_≥0²:a+b≤m}, m=n−3≥0. Map tournament tiles by

    (x,y) ↦ (a,b)=(y−1,n−x).

Thus |Q_m|=(m+1)(m+2)/2=binom(n−1,2). Grid reflection becomes (a,b)↦(b,a). The fixed base-path arcs lie just outside the free cells on a+b=n−2.

In Q_(m+2), there are two embeddings of Q_m:

1. **Leg growth:** E_L(a,b)=(a+1,b+1). Its complement L consists of a=0 or b=0. Each arm including corner has m+3=n cells; their overlap is one corner, hence |L|=2n−1.
2. **Diagonal growth:** E_D(a,b)=(a,b). Its complement D comprises a+b=m+1 and a+b=m+2. Their lengths are m+2=n−1 and m+3=n, so |D|=2n−1.

A concrete bijection F:Q_(m+2)→Q_(m+2) is

    F(a,b)=(a−1,b−1)                 if a,b≥1;
    F(0,b)=(b,m+2−b)                 if 0≤b≤m+2;
    F(a,0)=(a−1,m+2−a)               if 1≤a≤m+2.

The three disjoint domain pieces map bijectively onto Q_m, diagonal m+2, and diagonal m+1. This proves bijectivity, identifies the embeddings F∘E_L=E_D, and sends L onto D. It is a cell rearrangement, not a rigid planar motion; adjacency need not survive. This chosen shell bijection need not commute with reflection either.

In outer tournament coordinates the leg insertion sends (x,y)→(x+1,y+1). The original vertices become {2,…,n+1}, so all arcs and consecutive-path incidence restrict correctly. Conversely, diagonal anchoring sends an old tile (x,y)→(x+2,y): skip grows by two. This cannot arise from one common relabeling of old vertices in general, and no Hamiltonian-path preservation follows.

## H: proof, hostile, strongest survivor

For the zero L-boundary, new vertex n+2 is a source and 1 a sink. Every directed Hamiltonian path must start and end there. Removing them is a bijection to Hamiltonian paths in the old tournament. Therefore H_(n+2)(E_L(bits),0)=H_n(bits), for every n≥3 and every bit assignment. The difference polynomial contains at least one boundary variable in every monomial.

For a single bit at skip s=x−y≥2, the old descending path persists. Every additional path uses y→x. Each of the s−1 intermediate labels goes before the jump or after it, independently, giving 2^(s−1) additional paths. Exterior labels are forced before or after the interval. Thus H=1+2^(s−1).

Under zero diagonal growth by r layers, skip becomes s+r. On this **single-defect family only**,

    H_new−1 = 2^r (H_old−1).

One diagonal doubles H−1; two diagonals quadruple H−1. This is a genuine dyadic observable but not H→2H, nor a theorem for arbitrary multibit tilings. For n=3, bit (3,1)=1 has H=3. Leg insertion makes bit (4,2)=1 in order 5 and H=3. Diagonal anchoring makes bit (5,1)=1 in order 5 and H=9. Hence even the smallest informative example separates geometric equality from tournament/H equality.

Source: cell sets and bits. Target: same-sized cell sets and bits. Map: F. Preserved: cell cardinality, bit count, independence of coordinates, full Boolean-cube/Hamming structure. Destroyed: tournament edge incidence, generally spatial adjacency, and H. Sidecar needed: original pair labels plus induced-vertex embedding. Cheapest decisive test: the order-3 single-one carrier above.

## Four quantities that must stay separate

- Vertex count n: Mode B advances by 2.
- Free-cell count M_n=binom(n−1,2): advances by 2n−1.
- Number of framed tilings 2^M_n: multiplied by 2^(2n−1).
- H of one tournament: preserved on zero L-boundary, not determined by the geometric shell alone.

A literal factor-two recursion is impossible for H of nonempty tournaments if one invokes Rédei parity (H always odd). This report does not need that theorem for its claims; the single-defect formula already provides the honest factor-two law for H−1. Actual H+2 constructions require a separate graph operation and endpoint/path-count proof, not tile-cardinality arithmetic.

## Reproduction and scope

Run `python3 verify.py > verification.json` in this directory. `python3 -O verify.py` performs the same checks: explicit `require` calls remain active under optimization. Dependency-free Python uses exhaustive permutations rather than the proposed recurrences to count H.

FINITE-EXACT checks: all shell sets and bijections for n=3,…,100; all 2+8+64=74 framed inner tilings at n=3,4,5 for zero-boundary Mode B; every single-defect tile for n=3,…,7 (35 examples); the H=3 versus 9 hostile. These checks supplement, not replace, the all-order proofs. `shells.svg` illustrates n=5→7 (old 6 cells, new 15 cells, shell 9 cells).

## Additional hostile probe: strict H growth under diagonal anchoring

**FINITE-EXACT, not an all-order theorem.** `probe_monotonicity.py` applies the zero two-diagonal embedding to every one of the 32,768 framed tilings at order 7. Of these, 32,767 strictly increase H; the all-zero tiling alone keeps H=1. The smallest positive difference is 6. A separate deterministic sample of 3,000 nonzero masks at order 8 (with replacement, seed 2222026) also strictly increases H in every sample; the added zero control has equality. The order-8 sample is not exhaustive and its smallest observed gap is not a global bound.

The subset dynamic program counts paths by visited mask and final vertex. Its incoming adjacency matrix starts with descending arcs, flips exactly the specified nonconsecutive pairs, and sums predecessor endpoints over the smaller mask. Before this probe, it was independently checked against brute permutation counting on every framed tiling of orders 3, 4, and 5. `monotonicity.json` records the finite universes and results; reproduce with `python3 probe_monotonicity.py > monotonicity.json` (about 30 seconds here).

This was a targeted check of an intriguing signal supplied by the independent auditor at smaller orders. It does not support claiming that an all-order injection has been constructed. The first missing implication is precisely a map on arbitrary Hamiltonian paths: shifting only the high endpoint of each reversed arc does not map vertices consistently, and can change several path transitions at once. Stopping reason: the bounded hostile search was completed with a positive signal, but no incidence-respecting injection or general recurrence was established. Retain strict growth as a finite observation and potential future proof target, while the proved single-defect law remains the strongest all-order result here.
