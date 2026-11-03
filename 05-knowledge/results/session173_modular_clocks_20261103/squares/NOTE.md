# Square-root words, rigid placements, and optimal arithmetic connectors

**Status:** PROVED elementary placement and rigidity theorems, explicit optimal connector certificates, and FINITE-EXACT experiments. The known square-sum Hamiltonicity range is CITED, not improved. No claim of literature priority, uniform random-rotation termination, or a Collatz implication is made.

## 1. Constructive sources and inherited material

Robert Gerbicz's own entry in [OEIS A090461](https://oeis.org/A090461), dated January 21, 2018, states existence of paths for every n≥25 and cycles for every n≥32, with a deterministic construction linked at the forum. The linked forum page returned HTTP 403 during this work, so its proof was not inspected and is not re-proved here.

A primary implementation, [devisnik/square-paths](https://github.com/devisnik/square-paths), provides a 32-cycle and describes cycle extension by cutting and reinserting. Its [GerbiczExtender source](https://github.com/devisnik/square-paths/blob/master/src/main/java/extend/GerbiczExtender.kt) selects two neighbors of the incoming vertex and tests whether their successor ports can be rejoined. Its [CutAndInsertExtender source](https://github.com/devisnik/square-paths/blob/master/src/main/java/extend/CutAndInsertExtender.kt) keeps an exposed segment and searches compatible insertion ports. These are constructive prompts for retaining exact endpoints and occupied vertices, not just an unpointed existence count. We independently use ordinary endpoint rotations, not those implementations, below.

Repository inheritance: S172's `squares/NOTE.md` gives exact one-vertex Hamilton response via endpoint and two-path-cover counts. `crossroads_family_20260926_squares.md` already has reflection composition, the first square-sum 4-cycle at46, and designated-square-preserving affine lifts. Those mechanisms are not new here. The delta is the **full interval placement criterion, a sharp two-placement rigidity classification, and cardinality-optimal connector paths for the named anchors**. Only those targeted notes and THM-4505 were consulted; no forbidden canon file was read.

## 2. An exact integer lift of a square-root clock

For a positive integer root word w=(s₁,…,sₘ), let c₀=0 and cᵢ=sᵢ²−cᵢ₋₁. Every walk realizing that word is necessarily

    vᵢ(x) = (−1)^i x + cᵢ,             0≤i≤m.                    (1)

This is inherited reflection algebra. Its endpoint map is the affine isometry T_w(x)=ε_w x+d_w, with ε_w=(−1)^m, d_w=cₘ. Concatenation has the exact composition law

    (ε_v,d_v) ∘ (ε_w,d_w) = (ε_v ε_w, ε_v d_w+d_v).            (2)

Thus even blocks carry translations and odd blocks carry reflections. Reduction modulo8 and modulo9 is the reduction of this same action modulo72, by CRT. It forgets the integer multiple of72 in the translation, not a mysterious new arithmetic operation.

**Theorem 1 (all simple interval placements).** Put E={cᵢ:i even} and O={cᵢ:i odd}, retaining multiplicities for the initial test. If either list repeats an offset, there is no simple placement. Otherwise the exact legal integer starts for a simple path in [L,U] are

    [ max({L−e:e∈E} ∪ {o−U:o∈O}),
      min({U−e:e∈E} ∪ {o−L:o∈O}) ] ∩ Z

with every integer (o−e)/2 removed, for e∈E,o∈O of equal parity.

**Proof.** Equation(1) gives each interval bound separately. Same-parity vertices coincide exactly when their offsets coincide. Opposite-parity vertices coincide exactly when 2x=o−e. These conditions exhaust bounds and collisions. ∎

A path block therefore has an exact composable interface consisting of its affine endpoint map and its signed occupied offsets. On concatenation, replace each second-block affine expression a y+b by a(ε_w x+d_w)+b, remove the shared seam once, and test collisions by Theorem1. This is an explicit arithmetic specialization of graph-interface memory. It is sufficient, not claimed to be minimal, bounded-size, or efficient for Hamilton search.

### A real modular seam that fails over the integers

Our certified Hamilton path on [527] contains

    230 → 395 → 446,       square roots 25,29.

The integer endpoint displacement is 29²−25²=216=3·72. Both mod8 and mod9 clocks return to their start, but the exact ports differ. A modularly closed block is not an integer cycle. Retaining the integer translation closes this specific information gap.

At height44, mod72 is injective on vertex labels, so equal mod8/mod9 port labels **do** authenticate endpoint equality. That does not make modular squareness sufficient: 1+27=28 is a square residue modulo8 and modulo9, but is not an integer square. Root labels, integer lifts, and collision information are separate requirements.

## 3. Sharp rigidity: an entire exact root word can have two Hamilton placements

**Theorem 2 (at most two placements).** Fix an ordered root word of length n−1. It has at most one placement that is a permutation of [n] when n is odd, and at most two when n is even.

**Proof.** For odd n, summing(1) gives Σvᵢ=x+Σcᵢ, so the required sum n(n+1)/2 uniquely determines x. For even n, the first moment loses x, but

    Σvᵢ² = n x² + 2x Σ(−1)^i cᵢ + Σcᵢ².

Equating this with n(n+1)(2n+1)/6 is a nonzero quadratic, with at most two real roots and therefore at most two integer starts. ∎

The moment argument works for arbitrary prescribed edge sums, not only squares. Squareness supplies the concrete arithmetic witness below.

**Theorem 3 (exact classification of the double-placement case).** Suppose p=(p₀,…,pₙ₋₁) is a permutation of [n]. A second placement with the same edge sums and larger start p₀+t, t>0, exists if and only if

1. n is divisible by 2t;
2. the values in the even positions of p are exactly the alternating lower t-blocks

    {1,…,t} ∪ {2t+1,…,3t} ∪ … ∪ {n−2t+1,…,n−t}.

The second placement is qᵢ=pᵢ+(−1)^i t. Consequently, if the first path has square sums, so does the second, with the identical ordered roots.

**Proof.** Equality of edge sums gives qᵢ−pᵢ=(−1)^i t. Thus the map on the value set [n] sends each even-position value up t and each odd-position value down t. If q is a permutation, this is a permutation along the chains of integers with the same residue modulo t, moving every element exactly one step up or down. At the bottom of a chain the only possible move is up; bijectivity forces the second element down. Removing that pair and repeating shows that the entire chain is partitioned into adjacent transpositions. Every chain therefore has even length, and the up-moving elements are the first, third, fifth,… elements in it.

Write n=qt+r, 0≤r<t. The t chains have lengths q or q+1 (and t≤n is necessary). If r>0, both lengths occur and cannot both be even. Hence r=0 and q is even. Their alternating up-moving elements are precisely the lower t-blocks in the theorem. Conversely this block coloring makes the shifts adjacent transpositions in every residue chain, hence a permutation. ∎

This classification is stronger than the two-root moment bound and explains exactly which boundary information the full ordered root word can still lose.

### A sharp square-sum witness at n=36

The first placement is

    3,33,16,9,27,22,14,11,25,24,1,35,29,20,5,31,18,7,
    2,23,26,10,6,19,30,34,15,21,4,12,13,36,28,8,17,32.

Shifting alternate entries by +6,−6 gives

    9,27,22,3,33,16,20,5,31,18,7,29,35,14,11,25,24,1,
    8,17,32,4,12,13,36,28,21,15,10,6,19,30,34,2,23,26.

Both have root word

    6,7,5,6,7,6,5,6,7,5,6,8,7,5,6,7,5,3,
    5,7,6,4,5,7,8,7,6,5,4,5,7,8,6,5,7.

Theorem1 gives the start interval [3,9]. The numbers of distinct visited vertices for starts3,…,9 are respectively36,30,24,18,24,30,36. Exactly starts3 and9 survive. No minimality claim is made about n=36.

The twin paths do not have the same authenticated starting residue modulo8 or9. They demonstrate the need to retain a start port, even with the full exact edge word. The [527] subpath instead demonstrates a nonzero integer translation invisible to both clocks. These are different losses of information.

## 4. Exact optimal connection of 527,44 and the six-number motif

Let A={4,2,5,3,7,13,44,527}. The square-sum graph induced by A has exactly three components:

    44–5–4,       527–2–7,       3–13.

**Proposition (minimum connector count).** Exactly two additional positive integers suffice to place all anchors in a single simple square-sum path. One certificate is

    527–2–7–317–44–5–4–12–13–3.                              (3)

The roots on its edges are23,3,18,19,7,3,4,5,4. Thus the two extra vertices are12 and317, both below527.

**Proof of minimality.** Deleting k added vertices from any simple path leaves at most k+1 nonempty old path pieces. Each lies in one component of the old induced graph. Covering c old components requires at least c pieces, so k≥c−1. Here c=3, and(3) attains k=2. ∎

Including8 as another required anchor creates a fourth component, because8 has no square-sum neighbor in A. Exactly three additions are necessary and sufficient:

    527–2–7–317–44–5–4–12–13–3–1–8.                          (4)

The extra vertices are{1,12,317}, and the final roots are2,3. These certificates also minimize the maximum added vertex among paths using the minimum number of additions, as proved below.

The connectors come from actual difference-of-squares equations. To join7 and44 through z, solve (z+44)−(z+7)=37=19²−18², giving z=317. To join4 and13, use13−4=5²−4², giving z=12. To join3 and8, use8−3=3²−2², giving z=1. More generally a simple one-vertex connection between distinct a,b exists exactly when a positive integer z outside {a,b} satisfies z+a=u²,z+b=v². When splicing into an occupied path bank, z must additionally avoid every already occupied vertex. An appropriate same-parity factor pair of b−a supplies a candidate, which must pass those exclusions. For example a=2,b=7 gives the positive candidate z=2, but this repeats an endpoint and is not a simple connector. This is a constructive seam calculation, not a residue coincidence.

### Minimum connector height at the minimum count

In the equality case k=c−1 of the component lower bound, deleting added vertices gives exactly c nonempty pieces. Each old component is therefore traversed as one whole path, and each added vertex joins endpoints of two different old components. In particular the component527–2–7 must connect through one of its endpoints527 or7 to an endpoint of another component.

Factoring the corresponding square differences gives **all** positive candidates:

- from7 to3,4,13 (and optionally8): none;
- from7 to44:317;
- from527 to3:16897;
- from527 to4:68117;
- from527 to13: none;
- from527 to44:917,6197,58037;
- from527 to8, when8 is required:7217,67073.

Every valid minimal-count augmentation must therefore use a connector at least317. Paths(3) and(4) attain that lower bound. Thus317 is the exact least possible maximum added vertex, **conditional on using the minimum number of added vertices**. More connectors might lower the height; that optimization is not claimed here. The finite factor list is exhaustive because u²−v²=(u−v)(u+v), and both factors are positive and have the same parity. The checker independently tests every listed candidate and rejects occupied vertices.

## 5. Experiment and scope of further work

`check.py` starts from the cited32-cycle and uses deterministic seeded endpoint rotations and reversals to extend paths through every size33,…,527. Each resulting path is checked independently by its complete vertex set and exact integer square tests. Certificates retain full paths at32,44,72,527 and all extension rotation counts. This is a finite experiment; it is not Gerbicz's all-n proof and does not establish termination beyond the tested interval.

The checker also exhausts root words of length0,…,5 over roots2,…,7, comparing the placement theorem with direct walking. It checks the shift-rigidity criterion against every equal-size bipartition of [n] for even n≤14 and every positive t<n. All checks use explicit failures and survive optimized Python.

Reproduce:

    python3 check.py > certificate.json
    python3 -O check.py > certificate.optimized.json
    cmp certificate.json certificate.optimized.json

`twins.json` is an explicit input certificate, not trusted output: every edge, vertex, shift, and common root word is rechecked. `find_double.py` records its discovery search; `explore.py` records the initial rotation experiment.

Open constructive targets: classify which interval sizes admit double-placement square words; optimize connector height when additional connectors beyond the minimum count are allowed; replace repeated rotation search by a finite bank of root-word blocks with authenticated seams and disjoint occupation sets. These are sharply posed extensions. None follows solely from closing the residue8/9 clocks, and no consequence for arbitrary Collatz trajectories is asserted.

Verified finite totals: 9,331 root-word placement comparisons; 57,658 shift-partition comparisons; 19,691 rotation iterations in total, with at most271 for a single extension. Exact admissible starting vertices for the retained root words are40 at n=44,20 at n=72, and424 at n=527. The normal and optimized certificate files are byte-identical.
