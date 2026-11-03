# From the65 path chart to all three ternary outputs, with one shared port

**FINITE-EXACT** literal branch certificates and complete stated search;
**PROVED** the general native map to the inherited tribonacci family and
the subsequent zero-offset continuations. These do not establish a
repeatable three-branch return or universal odd-spectrum coverage.

## The projective decoder selects actual legal operations

Take E8 with its original labels0,...,7. Its rows are

    (126,252,248,240,224,192,128,1), H=65.

For each of its65 Hamilton paths, perform the native complement: preserve
the seven consecutive path arcs, reverse every other old arc. Then add
one vertex8, letting a mask specify which old vertices point into8.
This language preserves old labels, including physical port6, but changes
old edges. It has no helper deletion and uses exactly one birth.

An exhaustive search of65*256=16,640 labelled programmes gives37 strong
H195 targets,18 strong H197 targets, and2 strong H199 targets. These are
programme/labelled-target counts, not isomorphism-class counts. The
complete per-path counts are in `native_probe.json`.

The following choices have **the same retained output port6**, with J=0
on every branch. Field addresses refer to the fixed binary basis in P4;
the paths themselves authenticate legality without needing that basis.

| Target | Address | Legal path in the original E8 | H after complement | Birth mask | Unique in-neighbour of port6 |
|---:|---:|---|---:|---:|---:|
|195=3*65|60|(3,4,5,6,7,0,1,2)|27|229|5|
|197=3*65+2|41|(1,4,6,7,0,2,3,5)|115|254|4|
|199=3*65+4|32|(6,7,0,1,2,3,4,5)|17|140|8|

The final literal rows, in order0,...,8, are

    195: (258,4,257,23,39,335,415,319,26),
    197: (4,273,266,291,333,279,431,319,1),
    199: (2,4,265,275,39,15,191,319,115).

We retain the original vertex0 as well, with J respectively925,1287,1703.
It would be false to say that every old root response was preserved.
If the intended mark is6, it is the SAME physical vertex before and after
all three programmes; its unique predecessor changes and is recorded in
the table. The input's J at6 is484, so this is an actual marked improvement
under a changed old relation, not a phase-neutral source substitution.

For each target the zero-slack port and the inherited ear theorem give
strong induced continuations with counts

    195*3^j, 197*3^j, 199*3^j,  j>=0,

at order9+2j. At each ear retain its first new vertex as the next marked
zero-slack guard; the complete earlier vertices survive. The proof is
H'=3H when J=0, and the new guard again has J=0. These continuations are
inherited consequences of the guard theorem, not newly discovered scalar
values or further nonzero-offset branches.

## An unexpected exact meeting with the March tribonacci work

The199 intermediate is not an arbitrary17-path graph. Its first six
vertices have rows

    Q=(2,4,9,19,39,15),

which is literally the prime quotient of the incoming
[two-wing source family](../defect_offset_dynamics_20261103/WINGS.md).
The remaining old vertices form the source chain6→7→Q. Thus the finite-
field path address selects a native move into an actual earlier source
family with two explicit exterior vertices retained.

There is an all-size identity behind it. In E_n, n>=4, choose

    p=(n−2,n−1,0,1,...,n−3).

Its native complement is exactly the ordinal join

    (n−2) → (n−1) → T*_(n−2),                            (10)

where T*_m has forward consecutive arcs i→i+1 and every other arc
backward. Proof: the selected first edge n−2→n−1 stays; every other
incidence from these two vertices to the core is reversed to point into
the core, except n−1→0, which was already selected in that direction.
Inside the core precisely the consecutive increasing edges are protected.
This checks every pair and preserves the original vertex labels.

The H(T*_m) tribonacci recurrence was already proved in March, and is
recovered with its correction lineage in
[THM-4731, §6](../../../01-canon/theorems/THM-4731-tiling-line-metagraph-self-lines-by-twisted-burnside-sigma-twins-depth-and-the-tribonacci-partner.md).
It is inherited here, not re-proved as a new discovery:

    h_1=h_2=1, h_3=3, h_m=h_(m−1)+h_(m−2)+h_(m−3).

Equation(10) therefore sends the dyadic count2^(n−2)+1 to the actual
tribonacci count h_(n−2) without deleting a vertex. At n=8 this is65→17;
at n=9 it is129→31. The source-chain factor contributes1 to the raw
Hamilton count. The intermediate is not strong, which is allowed by the
declared native-plus-birth programme; the final three targets are strong.

This is a genuine explanation connecting a path-set field chart, a
metagraph operation, an old recurrence and the new wing interface. The
two recurrences remain different: the old conjecture identifying
tribonacci counts with2^(m−2)+1 already fails at m=7 (31 versus33).
The map changes H explicitly; it does not equate those sequences.

## Independent audit and stopping boundary

The main calculation uses subset Hamilton DP on each literal birth.
`native_audit.py` independently enumerates all old-vertex permutations,
cuts them into two valid directed segments, and assembles the complete
quadratic response to the birth mask. It checks every one of the16,640
counts and the strongness filter. It then independently counts rooted
ordered path pairs on the three displayed targets by literal9-vertex
permutations and cuts, verifying J(6)=0. No isomorphism quotient is used.

The search has a sound optional pruning rule H(intermediate)>199 by
induced-extension monotonicity; on this actual65-path bank no path is
pruned. It therefore checks all masks on all paths.

What is still missing: a family of sources closed under these three
programmes, or a justified return of the195/197/199 outputs to the
required input charts. Zero-slack continuation supplies repeated factor3,
not repeated offsets2 and4. The new target is to generalize the actual
operation in(10) and the authenticated birth interfaces, rather than
replacing a path-set projective action with an unproved native symmetry.

Reproduce from the repository root:

```text
python3 05-knowledge/results/six_roots_hidden_bits_20261103/native_probe.py
python3 05-knowledge/results/six_roots_hidden_bits_20261103/native_audit.py
```
