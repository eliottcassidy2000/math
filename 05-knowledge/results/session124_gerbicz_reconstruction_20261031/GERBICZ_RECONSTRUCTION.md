# Gerbicz's 2018 square-sum theorem: independently replayed proof data

**CITED classical theorem; independent certificate reconstruction, not a new
existence theorem.** The ordinary square-sum conjecture in the user's request
has already been solved. Here we recover and verify its actual construction,
not merely a bibliographic claim.

For the simple undirected graph G_n on {1,...,n}, with distinct a,b adjacent
iff a+b is an integer square, a Hamiltonian cycle exists for every n>=32.
Consequently a Hamiltonian path exists for every n>=32; supplied literal
paths also cover15,16,17,23 and25..31. The singleton n=1 is trivial if included
in the convention. This reconstruction does not itself prove nonexistence
at the remaining small orders; the existing repository small-order census
and endpoint-aware obstructions supply that separate classification.

## Source acquisition, attribution, and boundaries

Robert Gerbicz's author comment dated2018-01-21 at
https://oeis.org/A090461 states the path/cycle thresholds. Berlekamp–Guy,
*Fibonacci Plays Billiards*, arXiv:2002.03705, appendix and reference[2],
links Gerbicz's original code archive:
https://drive.google.com/file/d/1XpYRdJR00CldIfssjDW1clkh-HdPA708/view .
The source forum https://mersenneforum.org/showthread.php?p=477787 returned403
on2026-10-31; the Drive preview was inaccessible to the web tool. Ordinary
public download at
https://drive.google.com/uc?export=download&id=1XpYRdJR00CldIfssjDW1clkh-HdPA708
succeeded without escalation. The ZIP contains only squares.c.

ZIP SHA256:05ea9f972a8b470fee84efab821dd1fe93e9001e5d310dd36c3b1240be5bdede.
C SHA256:602c1919ee51101519a2f19f2c1817407006824f6f09c85a3b772a818a0b736e.
C size:17321973bytes. Its header attributes Robert Gerbicz,2018.
No explicit redistribution license was found. Keep source and archive in
research scratch, not a published packet. The independent verifier takes
this original file as a local input; it compiles/executes none of its code.
It extracts only three integer initializer arrays under a strict whitelist.
The review/publication package is NOT offline self-contained: reproduction
requires obtaining the external archive separately from the pinned public
source. The external-data dependency remains explicit even though the
arithmetic proof and certificate checker are independently reconstructed.

## The invariant

A nice pair at n comprises paths P_0 on[1,n] and P_1 on[1,n+1]. Both start1.
Their last vertex is3 for odd length and8 for even length. Every common
vertex occupies positions of the same parity. Positions begin at1.

Both paths are cycles after adding the closing edge, because1+3=4 and1+8=9.
This differs from the repository's later protected-seed B family (odd end10),
whose extra prefix/protection needs more elaborate final states.

## Exact lift

For each offset c in[-24,24], transform a chosen input path by

    x_j -> 49 x_j + (-1)^(j+1) c.

For each old vertex x<=n, the49 choices of c produce exactly
[49x-24,49x+24], independent of which input member is selected because the
pair shares positional parity. Internal square edge sums are multiplied
by49 and remain squares. The extra vertex n+1 has position parity n+1, not necessarily the final
position. Counting odd/even positions and subtracting the shared parity
counts for[1,n] forces this. Its lifted value is
49(n+1)+(-1)^n c.

For target length M=49n+d+e,24<=d<=72,e in{0,1}, select input member1
at c exactly when49+(-1)^n c<=d+e. This includes precisely the required
consecutive extra labels. Add singleton vertices1..24. All atoms partition
[1,M], with no repetitions. Each atom may be reversed.

## Finite universal transition certificates

The array `data` has196 words, one for every(d,p,e), where p=n mod2.
Each has49 complete blocks and24 singleton atoms. Metadata identifies its
state. The verifier checks all offsets and singletons occur exactly once,
the forced member choice, every connector square, and output endpoints1
and3/8 according to output parity. Endpoint formulas depend only on p,
c and selected member—not on the internal ordering or magnitude of n.

For each paired pair of words, all singleton position parities agree.
For a block starting in position parity s and having length parity L,
the output-position shift relative to its original input position is
s-1 for a forward block, or s+L for a reversed block (mod2).
These shifts agree for each c across the two output words. Thus every
common lifted vertex has matching output parity, regardless of where it
occurs in an arbitrary admissible input pair. This proves nice-pair
preservation for every n with the specified parity, not just tested n.
The finite check therefore certifies98 universal paired transitions.

## Bases and complete induction

Array `S` contains1992 pairs, at every n in[41,2032]. The checker validates
exact permutation coverage, every square adjacency and closure, endpoint
normalization, and common positional parity. Array `smallseq` supplies
cycles32..40 and the listed smaller paths.

For any N>=2033, put m=floor((N-24)/49),d=N-49m. Then
24<=d<=72 and41<=m<N. Strong induction applies the certified transition at
input m to construct the nice pair at N. Since2033=49*41+24, there is no
uncovered interval. Closing either member proves cycles at every N>=41;
the literal small cycles cover32..40. This closes the full original
threshold without a separate enormous finite bridge to128453.

## Replay and relation to repository proof

Run `python verify_gerbicz_data.py /path/to/squares.c` and repeat with `-O`.
The checker uses explicit exceptions, not removable assertions. It performs
4205601 checks. Normal and optimized receipts agree. The mathematical
universal conclusion uses the induction and parity/coverage proof above,
not extrapolation of finite constructions.

This packet is a primary-source certificate replication, distinct from the
repository's separately verified ordinary/seeded construction chains. It
does not depend on copying those historical bundles.

## Independent audit

The separately authored `check_independent.py` and `AUDIT.md` validate the
base and universal transition obligations,206 expanded paths and five
negative controls. Both producer and auditor reject a source file whose
SHA256 differs from the pinned original before parsing it. Current exact
counts appear in the accompanying receipts.
