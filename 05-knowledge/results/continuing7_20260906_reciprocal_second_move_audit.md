# Independent audit of the second reciprocal transposition

**Verdict: PASS.** The universal N_p>=2 obstruction, the retained-fixed-point
obstruction for N_p>=1, and the all-prime p=37 mod360 failure family are
accepted as analytic statements in their declared scopes. The exact p37,
p43, p113 and p197 controls pass by a separate geometric implementation.
No mathematical repair was requested. The producer's final direction
clarification in the last family witness was read and accepted.

Audited producer:
[proof](continuing7_20260906_reciprocal_second_move.md),
[source](../../04-computation/continuing7_20260906_reciprocal_second_move.py),
[certificate](continuing7_20260906_reciprocal_second_move_certificate.json).
These concern the p-point integer permutation board after the prescribed
first output exchange0↔1, followed by one ordinary output transposition.
No intervening affine move or different first swap is included. This is
not a2p-point no-three-in-line theorem.

## 1. Independent proof reconstruction

For the involution g with g(0)=1, g(1)=0 and g(x) the native reciprocal for
2<=x<=p−1, every old triangle has exactly one special anchor P=(0,1) or
Q=(1,0). Three remaining points cannot be collinear, since an affine line
over F_p meets xy=1 in at most two points. Both anchors lie on x+y=1,
which contains no remaining point. Integer collinearity implies modular
collinearity here, while the converse is never used.

The P secants produce disjoint endpoint pairs X_i on the row shore and
Y_i=g(X_i) on the output shore; their transpose supplies the Q triangles.
The output supports are {1}∪Y_i and {0}∪X_i. The sets X_i are disjoint
for distinct i, and so are the Y_i. The fixed point F=(p−1,p−1) cannot
lie in one of these triangles: the integer line through P,F has equation
(p−1)(y−1)=(p−2)x, forcing x=p−1 for any remaining point in the native
square. The other modular fixed point(1,1) is absent.

If X_i meets Y_i, the corresponding reciprocal endpoint would either
be fixed or the two endpoints would be transposes. The latter makes its
line equal its transpose and hence puts P,Q and a positive reciprocal
point on x+y=1, impossible. Similarly, for i!=j, X_i∩Y_j has at most
one element: its points lie on two distinct integer lines, the P secant
i and the transposed secant j. Inversion sends X_i∩Y_j bijectively to
Y_i∩X_j. These are the exact geometric ingredients needed by the proof;
arbitrary triangle supports would not suffice.

A successful swap must hit every old triangle's output support. For N>=2,
a swap involving exactly one special label cannot hit the disjoint
opposite-anchor endpoint pairs. Swapping both special labels restores
the original diagonal triple(0,0),(1,1),F. With neither special label
moved, N>=3 is impossible by disjointness. For N=2, a two-label cover
must choose one element from each cross intersection. The intersections
have size at most one and inversion exchanges them, so the two labels
must be v,g(v). That move creates(v,v),(g(v),g(v)),F, three distinct
diagonal points. Thus all candidates fail.

This proves failure of every second transposition, whether or not it is
required to be monotone. It does not constrain longer sequences or
completed affine moves. When N>=1, a swap involving p−1 misses one of
the two disjoint supports of any packet, giving the separate fixed-point
obstruction. In particular N=1 is necessary, but not sufficient, for a
successful second move of an initially failed board.

The independent finite incidence referee exercises a subtle branch not
present in the named prime examples:144 admissible abstract two-packet
patterns contain96 non-anchor inverse-pair covers. Every such cover
creates the diagonal trap. These are abstract lemma controls, not claimed
reciprocal boards or a replacement for the universal proof.

## 2. Uniform family and normalization

For p=37+360r, r a nonnegative integer, the displayed coordinates

```text
X1=24+216r, Y1=17+160r,
X2=30+270r, Y2=21+200r
```

give the named transposed pair of old triangles at prime p. Their strict
ordering supplies nine and only nine possible covers of those two old
triangles. This argument allows additional packets throughout the
progression; it never asserts N_p=1 for the family.

The inverse-pair and undo moves account for three covers. The independent
symbolic calculation verifies all three displayed remaining witnesses,
their conjugates, and the entire nine-cover union. It checks exact
determinant zero, integer polynomial quotients(xy−1)/p, positive native
coordinate bounds, distinct row coordinates for every integer r>=0,
and that the two claimed unchanged points have neither moved output.
This is a polynomial identity and inequality proof, not finite parameter
sampling. In the third representative, the moved point is the first
unchanged point plus(1,1), as clarified in the final producer report.

The conjugation is typed correctly: the inverse/transpose of tau g is
g tau=(g tau g)g. Thus the transposed move exchanges g(u),g(v), not
necessarily the same two output labels. The direct geometry referee
also verifies this identity on all declared finite covers.

Infinitude is separately **CITED** Dirichlet, with the producer's primary
[1837 proof in English translation](https://arxiv.org/pdf/0808.1408v2).
The arithmetic gate gcd(37,360)=1 is checked. The coordinate obstruction
itself is self-contained and holds at every prime in the progression.
The inherited negative character(5/p) is consistent but is not used as
a substitute for the second-move proof.

## 3. Independent exact implementation

The referee imports no producer implementation. It reconstructs the native
reciprocal permutation by searching integer modular inverses. Its full
triangle engine groups all later-row vertices by reduced integer direction
from each least-row vertex. Consequently each actual triangle is listed
once. This is distinct both from the producer's literal three-subset loop
and from its modular-conic/off-anchor decoder.

It independently reconstructs1571 moved boards: all666 transpositions at
p37, all903 at p43, and the sole old-support cover at each of p113 and
p197. For every other transposition at the two larger primes, it checks
an explicit unchanged old triple directly on the moved board. Every
declared count, complete triangle set at a cover, safe permutation, and
packet is compared with the frozen certificate.

Results:

| Prime | Native packets | Covers of every old triangle | Successful second moves |
|---:|---:|---:|---|
|37|1|9|none|
|43|1|9|(0,13),(1,10)|
|113|2|1|none|
|197|3|1|none|

The positive p43 pair is retained and checked as complete safe boards.
It prevents overstatement of the obstruction as applying to every
one-packet case. The p113 and p197 controls verify separate branches of
the N>=2 theorem, rather than expanding a prime census.

For the all-r family the independent implementation uses SymPy determinant
and polynomial division identities, rather than the producer's own
Fraction coefficient operations. It also separately enumerates the
abstract two-packet incidence controls described above.

The resulting verifier passes **29,274 always-active exact gates** in
normal and optimized modes, with byte-identical raw LF output. It needs
the frozen producer source, output and JSON only as pinned input data.
Its default path works from the outside moments folder or when filed
beside the producer under04-computation; `--producer-dir` explicitly
selects another input directory.

## 4. Reproduction and pins

- [Independent referee source](../../04-computation/continuing7_20260906_reciprocal_second_move_audit.py)
- [Normal output](continuing7_20260906_reciprocal_second_move_audit.out)

```text
python continuing7_20260906_reciprocal_second_move_audit.py
python -O continuing7_20260906_reciprocal_second_move_audit.py
```

SHA256:

```text
producer source 1e1571440d0ea569b7e3f65ffdc1ebe9ae61b2ae3d1f11cc39a7a6667d4b708c
producer output 060493fb23ff746151420891fd9cce0df5eb3e0a7ebd7b4e0a56d3a84f2d328c
producer JSON   f5117590d335c3df7579b8bb483d2a0a9d0906ab67c06bbb1cfd01ff7a579c3f
producer report 36aec99ba703ccea5892262a0f7aa8a060921eeabeee5e336df8eb61c08863b3
audit source    1629c1121e90e3674b9574071ef493e522c67cdc81dae997f1750dc4417aab3f
audit output    66fea8b69a1ff6f99a5f5a0d24d7d3a1886a9fd477ccff976fc0ae5064c7fd8c
```

The producer files were not edited. Parent owns repository filing and
status promotion. No general no-three-in-line conjecture is closed.
