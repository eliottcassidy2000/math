# Recognizing skew conference tournaments from an unlabelled six-path histogram

2026-11-27. For order 24, the frequency of a single exact
path count distinguishes the nine classes up to switching,
permutation and converse. The corresponding order-28 test
begins with the same third shell. At order 24 the two codes
from the extremal cuts also give exact partner-orbit counts
for a fixed Golay code.

## 1. The third shell at every order at least 24

Use the [six-run moment identity and minimum edges](../pants_extremal_cuts_golay_20261127/README.md).
Write $P_6$ for the number of ordered directed paths on six
distinct vertices (no singleton ordering factor), and put
$$
 L=n(n-2),qquad S=\sum_i ((A\varepsilon)_i^2-1)^2,qquad
 P_6=\frac{(n)_6+n(3(n-1)^2-2(n-1))-S-n(2n-3)}{32}.
$$
Cuts are counted modulo antipodes. Let $s_0,s_1,s_2$
be the values of $S$ in the first three shells, and
$U=s_0-3s_1+3s_2$.

Fix one minimum edge, normalize its source endpoint to the
all-positive vector, and write $r$ for its direction. At a
third-shell cut let $B=\{a,b,c\}$ be the three flipped
coordinates off $r$. Put $Q=Q_{rabc}$ in this normalization
and $K=A_{ac}A_{bc}+A_{ab}A_{cb}+A_{ba}A_{ca}$. Then
$$
 S=U-192D,qquad D=n-1+K\pm Q,
$$
where $K=-3$ for a cyclic triple and $K=1$ otherwise.
The source endpoint takes the plus sign. Indeed, the discrete
third difference of the quartic is $-192\sum_{x\notin B}Q_{B\cup\{x\}}$.
In the sum over columns, column $r$ contributes $n-4$.
The remaining columns of the source-normalized matrix sum to
1, giving $Q$; removing the rows of $B$ contributes $3+K$
by column orthogonality. For example, the excluded $b$-row
sum over columns other than $r$ is
$\sum_{j\ne b,r} A_{aj}A_{cj}=-1-A_{ab}A_{cb}$.
Toggling coordinate $r$ in the base cut negates only the
$Q$ term.

**No other shell contributes to these bins.** Distinct minimum
edges are at distance $(n-2)/2$, so radius-three shells are
disjoint. For $d\ge4$ the largest absolute row sum $M$ is at
most $n-9$. When $M\le n/2$, $S\le(n^2/4-1)L$; otherwise
$S\le g(M)=(M^2-1)^2+((n-M)^2-1)(L-M^2+1)$ as in the
previous note. This cubic is convex on the relevant interval.
Writing $x=n-24\ge0$, the third-shell lower bound
$((n-7)^2-1)^2$ exceeds $g(n/2)$ by
$3x^4/4+89x^3/2+905x^2+6670x+7440$, and exceeds
$g(n-9)$ by $8x^3+384x^2+4864x+8448$. The earlier
second-shell separation covers cuts of distance at most two.

**The unlabelled histogram recovers a joint local census.**
Write $p=|\operatorname{Pf}(A[S_4])|$ for a four-set $S_4$.
For each of its four roots, the complementary triple is
cyclic exactly when $p=3$, and $|Q|$ is unchanged.
Thus the four-set contributes four times at each of
$D=n-4\pm|Q|$ for $p=3$, or at $D=n\pm|Q|$ for
$p=1$, double when the bins coincide. By the
[quartic residue and bound](../pants_skew_conference_quartic_detection_20261127/README.md),
all the indices $j=D/8$ are integral. The $p=3$ bins
are paired by $j\mapsto (n-4)/4-j$ on
$0,\ldots,(n-4)/4$. The $p=1$ bins are paired by
$j\mapsto n/4-j$ on $1,\ldots,(n-4)/4$. The two
reflections form a path in the order
$0,(n-4)/4,1,(n-4)/4-1,\ldots$, with one terminal loop.
After dividing frequencies by four, peel from $j=0$:
the remaining frequency at a bin is the contribution of its
yet-unknown edge, or twice the loop contribution. This
inverts the histogram to the exact joint distribution of
$(p,|Q|)$, without naming any cut.

For $8\mid n$, the first third-shell bin $D=0$ has
frequency four times the number of tetrads in the two
[minimum codes](../pants_paired_golay_leech_codes_20261127/README.md)
combined. Moreover the two cut parities split that count
exactly. The $p=1$ bound precludes zero; at $p=3$ the local
row products of $I+A$ on the four-set are all $+1$ in a
source normalization with root in the set. The four-set is
a word of the source code iff $Q=n-4$, by self-duality and
the constant product test from the previous note. In the sink
code the sign is reversed: its normalized generators
$h_i+e_i+e_r$ toggle the product on rows outside the set,
but not on any of the four rows of the set;
the resulting product is constantly positive iff $Q=-n+4$.
Thus a qualifying four-set is counted from the *opposite*
endpoint at each of its four roots. Flipping three coordinates
changes parity, so the cut itself has the parity of the
qualifying code's endpoints. With all-positive cut a source
minimum those are even for the source code and odd for the
sink code.

## 2. Order 24: one frequency and two partner orbits

[Spence's published representatives](https://www.rangevoting.org/SpenceMat.html)
exhaust ordinary Hadamard equivalence. The finer signed
congruence classification also has 16 classes
$[Araya et al., Table 1](https://arxiv.org/pdf/2402.11074)$,
so these representatives exhaust switching and permutation.
The replay normalizes every target root and tests adjacency
isomorphism to determine the converse pairs.

Here $U=94464$, and $D=8$ gives $S=92928$ and
$P_6=3026628$. The replay verifies the count also by
enumerating actual six-vertex paths for one cut. The nine
frequencies are:

| Spence blocks | Frequency at 3,026,628 |
|---|---:|
| 1, 6 | 800 |
| 2, 4 | 684 |
| 3, 5 | 864 |
| 7 | 544 |
| 8, 9 | 696 |
| 10, 11 | 528 |
| 12, 13 | 504 |
| 14, 15 | 1320 |
| 16 | 0 |

The code pairs from blocks 14, 15 and 16 are the only ones
having a Golay constituent. The first two are converses of
one another, and their non-Golay code has 66 tetrads.
The Paley block 16 has two Golay constituents. Consequently,
for a fixed Golay code on the 24 contacts there is one
compatible Golay-partner orbit and one other Type II partner
orbit. The compatibility test here includes intersection
equal to repetition and the weight-12 projection condition
of the previous note. Orbit-stabilizer in $M_{24}$ gives
$244823040/6072=40320$ Golay partners and
$244823040/1320=185472$ mixed partners, totaling 225792.

In the mixed pair the 66 tetrads are precisely the unions
of two out of twelve matched pairs. Their incidence recovers
the pairs. Enumerating normalized core isomorphisms counts
1320 permutations preserving the ordered codes; their
induced action on pairs is faithful and, in the field order
given by the replay, equals the fractional-linear action
of $\mathrm{PGL}_2(11)$ on the projective line. The
double-Golay ordered group is $\mathrm{PSL}_2(23)$ by the
[previous result](../pants_paired_golay_leech_codes_20261127/README.md).

The same integral-coordinate construction therefore gives a
Leech lattice and a Niemeier lattice on the same contacts.
For any length-24 Type II code this construction is even and
unimodular: write $x=m\mathbf1+2c+4k$ with $m=0,1$,
and the sum condition is $\sum k_i\equiv m\bmod2$. Thus
$\|x\|^2\equiv0\bmod16$; there are $2\cdot4096=8192$ cosets
of $\sqrt2 D_{24}$, whose covolume is 8192. The two lattices
have the same common sublattice of index 2048 as the
double-Golay pair. On the non-Golay side the roots are four
entries $\pm2$ on a tetrad, scaled by $1/\sqrt8$, with an
even number of minus signs. A lone $\pm4$ fails the sum
congruence; odd vectors cannot have norm two. Rotating
each matched coordinate pair to its normalized sum and
difference axes sends each root to a two-coordinate root
in one sector only: the two occupied pairs either both
have equal signs or both have opposite signs. Hence the
root system is $D_{12}\oplus D_{12}$, with 528 roots.

## 3. Order 28: one collision, resolved by the fourth shell

The same published page contains 65 order-28 blocks. They are
pairwise signed-inequivalent by the replay, so they exhaust the
65 signed classes in Araya et al. The replay also shows
that they give 41 classes up to converse.
The complete third-shell histogram distinguishes 40 of
these groups, with just one two-way collision, blocks
44 and 48 (one-based numbering among the order-28 blocks).
At order 28 the fourth shell is also disjoint in score
from all others. Its largest row magnitude is 19, so
the leading contribution is 129600. The other 27 values
$y_i=(A\varepsilon)_i^2-1$ add to 368. For the allowed
values $0,8,24,48,80,\ldots$ one has $y_i^2\ge32y_i-192$.
Therefore the fourth shell lies at $136192\le S\le159040$,
the upper bound being $g(19)$. The third shell starts
at $S\ge193600$. At distance at least five, odd $M\le13$
gives $S\le168\cdot728=122304$, and $M=15,17$
gives $S\le\max(g(15),g(17))=135744$.
Radius-four balls about different minimum edges are disjoint.

The exact fourth-shell score $S=139840$ gives
$P_6=8474100$; the two frequencies are 110760 and 111064.
Thus the unlabelled six-path histogram identifies the class
also at order 28, even though the entire third shell no
longer does so. In fact **two frequencies suffice**: the single
third-shell bin at 8,472,180 paths (D=32) already gives
40 different values, so its only class collision is the same
one. Use the 8,474,100-path frequency only to break that tie.

Replay with Python 3.10+:
```sh
curl -L https://www.rangevoting.org/SpenceMat.html -o SpenceMat.html
python replay.py SpenceMat.html
```
The replay uses the published matrix order and source-border
convention, and counts each cut once modulo antipodes.
