# Frieze coordinates, legal Collatz words, and the first-reset-2 boundary

**Status.** PROVED: the elementary marked-coordinate embedding, single-ear
obstruction, guarded positive-quiddity descent theorem, excluded dyadic cell,
and explicit common-future families below. CITED: the classical polygon/frieze
and Grassmannian cluster interfaces. FINITE-EXACT: the declared experiments.
No new universal coverage, completed basin, or solution of the first-reset-2
boundary is claimed.

## 1. Inheritance and the connection being tested

The closest arithmetic mechanisms are ordered carriers and native source
guards in [Collatz carry interfaces](collatz_carry_interfaces_20261004.md),
the simultaneous endpoint fiber in
[the 223/233 prime-partition note, section 4](collatz_prime_partition_223_20261007.md),
and the variable initial-one-run interface in
[reset-two debt families](reset_two_debt_families_20261004.md). The canonical
hostile is source 7: its actual initial valuations are (1,1,2), and its
first-reset-2 boundary does not inherit the simpler reset-at-least-3 rule.
The corrected near miss is replacing a guarded arithmetic map by a matrix
quotient. The least-used useful sidecar here is the chronological list of
prefix carries, together with its terminal halving cost.

Historical `tropical_cluster_H.out` and `tropical_cluster_dynamics.out` propose
cluster language for tournament polynomials; they are not dependencies of
the results here. A separate current positivity warning is
[THM-2873 — two-ray factorial-response TP3 curvature](../../01-canon/theorems/THM-2873-two-ray-factorial-response-tp3-curvature.md):
its positivity is created by a specified weighting, not inherited from an
arbitrary positive output measure. That theorem concerns a different object.

Our concept board is **marked prefix carries / positive minors / frieze ears /
ternary endpoint labels / dyadic source guards / source height**. The anchor
is actual reset-2 receipts; the niche is a lossless coordinate change; the
wildcard is whether a positive frieze mutation can be a legal word rewrite.

The primary interfaces are explicit. Fomin and Zelevinsky's
[Cluster algebras I, section 1](https://www.esi.ac.at/preprints/esi1023.pdf)
describes the Gr(2,n) relation
`[ik][jl]=[ij][kl]+[il][jk]` and triangulation coordinates. Morier-Genoud,
Ovsienko and Tabachnikov's
[SL2(Z)-tilings of the torus, Coxeter–Conway friezes and Farey triangulations](https://sophie-moriergenoud.perso.math.cnrs.fr/Publi/TorusTiling5.pdf),
introduction and section 5, relates positive integral polygon friezes,
triangle-incidence quiddities and the recurrence `V_(i+1)=a_i V_i−V_(i−1)`.
The calculations below are proved directly. We do not import a general
cluster-positivity claim as a Collatz guard or re-audit those papers.

| Source | Target and exact map | Preserved predicate | Lost data / required sidecar |
|---|---|---|---|
| Actual valuation word | Marked columns `(1,B_i/3^i)` | Positive ordered minors and exact Plücker relations | Terminal exponent requires `Q_r`; source requires its native cell |
| A fixed marked configuration | A different triangulation coordinate set | Every marked minor and hence the same arithmetic word | Do not discard vertex order or frozen boundary gaps |
| Polygon ear insertion | Same frieze transfer matrix | Integral recurrence transfer | It changes length, cost, carry and legal endpoint class |
| Positive quiddity, length at least 5 | A dyadic source progression above an explicit cut | Actual finite descent | Descending endpoint still needs a ROOT proof for standalone completion |
| Two endpoint-compatible words | A common-endpoint progression | Actual common future and smaller child | The ternary label and source/height guards remain explicit |

## 2. A lossless positive-coordinate realization

For `U(n)=(3n+1)/2^v2(3n+1)` and a word `w=(a_1,...,a_r)` of positive
integers, use chronological prefix carriers

    F_i(n)=(P_i n+B_i)/Q_i,
    P_i=3^i, Q_i=2^(a_1+...+a_i),
    (P_0,Q_0,B_0)=(1,1,0),
    B_(i+1)=3B_i+Q_i.                              (1)

Put `z_i=B_i/P_i`. Then

    z_0=0,
    Delta_i=z_(i+1)−z_i=Q_i/3^(i+1)>0.             (2)

For a nonempty word the columns `v_i=(1,z_i)` form a rank-two real matrix
whose every ordered 2-by-2 minor is positive:

    d_ij=det(v_i,v_j)=z_j−z_i>0,  i<j.             (3)

Direct expansion proves, for `i<j<k<l`,

    d_ik d_jl=d_ij d_kl+d_il d_jk.                 (4)

Thus a positive Plücker flip is an exact coordinate change on this same
configuration. This is a special marked locus in the positive Grassmannian,
not an assertion that all positive Grassmannian data represent Collatz words.
For example, gaps `(1/3,1/3)` are positive but would require `2^a=3`.

There is an all-length inverse, with an essential terminal coordinate:

    Delta_0=1/3,
    3 Delta_i/Delta_(i−1)=2^(a_i),  1<=i<r,
    a_r=log2(Q_r)−sum_(i<r) a_i.                   (5)

The right sides must be exact powers of two with positive exponents.
Equations (2) and (5) prove that the marked `(z_0,...,z_r,Q_r)` is lossless.
The script's `decode_chart` verifies all these conditions and re-encodes the
full chart. Empty words are a separate identity case. The rational data use
exact integers and fractions; booleans and floating-point aliases are rejected.

The terminal coordinate cannot be dropped: every one-letter word has chart
`(0,1/3)`. Word (1) is actual at source 3 and increases it to 5; word (2)
is actual at source 9 and decreases it to 7. More generally `B_r` omits the
last valuation. A frieze or positive character forgetting that coordinate
cannot decide even this one-step drift.

The frozen boundary coordinates of a polygon contain every adjacent gap.
Consequently a triangulation flip preserves all of (5), and hence `P,Q,B`,
when the terminal `Q` and marked order are retained. The
five-flip pentagon relation is checked exactly in the script. This is a
lawful way to re-express a word; it is not a compression of the frozen gaps
or an added source-coverage theorem.

To restore actual arithmetic, retain the original positive odd source and

    P_r n+B_r = Q_r mod 2Q_r.                      (6)

This is necessary and sufficient for the full formal word to be actual:
working backward from its odd final endpoint, every earlier inverse numerator
has its required factor of 3, because the cleared final identity reduces
successively modulo 3. Every intermediate integer is odd; positivity also
follows directly from the positive forward carriers. Equivalently, forward
induction peels the exact factors of two, since each remaining carrier has
odd multiplier and odd carry. For first-hit receipts, impose that all proper
prefix states exceed 1. A flip preserves these guards because it preserves
the word; positivity of the minors alone does not imply them.

## 3. Why a frieze ear is not a Collatz rewrite

Define the recurrence matrix

    M(a) = [[a,−1],[1,0]].

Matrix multiplication gives the elementary ear identity

    M(a) M(b)=M(a+1) M(1) M(b+1).                  (7)

In a triangulated polygon this inserts a new ear and increases the two
neighboring incidence numbers. Interpreting those entries as valuations
instead gives words `(a,b)` and `(a+1,1,b+1)`: length increases by one and
halving cost by three, so the Collatz slope is multiplied by `3/8`.

There is a stronger integer obstruction, valid for arbitrary positive `a,b`.
Every actual word whose last valuation is `b` has an odd endpoint `y` with

    2^b y=1 mod3,   hence y=(−1)^b mod3.           (8)

The two blocks in (7) therefore have disjoint legal endpoint classes modulo
3. **They cannot reach the same integer endpoint even from different integer
sources.** Common prefixes do not change this block statement. Appending the
same formal suffix cannot repair it: that suffix is injective as an affine
map, so equal final endpoints would require equal block endpoints.

The smallest example is

    F_11(n)=(9n+5)/4,
    F_212(m)=(27m+29)/32.

Equality would require `27m=72n+11`, impossible modulo 3. A positive
subtraction-free identity in recurrence coordinates is not a legal arithmetic
substitution. Multiple ears can return to compatible arithmetic fibers;
section 5 supplies a fully guarded example.

## 4. A positive all-height subclass and an excluded whole cell

A positive polygon quiddity records the number of triangles incident to each
of the `N>=3` marked vertices. There are `N−2` triangles, so its valuation
cost, when interpreted as a word, is

    A=3N−6,
    P/Q=3^N/2^(3N−6)=64(3/8)^N.                  (9)

This slope is strictly below one for every `N>=5` and above one for `N=3,4`.
Given such a contracting word with carrier `(P,Q,B)`, let

    C=max(1, floor(B/(Q−P)),
          floor((Q_i−B_i)/P_i) for 0<=i<N),
    r=(Q−B)P^(−1) mod2Q,
    n_0=the least integer n>C with n=r mod2Q.      (10)

Then for every integer `t>=0`,

    n=n_0+2Q t,
    U_w(n)=(Pn_0+B)/Q+2P t < n,                  (11)

the word is actual, and no proper prefix reaches ROOT. Formula (6) gives
the exact word, the prefix terms of (10) rule out early ROOT, and the term
`B/(Q−P)` gives strict descent. Thus `descent_row` is a total constructive
compiler for every supplied positive quiddity of length at least five.
It excludes only a finite head of that word's native cylinder. This proves
a descent rule, not ROOT for the union of its cylinders: induction or an
authenticated endpoint proof is still required.

This positive subclass has an exact, unavoidable coverage boundary. For a
polygon of at least four vertices, two adjacent vertices cannot both be
ears. Otherwise their ear diagonals cross. Thus no positive quiddity of
length `N>=4` contains adjacent entries 1, including across its cyclic seam.
Every positive odd `n=7 mod8` has its first two actual valuations equal to 1.
Consequently:

> No source `n=7 mod8` has an actual contracting positive-quiddity prefix.
> Extending the actual word to greater quiddity depth does not remove this
> obstruction.

For source 7, even the only length-three quiddity `(1,1,1)` fails, since its
third valuation is 2. Larger initial one-runs can realize that triangle word,
but it increases the source and cannot supply (11). In particular the direct
positive-quiddity scheme cannot pay the first-reset-2 problem at
`2^1459−1`, whose initial one-run has length 1458. This does not obstruct
other guarded Collatz rules or the marked coordinate representation of
section 2, which accepts every word.

**Monodromy hostile.** Requiring only `M(a_1)...M(a_N)=−I` is too weak:
`(1,1,1,1,1,1,1,1,1)` has that product, but is not the quiddity of a polygon.
The script uses constructive ear deletion to a triangle, not this insufficient
matrix test. This is the concrete distinction between recurrence transfer and
the positive polygon/non-oscillation condition in the cited frieze interface.

## 5. Guarded multi-ear common futures

For a word w write

    c(w)=B_w Q_w^(−1) mod 3^len(w).                (12)

Its integer odd endpoints occupy precisely that ternary class, with positivity
enforced separately. Two words admit arbitrarily large simultaneous positive
odd endpoints iff their labels agree modulo the smaller power of 3. This is
the nested-modulus form of the inherited simultaneous-endpoint theorem.
If `R=max(len(w),len(v))`, a compatible base pair extends by

    H(t)=H_0+2·3^R t,
    n(t)=n_0+2Q_w 3^(R−len(w))t,
    h(t)=h_0+2Q_v 3^(R−len(v))t.                  (13)

`join_row` also imposes a strict finite height cut making both sources and
proper-prefix states positive and avoiding early ROOT, and orients the row
so `h<n`. It constructs the forward ray `t>=0`; it does not claim this is
the complete set of allowed negative parameters.

An ear-generated block gives

    w=(2,2,2,3,1,2,5), v=(1,1),
    n=66817+262144t,
    h=495+1944t,
    U_w(n)=U_v(h)=1115+4374t.                     (14)

The marked first-reset-2 example is especially explicit:

    w=(1,2,2,2,2,2,2,1,7), v=(1,1,1),
    n=2583211+4194304t,
    h=7183+11664t,
    U_w(n)=U_v(h)=24245+39366t,  t>=0.            (15)

Both words in (15) are genuine positive polygon quiddities with matrix `−I`;
the nine-letter word is obtained from the triangle by six ears. The
endpoint congruence and actual source cell, rather than matrix equality,
authorize the join. A supplied ROOT certificate for the smaller h transfers
to n via the common future; the row does not supply that certificate itself.

This is deliberately not counted as added stopping coverage. The source
already has

    U_(1,2,2)(n)=(27n+23)/32<n                   (16)

throughout (15). Likewise (14) starts with valuation 2 and already descends.
These examples demonstrate that the arithmetic-compatible part of a frieze
fiber is nonempty, while a single elementary ear cannot traverse it.

## 6. Finite experiment and remaining obligation

The companion script uses exact rational arithmetic and no external package:

    python -B -X utf8 04-computation/experiments/frieze_collatz_reset2_20261007.py
    python -B -O -X utf8 04-computation/experiments/frieze_collatz_reset2_20261007.py

Its declared universe is: all 1,365 words of lengths 0..5 over valuations
1..4 for marked round trips; 1,280 pentagon flips; all 256 single-ear pairs
`1<=a,b<=16`; and every descendant of block (1,1), by adjacent interior ear
insertions, through block length 9. The last universe has 626 blocks, with
counts `1,1,2,5,14,42,132,429`. Prepending a fixed 1 makes a marked quiddity
in every case; these are not all marked quiddities of their sizes. Contracting
rows are independently replayed at parameters 0,1,7. The first compatible
block and first compatible prefixed reset-2 word within this declared universe
are (14) and (15), respectively. The script additionally replays 160 join
instances, 512 members of the excluded cell, and 32 variable one-runs.
Malformed source/parameter/record and early-ROOT controls remain active under
optimization. Normal and optimized output must match the saved `.out`.

The useful survivor is a lossless positive coordinate calculus with explicit
arithmetic labels, plus a genuine guarded descent subclass. The missing
coverage coordinate is not positivity of cluster variables: it is the actual
word's source guard and a paid legal successor outside the no-adjacent-ears
subclass. A future rewrite that claims residual coverage must preserve those
coordinates or explicitly compile a new common-future incidence relation.

**Draft correction.** The initial one-letter hostile used source 5 while
calling its valuation 2; its actual valuation is 4. The corrected pair is
`3 --(1)--> 5` versus `9 --(2)--> 7`, and the checker now replays the stated
words, rather than merely comparing unlabelled one-step drift signs.
