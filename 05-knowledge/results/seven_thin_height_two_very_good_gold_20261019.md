# Every non-total seven-thin height-two poset has a very good pair and weak Gold

Date: 2026-10-19. Status: all-size structural theorem, with an elementary set-family proof. The Gold consequence uses the explicitly attributed good-pair synthesis below. Ordinary balance for height-two posets is already known; the conclusion here is the stronger very-good-pair structure and its weak Gold consequence. No historical priority claim is made.

## 1. Statement, conventions and isolated vertices

Let P be a finite, non-totally-ordered poset with no three-element chain. If every vertex has at most seven incomparable neighbors, then P contains a very good pair. Consequently P admits weak Gold in the convention of the [good-pair Gold synthesis](good_pair_gold_synthesis_20261019.md): after two comparisons the initial extension count E is at least the sum of the two successive surviving counts on every realized path. Determined or repeated comparisons are allowed there.

A very good pair is an incomparable pair which, in P or its dual, has equal strict downsets and both directed upset differences are chains. The definition and balanced-pair consequence are attributed to Zaguia, Definitions 1,5 and Theorem 2, [author-deposited primary paper](https://arxiv.org/pdf/1610.00809v3). The linked synthesis proves the weak Gold consequence, combining the attributed decreasing insertion-bin and adjacent-threshold arguments. We use that result, rather than claiming that Zaguia's theorem itself states Gold.

Partition P into A, all minimal vertices including every isolated vertex, and B, the remaining vertices. Every vertex in B is maximal, and every relation goes from A to B. Both levels are antichains. More generally isolated vertices may be assigned to either level; the argument only needs such a two-antichain partition. By passing to the dual and exchanging levels if necessary, assume a=|A|≥b=|B|. This exchange remains legitimate with isolated vertices: the proof does not require every B-vertex to have a predecessor.

Empty and singleton posets, and the two-element chain, are totally ordered and excluded from the theorem. These are precisely the total orders of height at most two; no two-comparison Gold assertion is made for already sorted total orders. If b=0, P is an antichain, and because it is non-total it has at least two vertices, which are exact twins. If a=b=1, non-totality says the two vertices are incomparable, hence again exact twins. These cases include the isolated-vertex edge cases. Henceforth a≥2 and b≥1.

## 2. Missing-incidence families

For each v in A let Z_v={w in B : v is incomparable with w}, and define missing columns dually. The exact incomparability degrees are

    d_inc(v)=a−1+|Z_v|,
    d_inc(w)=b−1+|Z_w|.

Thus missing rows have size at most 8−a and missing columns at most 8−b.

Two rows R,S are called near when |R\S|≤1 and |S\R|≤1. Their vertices have equal empty downsets, and the two upset differences have exactly these cardinalities. Since B is an antichain, near rows give a very good pair. Near columns give a very good pair by duality. Repeated rows or columns are near.

Suppose, for contradiction, that neither level has such a pair. All rows and columns are distinct, and every two rows (or columns) have at least one directed difference of size at least two. Call such a family separated.

## 3. Reduction to at most five rows, with no circular bound

We have chosen a≥b. If a≥7, the row-size bound gives |Z_v|≤1 (and is negative for a>8, already impossible). Any two sets of size at most one are near, contradicting a≥2. Hence a≤6 and therefore b≤6.

If a=6, every row has size at most two. Separated two-element sets are pairwise disjoint, so there are at most floor(b/2)≤3 of them. Among empty sets and singletons there can be at most one row, since any two are near. Thus there are at most four rows, contrary to a=6. Consequently a≤5. This uses no prior bound on b beyond b≤a.

## 4. Complete small-family argument

### At most three rows

On a universe of at most one element every two subsets are near, so b=1 is impossible for a≥2. On at most two elements, the only separated pair is the empty set and the full two-element set. Thus a=2 forces b=2 and these two rows; their columns are equal, a contradiction.

On a universe of at most three elements, all subsets of size at most one are mutually near, and all subsets of size at least two are mutually near. A separated family therefore has at most two members. This excludes a=3 (and also larger a with b≤3).

### Four rows

Pad the column universe to four elements if necessary. Among rows of size at most one, at most one is possible. Separated two-element rows are disjoint, so there are at most two. Among rows of size at least three, at most one is possible. Four separated rows must therefore consist of one small row, two disjoint doubletons, and one large row.

The doubletons partition the four-element universe. Every triple contains one doubleton and is near to it, so the large row must be the full four-element set. Any singleton belongs to one doubleton and is near to it, so the small row must be empty. Up to column permutation the unique family is

    ∅, {1,2}, {3,4}, {1,2,3,4}.

Columns 1 and 2 are equal, contradicting the absence of a near column pair. If padding was used, the forced full row already shows that no padded point can actually be absent, so this argument also covers b<4.

### Five rows

Now a=5, so every row has size at most three. Pad to a five-element column universe if necessary. There is at most one empty-or-singleton row. The doubletons are disjoint, so there are at most two. Two triples are separated exactly when their intersection has size one, equivalently their two-element complements are disjoint. There are therefore at most two triples.

Five separated rows must have sizes

    (0 or 1), 2, 2, 3, 3.

Write the doubletons as {1,2} and {3,4}, with remaining point 5. A triple cannot contain a doubleton, since then they would be near. Therefore each triple contains point 5 and exactly one point from each doubleton. The two triples must intersect in exactly one point, so their choices from both doubletons are complementary. Relabeling within doubletons makes them {1,3,5} and {2,4,5}.

If the small row is a singleton, it cannot belong to either doubleton, again because it would be near to that doubleton. Thus it must be {5}. We obtain precisely two families up to column permutation:

    ∅,   {1,2}, {3,4}, {1,3,5}, {2,4,5};
    {5}, {1,2}, {3,4}, {1,3,5}, {2,4,5}.

In either family columns 1 and 2 agree in the doubleton row and differ only by exchanging the two triple rows. Each directed column difference has size one. They are near, the final contradiction. Every point occurs in one of the forced doubletons or triples, so padding cannot introduce an unused point; b<5 is included.

This exhausts all possibilities. Some same-level very good pair exists unless the already handled two-point antichain supplies a cross-level twin pair. The theorem follows, and the cited good-pair synthesis gives weak Gold.

## 5. Exact finite controls and scope

The accompanying standard-library program `../experiments/verify_height_two_seven_thin_gold_20261019.py` independently enumerates separated row families for 2≤a≤6 and 1≤b≤a, subject to the seven-thin row-size budget. Only (a,b)=(2,2),(4,4),(5,5) survive, with respectively 1,3,60 labeled families. It checks a near column pair in every survivor and independently verifies that the four-row families form one orbit and the five-row families form precisely the two displayed orbits, each of size 30. Normal and optimized executions are recorded separately; all guards remain active under `-O`.

These small checks are regression controls, not a replacement for the all-size elementary proof. The column degree bound is not even needed after a≥b is chosen: row thinness and the symmetric near-pair condition suffice. The statement allows arbitrary isolated vertices and imposes no c=0 critical-pair assumption, strengthening the scope of the [earlier c=0 height-two result](seven_thin_height_two_c0_very_good_pair_20261019.md). It makes no assertion for higher-height posets or thinness eight or nine, and does not settle the general seven-thin conjecture. No novelty claim is made for the structural argument or the inherited Gold mechanism.
