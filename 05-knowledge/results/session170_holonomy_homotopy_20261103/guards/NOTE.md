# Two reset words, a Tribonacci zero, and recoverable observer memory

**Status.** An exact rank classification for the inherited C/D transfer semigroup, a second universal reset, and a marked-observer obstruction. The classification uses a classical theorem on the integral zeros of the Tribonacci sequence, explicitly credited below. The full pointed profile does not lose the same information: its transfer is injective, and the two reset gaps have a first-derivative recovery rule. These results concern this particular physical attachment language, not arbitrary tournament extensions or the odd Hamilton-path spectrum.

## 1. Inherited physical interface

[S164's guard algebra](https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session164_emergent_20261103/profile_algebras/GUARD_ALGEBRA.md) supplies the following construction. Let q have unique in-neighbor rho in a tournament X. Write

    s = H_start(X,q),
    u = H(X) − s,
    d = H_start(X−q,rho).

Add a new vertex a with q→a and a→X−q. Either mark a as the next guard (C) or keep q (D). Both operations preserve strongness when X is strong. Their endpoint transfers, on column vectors, are

    C = [[0,1,0],[1,1,1],[1,0,0]],
    D = [[1,0,1],[0,2,0],[0,0,0]].                         (1)

Thus C(s,u,d)=(u,s+u+d,s) and D(s,u,d)=(s+d,2u,0). Every letter adds one vertex. Words are written in temporal order: DCD means first D, then C, then D, and its matrix is D C D. For a word w, M_w denotes its temporally ordered transfer.

For an old observer p≠q, let X_0=X, X_1=X^a_p and X_2=X^aa_p be the two successive singleton queries: the first new vertex beats only p, and the second beats only the first. Set V to be the 3-by-3 matrix whose columns are the endpoint states of X_0,X_1,X_2 at q. The same C/D matrix acts on every column. This is inherited query commutation, because every new guard vertex beats every query vertex. Recover the observer from h_i=(1,1,0)V_i by

    v(X,p)=(h_0,(h_1−h_0)/2,(h_2−h_1)/2).                 (2)

Our subject is the information retained by M_w, rather than just its value on the canonical ray (1,2,0).

## 2. Exactly two primitive rank-losing gaps

Let P be the plane d=0. After a D move the state belongs to P. Between consecutive D letters, a run of k C letters therefore induces on P the two-dimensional map

    A_k = first two rows and columns of D C^k.

Define the bilateral Tribonacci sequence by

    T(0)=0, T(1)=T(2)=1,
    T(n+3)=T(n+2)+T(n+1)+T(n),  n in Z.

**Exterior-power identity.** For every k≥0,

    det A_k = −2 T(−k−3).                                  (3)

**Proof.** In the oriented-pair basis (01,02,12), the second exterior power of C is

    E = [[−1,0,1],[−1,0,0],[−1,−1,0]].

Direct multiplication gives E³+E²+E−I=0. The determinant det A_k is a fixed linear functional of E^k, so it obeys

    delta_(k+3)=−delta_(k+2)−delta_(k+1)+delta_k.

Its initial values are 2,0,−4. Reversing the defining Tribonacci recurrence shows that −2T(−k−3) obeys the same recurrence with these initial values. This proves (3). ∎

We now use the **classical integral-zero theorem**

    T(n)=0 exactly for n in {0,−1,−4,−17}.                  (4)

The indexing and statement are given explicitly on page 1333 of Bilu, Luca, Nieuwveld, Ouaknine and Worrell, [*On the p-adic zeros of the Tribonacci sequence*](https://people.mpi-sws.org/~joel/publications/tribonacci22.pdf), Mathematics of Computation 93 (2024), 1333–1353. They cite Mignotte and Tzanakis, *Arithmetical study of recurrence sequences*, Acta Arithmetica 57 (1991), Example 2 on page 360, [DOI](https://doi.org/10.4064/aa-57-4-357-364). The infinite exclusion in (4) is an imported theorem, not a conclusion of our finite search or a new Tribonacci result.

**Rank-language theorem.** Over Q or R, every finite C/D word satisfies:

1. rank M_w=3 if w has no D;
2. rank M_w=1 if w contains DCD or D C^14 D as a contiguous subword;
3. rank M_w=2 otherwise.

**Proof.** C is invertible with determinant one. The first D drops rank to two and has image P. Until another rank drop, this image is transported through the invertible C runs and then through the maps A_k. By (3)–(4), A_k is singular exactly for k=1 or k=14. It cannot have rank zero: C^k maps a two-dimensional plane injectively, while ker D is one-dimensional. Thus the first such gap drops the total rank to one. A later C or D cannot raise rank. Nor can it kill this image: C and D are nonnegative with no zero column, so neither annihilates a nonzero nonnegative vector. Their product sends the strictly positive vector (1,1,1) to a nonzero vector at every stage. So the rank remains one. Initial and final C runs do not change it. ∎

This identifies a second reset that cannot be detected by looking only for the short DCD pattern. It also says exactly when the universal transfer loses two independent endpoint directions. A particular seed or query matrix can already have smaller rank; the theorem classifies the operator, not the rank of every particular input.

For completeness, the number r_n of rank-two words has generating function

    sum_(n≥0) r_n z^n
      = z / [(1−z)(1−2z+z²−z³+z^15−z^16)].               (5)

Indeed a rank-two word has arbitrary initial and final C runs and at least one D, with every internal C-run length avoiding 1 and 14. Summing these runs gives

    (1−z)^−2 z / [1−z((1−z)^−1−z−z^14)].

Equation (5) is a count of operator words, not a count of distinct tournaments or distinct matrices.

## 3. The long reset and an exact source return

The two reset matrices have particularly simple factorizations:

    M_DCD       = (1,2,0)^T (1,2,1),
    M_(D C^14 D)= (47,112,0)^T (47,112,47).               (6)

For the second identity, an exact intermediate certificate is

    C^14 = [[1431,1705,927],
            [2632,3136,1705],
            [778,927,504]].

Multiplying by D on both sides gives

    [[2209,5264,2209],
     [5264,12544,5264],
     [0,0,0]].

Consequently the 16-letter word D C^14 D sends every guarded source to

    (s,u,d) ↦ (47s+112u+47d)(47,112,0),                  (7)

with Hamilton count 159(47s+112u+47d). It is a literal induced extension and preserves strongness, not just a formal map.

Write R=M_(D C^14 D). Then

    R²=14753 R,     14753=47²+112².                      (8)

Thus the same word is a reusable source return with multiplier 14753 on the new ray (47,112,0). Starting with the directed triangle, whose state is (1,2,0), one reset gives

    (s,u,d)=(12737,30352,0)=271(47,112,0),
    H=43089=159·271,      order=19.

Repeating the reset m further times gives actual strong sources

    H=43089·14753^m,       order=19+16m,       m≥0.        (9)

This is a sparse multiplicative family, not a new all-count theorem or an efficient minimum-order construction. S164 already classified canonical-ray return words; (7) identifies a different universal reset ray and explains its arithmetic origin.

## 4. Marked-observer locking and a concrete path-order effect

**Observer-locking theorem.** If a C/D word contains either reset pattern, then every subsequent C/D extension, using the same physical old observer p, preserves the ratios x/h and y/h of the observer (h,x,y) in (2).

**Proof.** After the reset, the query matrix has the form b lambda^T, with one shared nonzero endpoint column b and three scalar channel weights lambda_i. Every future word replaces b by M_w b and leaves lambda unchanged. Its three query counts are therefore all multiplied by the same positive scalar. Equation (2) gives the same scaling for (h,x,y). ∎

In particular, after a reset no further word in this language can realize

    v_new = a v_old + b(2,−1,0),       b≠0,              (10)

for any scalar a. A physical observer has h>0 and x,y≥0, so it is not parallel to (2,−1,0). Here nonnegativity can be proved without a transfer calculation: any old directed Hamilton path admits insertion of an added tournament vertex somewhere along the path. Choose a deterministic valid insertion; deleting the added vertex recovers the old path, giving an injection. Thus each successive query can only increase H. Since the locking theorem requires v_new to be parallel to v_old, (10) is impossible.

This is a genuine obstruction for post-reset **marked** return searches. It does not prevent an offset before the reset, moving or replacing the observer, adding an operation outside C/D, or using a larger interface.

A small literal example separates path order from scalar output. Start with rows (2,4,1), guard q=0 and observer p=1. Its three query columns are

    V=[[1,2,3],[2,3,6],[0,1,2]],       det V=−2.

The words DCD and CDD each use one C and two D letters. Both produce strong six-vertex tournaments with H=15, but their observers differ:

    DCD: rows (26,4,1,6,46,15), guard 4, observer (15,6,12);
    CDD: rows (10,4,1,54,7,23), guard 3, observer (15,7,12).

The first operator has rank one; the second has rank two. Thus even identical letter counts, vertex cost, and unmarked Hamilton count do not determine the marked response. This is noncommutative path dependence, not a claimed fundamental-group action: D is singular and irreversible in this observer.

For the long reset, the same physical seed and observer give endpoint query columns

    (12737,30352,0), (22419,53424,0), (42629,101584,0),

and observer (43089,16377,34185). The independent C++ certificate counts these directly on the literal tournaments of orders 19,20,21.

## 5. The full profile remembers what H forgets

The [S164 full-slot refinement](https://github.com/eliottcassidy2000/math/blob/934faa959d0fa95aab888be9de05ecd0bf18e7e6/05-knowledge/results/session164_emergent_20261103/profile_algebras/POINTED_PROFILE_ALGEBRA.md) replaces (1), at an integer number r of labeled path slots, by

    C_r=[[r−1,r,0],[1,1,1],[1,0,0]],
    D_r=[[r,0,1],[0,r+1,0],[0,0,r−1]].                   (11)

These act on the three endpoint-restricted path-cover polynomials S(r),U(r),d(r). Their determinants are

    det C_r=r,             det D_r=r(r²−1).

Hence for a word with length n and j occurrences of D,

    det M_w(r)=r^n(r²−1)^j.                              (12)

Every word is therefore invertible over Q at every integer r≥2 and injective on polynomial triples over Q[r]. The reset does not erase the full pointed profile. Its rank loss is specific to Hamilton evaluation r=1. The same determinant for all words with equal letter counts also shows why determinant alone misses the noncommutative rank distinction at r=1.

There is a sharper derivative statement. Fix the two-D gap word D C^k D and put t=r−1. Its determinant has a zero of order exactly two at t=0. If k=1 or 14, its value at t=0 has rank one, so all its 2-by-2 minors vanish at t=0. The adjugate formula then shows that its inverse has at most a simple pole in t. Therefore the output value and first derivative at r=1 determine the **input value** at r=1, even when the input polynomial triple itself is not constant.

For every other k≥0, the specialization has rank two, so at least one 2-by-2 minor is nonzero at t=0. The inverse then has a genuine double pole. Algebraically the local elementary-divisor orders are

    (0,1,1) for k in {1,14};      (0,0,2) otherwise.      (13)

For the latter case, two derivatives suffice in the unrestricted linear polynomial interface, while one derivative cannot suffice for every polynomial input. To see necessity, take a column with nonzero leading double-pole coefficient of the inverse, multiply that column by t², and truncate the resulting regular input through degree one. It has a nonzero input value but output vanishing to order at least two. Thus its zero and first output jets are indistinguishable from zero. This is an algebraic reconstruction statement, not a claim that every such formal input difference is realized by two tournaments.

Equations (12)–(13) expose the relevant boundary: the Hamilton observer can become permanently one-dimensional while the richer cover-profile interface retains the missing directions just above the specialization.

## 6. Reproduction and scope of verification

Run `python3 check.py` and `python3 -O check.py` in this directory. The standard-library checker has explicit failure tests that survive optimization. It verifies:

- the exterior-power matrix and recurrence, the gap identity through k=200, and both exact reset factorizations;
- the rank classification by direct integer minors for all 524,287 words of lengths zero through eighteen, including the new length-sixteen pattern;
- the rank-two generating function against that exhaustive census;
- 378 full-slot determinant controls and 41 gap-jet checks, including the two different local elementary-divisor patterns;
- 63 literal strong guard itineraries through depth five, each with all three independently counted physical query states;
- observer-ray locking for every suffix through length six, and the long-reset transfer and graph construction.

The all-depth exclusion of other singular gaps depends on the cited classical zero theorem. The finite tests are controls, not a replacement for that result.

For a separate literal computation of the large example, run

    g++ -O2 -std=c++17 check_literal.cpp -o check_literal
    ./check_literal

It uses subset/endpoint dynamic programming directly from adjacency rows, with no C/D transfer formulas, and verifies all three endpoint query states at orders 19,20,21 and the necessary deleted-guard counts. Every counted graph has at most 21 vertices; all counts fit in unsigned 64-bit integers because 21! exceeds that limit but the program's actual subset counts in this particular sparse family are small (the checked outputs are recorded in `literal.json`). The implementation uses an overflow check before each addition, so this is verified rather than assumed.

The packet contributes the reset-language reduction, the second physical reset, its marked locking consequence, and the profile-derivative distinction. The physical C/D calculus, full-slot matrices, and integral Tribonacci zero theorem are inherited. No homotopy group, autonomous tournament dynamics, or general odd-spectrum classification is claimed.
