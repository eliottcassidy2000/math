# Fresh-block apex migration: a closed interface and an all-length strip obstruction

**Status:** PROVED algebraic constructor/finite-state reduction; FINITE-EXACT complete restricted-language exclusions and explicit repairs. The unrestricted claim that every positive odd H except 7 and 21 occurs remains **OPEN**. No literature-priority claim is made.

## Inheritance and scope

Baseline: `04cdf182d9e4543800b421a8893232ccf69b2399`.
The closest inherited mechanisms are:

* `05-knowledge/results/session124_digital_constructor_inverse_20261031/PROOF.md`: actual marked state, singleton maps, source ancestry, frozen fibers, and the distinction between a complete inverse and universal acceptance.
* `05-knowledge/results/session124_digital_constructor_inverse_20261031/global_reduction/README.md`: unrestricted reduction to non-3-divisible targets, not equivalence with any particular grammar.
* `05-knowledge/results/square_sum_spectrum_overnight_20261031b/README.md`, sections 3 and 5: ordinal block action and the positive cut/interior coordinates. The finite scalar-atlas obstruction expressly leaves lawful re-marking outside its scope.
* `05-knowledge/results/square_sum_spectrum_overnight_20261031h/README.md`: bounded source automaton, actual source entry, and the warning that numerical four-state equality does not justify arbitrary re-marking.
* `01-canon/theorems/THM-1370-h-spectrum-omits-7-21-all-n.md`: inherited unrestricted exclusions. Neither exclusion is new here.

The live concepts are source-preserving migration, positive cut/interior coordinates, finite capped semirings, and full boundary response rather than only scalar H. The least-used inherited sidecar is the **position of the marked vertex**: it now determines which fresh edges can return to the old parent. The source graph, marked cut, local recipe, and vertex clock are retained. A scalar count alone is never an input constructor.

## 1. State and the elementary migration

Let A be a nonempty tournament and x an external marked vertex with fixed cut S=N+(x) intersect A. Use

    h=H(A), F=H(A+x),
    q=# ordered path covers (P,Q) of A with Q nonempty and first(Q) in S,
    g=# Hamiltonian paths of A ending in S.

P may be empty. The inherited positive coordinates are

    C=(F-h)/2, L=(q-F+h-g)/2, t=h-g,
    (h,F,q,g)=(g+t, 2C+g+t, 2C+2L+g, g).

Append z with A->z->x, and **re-mark z**, leaving parent A+x and cut {x}. This gives

    R1(h,F,q,g)=(F, F+q+g, q+h, h-g),
    R1(C,L,g,t)=(C+L+g, 0, t, 2C+g).

The total count is the inherited zero-singleton count. The new q counts covers whose second path starts x: a nonempty suffix after x is exactly an old q cover; the empty suffix Q=(x) contributes h. The new g counts old full paths ending x, namely h-g. Dropping the empty suffix is already false for the one-vertex empty-cut source.

This is not arbitrary re-marking and not a claimed converse symmetry. The new mark has no out-neighbor in old A. That restriction is precisely what allows closure.

## 2. General fresh-block migration theorem

Take any finite local tournament G on old port x and a nonempty set of fresh vertices. Choose a fresh vertex z as the new mark. Orient **every edge from A to every fresh vertex**; keep the old A-x cut and all edges of G. Put

    B=V(G) minus {x,z}, R=N+_G(z) intersect (B union {x}).

Then the four coordinates of the new parent A+x+B with marked vertex z are an **integer-linear function of (h,F,q,g)**. The following finite local formulas explicitly compile that matrix. No claim that every integer vector is an actual source is used.

For a subset U of fresh vertices, let H(U) count paths of G[U], with H(empty)=1. Let s(U) count full paths of U starting at an out-neighbor of x. Let p(U) count ordered two-path covers (P,Q) of U with P nonempty and last(P)->x; Q may be empty. Define

    Phi_U = H(U)F + p(U)q + (H(U)-s(U))g
            + (H(U+x)-H(U))h                 (U nonempty),
    Phi_empty = F.

This is exactly the inherited ordinal-block action: Phi_U=H(A+x+U). Thus

    h_new=Phi_B, F_new=Phi_(B union {z}).

For g_new, if B is empty, use g_new=h-g when z->x, and zero otherwise. If B is nonempty, define local integers:

* e: full paths of B ending in R;
* se: full paths of B starting at an out-neighbor of x and ending in R;
* pe: ordered covers (P,Q) of B, both nonempty, with last(P)->x and last(Q) in R;
* k: full paths of B+x ending in R.

Then

    g_new = e F + pe q + (e-se)g + (k-e)h.       (1)

For q_new, sum over every partition B=U disjoint-union V. Let

* s_R(V): full paths of V starting in R (zero for V empty);
* b_R(V): full paths of V+x starting in R;
* a_R(V): ordered covers (K,J) of V+x with K nonempty, K ending x and starting in R, and J arbitrary, possibly empty.

Then

    q_new = sum_(U subset B) [ s_R(V) Phi_U
                             + H(U) a_R(V) q
                             + H(U) b_R(V) h ].             (2)

These formulas are implemented by local permutation/partition enumeration in `migration.py::compiler`. The program's pair-code convention is lexicographic pairs i<j, bit 1 meaning i->j; old port is label 0. `extend` reconstructs every edge and the new marked cut.

### Proof of the endpoint and cover formulas

Every fresh vertex loses to A. Consequently a path can return from fresh vertices to A only through the single old port x, used once.

For (1), separate paths with no fresh vertex before x from paths with a nonempty fresh segment immediately before x. In the first class the old A+x path is followed by a B path ending in R. An old path ending x is admissible exactly when the B path starts at an out-neighbor of x; this gives eF+(e-se)g+(se-e)h. In the second class, a nonempty old A suffix after x gives pe*q. An empty A suffix leaves a local B+x path ending in R; subtract the se paths with x first to avoid the first class. This contributes (k-se)h. Adding gives (1).

For (2), delete the new root and consider a cover (P,Q) of the new parent whose Q starts in R. If x lies in P, Q contains **only** fresh vertices: without x it cannot return to A. For its vertex set V there are s_R(V) choices; P is an arbitrary full path of A+x+U, giving the first term. If x lies in Q, all A vertices in Q form one suffix after x, possibly followed by a fresh suffix. When that A suffix is nonempty, the old A pieces are exactly an old q cover; the local Q pieces form a cover counted by a_R(V), and the fresh P part has H(U) choices. If the A suffix is empty, all A lies in P, contributing h, and Q is a local path counted by b_R(V). These cases are disjoint, exhaustive, and invertible. Empty P and empty local suffixes are explicitly allowed.

The restriction on fresh-to-A edges is essential. If a fresh vertex could return directly to A, a cover could have more old-A pieces; these four coordinates would no longer be justified as complete.

No universal nonnegativity theorem for the matrix in (C,L,g,t) is asserted here. For the explicitly used finite bank below, nonnegative integer entries are checked exhaustively. The strip matrices have separate displayed nonnegative formulas.

## 3. Arbitrarily long zero strips, with the correct cover coefficient

Append a transitive fresh strip z1->...->zr with A->zi->x for every i, and mark its **first** vertex z1. For r>=2 put

    a=2^(r-1), b=2*3^(r-1)-2^r+1.

The exact result is

    h'=F+(a-1)q+g,
    F'=F+(2a-1)q+g,
    q'=a h+(a-1)F+b q+(a-2)g,
    g'=F+(a-2)q+g.                              (3)

For r=1 use R1 above. In positive coordinates (3) becomes

    C'=a C+a L+(a/2)g,
    L'=(b-a)C+(b-2a+2)L+((a+b-3)/2)g+(a-1)t,
    g'=2(a-1)C+2(a-2)L+a g+t,
    t'=2C+2L+g.                                 (4)

Every coefficient in (4) is a nonnegative integer for r>=2.

### Derivation

The old singleton formula gives h' and F'. A full path of the new parent ending in A has all r-1 fresh vertices consecutively before x, and is specified exactly by an old q cover. Hence g'=h'-q. For q', apply (2) with m=r-1 fresh vertices, all selected, all beating x, and all mutually transitive. The coefficient of q from x in P is

    sum_(U proper subset B) (2^|U|-1)
    =3^m-2^(m+1)+1.

The coefficient from x in Q is sum_(U subset B) 2^|V|=3^m. Their sum is b. The F,h,g coefficients are respectively 2^m-1, 2^m, 2^m-2. This proves (3); conjugation gives (4).

**Retained pre-release correction.** The tempting coefficient 3^m counts only the second allocation of x, and is false from r=3. For actual old state (h,F,q,g)=(1,1,1,1), the new state is (5,9,20,4), not a state with q=18. The missing covers allow fresh vertices after the old-A suffix. Earlier exploratory claims based on coefficient 3^m are discarded and are not inputs to any result here. The corrected formula agrees with the general compiler through r=6 and with the independent literal consumer; the proof above is for all r.

## 4. Exact all-length decision by saturation

For M>=1 let cap(n)=min(n,M+1), coordinatewise on vectors and entrywise on nonnegative integer matrices. For every such matrix A,

    cap(A cap(v))=cap(Av).

This follows from cap(a+b)=cap(cap(a)+cap(b)) and cap(ab)=cap(cap(a)cap(b)) for nonnegative integers. Zero coefficients may erase overflow; no illegal monotone-q pruning is used. Composition preserves the identity. Thus breadth-first search of capped states from an **actual** seed exactly decides every output F<=M in a fixed nonnegative matrix grammar, with no word-length cutoff. Every stored edge retains an actual recipe, which can be replayed without saturation. Output is F=2C+g+t; L overflow is retained, not thrown away.

For the all-strip grammar use S0,S1, Rr for **every r>=1**, and reset Z. The inherited singleton maps and lawful reset are

    S0=(2C+L+g,L,0,g+t),
    S1=(C,C+2L+g,g+t,0),
    Z=(0,0,0,2C+g+t).

Z appends a fresh universal sink and re-marks it, so it must not be confused with a fabricated zero state. The start is (0,0,0,1), the one-vertex parent and empty cut.

The apparently infinite alphabet also has an exact finite quotient. Let

    r0=ceil(log2(M+1))+2.

For every r>=r0, all growing positive entries in (4) are at least M+1; all other entries are the same fixed zeros, ones, and twos. Thus cap(Rr)=cap(Rr0) for every r>=r0. The finite alphabet S0,S1,Z,R1,...,Rr0 therefore covers **all strip lengths and all word lengths** for targets at most M. This is a mathematical stabilization argument, not a sampled eventual-period claim.

Each operation adds vertices without changing old induced edges, and the formulas give nondecreasing F directly. The implementation may therefore discard states already having F>M. Alternatively retaining every capped state gives the same target language.

## 5. A genuine obstruction and its actual mixed-cut repair

At M=173, r0=10. The exact 13-generator closure has **5,373** capped states and misses exactly

    7, 21, 67, 173

among positive odd counts <=173. Thus **67 and 173 are impossible at every depth and every strip length in this specified grammar**, even with reset. This is not an unrestricted exclusion. It closes a concrete enlargement of the pure singleton route, rather than merely testing its first few levels.

As an alternative authenticated grammar, allow any fresh **transitive** block of at most three vertices, with arbitrary edges between old x and that block, and re-mark any fresh vertex. Relabeling puts the new mark last; the compiler enumerates all local tournaments of orders 2..4 whose fresh induced subtournament is transitive. This fixed bank is not claimed to contain every arbitrary-length strip operation. Together with S0,S1 there are **32 distinct nonnegative integer positive-coordinate matrices**. At the same output cap, complete closure has **6,234** states and misses only 7 and 21. This small comparison is a mechanism control, not an improvement of the inherited 100-million finite coverage and not a universal acceptance theorem.

Two actual repairs are retained, counting raw Hamiltonian paths, not paths modulo automorphisms:

* 67: start with S1, then local recipes (n,code,z)=(4,2,3),(3,0,2). The final parent has masks (126,120,82,68,72,28,32), cut 96, and state (31,67,73,13). Its full tournament has 8 vertices.
* 173: local recipes (3,2,2),(4,4,3),(4,2,3). The states are (1,3,2,0), (9,19,26,6), (89,173,325,63). The final parent masks are (510,504,498,452,392,408,304,64,128), cut 448. Its full tournament has 10 vertices.

`verify.py` rebuilds all adjacency matrices, counts their paths/ordered covers independently of the transfer coefficients, and checks every intermediate state. The independent auditor additionally checks the compiler against separately constructed graphs and literal covers. A mixed old-x cut on a fresh transitive block is enough to leave the zero-strip language; a nontransitive fresh block is not required for these repairs.

## 6. Reproduction and remaining obligation

From any working directory, with these files together:

    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/verify.py
    PYTHONDONTWRITEBYTECODE=1 python3 -O /path/to/verify.py

The retained normal and optimized JSON outputs are identical. A fresh-directory replay from working directory `/` in both modes also gives the identical receipts (`relocated.normal.json`, `relocated.optimized.json`). No asserts, external data files, third-party Python modules, random seed, or arbitrary depth bound is used. The complete finite universes and actual witnesses appear in those receipts. The independent audit scripts/receipts are separate from the producer.

The new source-preserving four-state migration interface is a lawful candidate for a uniform construction. What is still missing is a theorem that the fixed mixed-chain bank, a larger authenticated bank, or some controlled family of local blocks accepts **every** eligible target. A complete finite-state decision at each fixed M does not imply that statement. No global attainment, density, polynomial-time, or unbounded consecutive interval theorem is claimed.
