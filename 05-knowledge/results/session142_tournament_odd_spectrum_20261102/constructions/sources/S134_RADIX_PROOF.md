# Fresh-block affine rays: a finite obstruction and an actual radix-80 survivor

**Status:** FINITE-EXACT no-two-nonzero-digit result for words of at most four operations; PROVED explicit infinite family. The unrestricted Hamiltonian-path spectrum claim remains **OPEN**. No literature-priority claim is made.

## 1. Inheritance, scope, and the sidecar

This appendix uses the authenticated four-coordinate fresh-block compiler in the sibling Session134 constructor packet, `constructor/PROOF.md`, Sections 1–2 and 5. It does not use an arbitrary numerical vector as a source tournament. The old port, the marked cut, and the actual graph recipe are retained.

The state is `v=(C,L,g,t)`, with full Hamiltonian-path count

    c v = 2C+g+t,             c=(2,0,1,1).

The 32-generator bank consists of the inherited singleton operations S0,S1 and every fresh transitive block of at most three vertices, with arbitrary edges to the old port, followed by re-marking a fresh vertex. By relabeling fresh vertices the new mark may be put last. Local recipes are `[n,code,n-1]`, using lexicographic pair bits, bit 1 meaning the smaller label points to the larger. `bank.json` freezes the distinct matrices and a recipe for each. `migration.py` is an unchanged copy of the audited constructor compiler; `verify.py` regenerates and compares the complete bank.

**Word length here counts operations, not added vertices.** One operation can add up to three vertices. This is neither the all-strip grammar nor the language of every possible fresh block. It is not the earlier fixed-two-ear model. The earlier eigenvalue-{1,2} probes only motivated this test and are not inputs to the result.

The local concept board is: common affine rays; the integer spectral gate; preservation of actual source/cut; lower-dimensional invariant faces; and lattice-phase information lost by scalar output. The positive survivor below shows why genuinely migrating words must not be discarded merely because elementary migration lacks eigenvalue 1.

## 2. Why nonzero affine digits require two eigenvalues

Suppose a fixed affine profile

    v(H)=a H+b

has output `c v(H)=H` and is componentwise nonnegative for an unbounded set of positive H. Suppose a word matrix M implements an affine digit

    M v(H)=v(BH+d)

for infinitely many such H, with fixed integer B>=2. Comparing coefficients gives

    c a=1, c b=0, a>=0,
    M a=B a, (M-I)b=d a.                         (1)

The nonnegativity of a follows from unbounded H, not from any sign assumption on b. **The offset b may be signed.**

If d is nonzero, M must also have eigenvalue 1. Otherwise M-I is invertible, and its action on a gives

    b=d a/(B-1),

contradicting c b=0. This is only a necessary condition: neither a positive eigenvector nor (1) proves that a profile is realizable by tournaments.

If d=0, the eigenvalue-1 gate need not hold, for example when b=0. Zero-digit maps are not discarded from the underlying grammar and are outside this necessary gate.

## 3. Complete finite spectral certificate

Enumerating **all** words at each indicated length and identifying equal matrices gives:

| Operations | Distinct matrices | Matrices having eigenvalue 1 |
|---|---:|---:|
| 1 | 32 | 7 |
| 2 | 690 | 15 |
| 3 | 14,396 | 49 |
| 4 | 297,239 | 143 |

The union contains 166 distinct matrices with eigenvalue 1. For this union, the verifier computes the characteristic polynomial of every matrix with eigenvalue 1, removes its zero roots, and uses the monic rational-root theorem to test **every integer eigenvalue B>=2**. It checks that every such B-eigenspace is one-dimensional; no higher-dimensional positive cone is silently omitted. It normalizes the eigenvector by c a=1 and retains exactly those with a>=0.

There are 88 integer-radix eigenspaces in this union; all are one-dimensional. Four have output zero and none has a sign-changing output-normalized vector. There are exactly **84 distinct retained matrices**. The possible integer radices are

    2, 4, 8, 16, 80, 216, 508, 676, 728, 872,
    1544, 1600, 1672, 3200, 6400.

The smallest retained radix from a matrix not represented by a pure-singleton word is 80, and Section 4 supplies an actual source on such a ray. This minimum is scoped to this finite universe.

Most importantly, the map

    M -> (B, its c-normalized nonnegative B-eigenvector)

is **injective** on this retained set. The complete list, exact rational rays, representative words, matrices, and characteristic coefficients are in `normal.json`. Different words sometimes give the same matrix; deduplication occurs across all four lengths, not just within each layer.

**Finite obstruction.** No fixed affine profile can support two distinct nonzero digits at a common integer radix B>=2 using words of at most four operations from this bank.

Indeed, (1) puts both digit matrices in the retained set with identical `(B,a)`, so injectivity makes the matrices equal. Their two equations `(M-I)b=d a` then imply equal digits since c a=1.

In particular, for B=8 the only matrices with eigenvalues 1 and 8 are the eight pure-singleton words of length three. The all-S1 word has output-zero 8-eigenvector and is not retained. The other seven have pairwise different normalized rays. No length-four matrix adds a radix-8 candidate.

A residue-complete odd-output fixed-radix scheme would require B distinct odd residue classes modulo 2B. For even B all offsets are odd and nonzero. For odd B>=3 at most one offset can be zero, leaving at least two nonzero offsets. Therefore the finite obstruction excludes such a complete fixed-profile bank within this operation-length budget. It does **not** exclude a bank whose profiles change with a finite phase, higher-dimensional invariant sets, longer words, variable radices, or more general actual constructions.

## 4. A positive migrating affine ray with an actual seed

The genuine migration word

    W = [(3,0,2), (4,29,3), (3,5,2)]

has matrix

    [64 44 46 18]
    [ 0  0  0  0]
    [ 0  0  0  0]
    [56 38 41 17].

It preserves the face `L=g=0`; on that face

    (C,t) -> (64C+18t, 56C+17t).

The left invariant and output law are

    I=8C-9t,
    I'=I,                F'=80F+3I, F=2C+t.       (2)

This is checked directly: `8C'-9t'=8C-9t`, while `2C'+t'=184C+53t=80(2C+t)+3(8C-9t)`.

Start from the **actual** one-vertex parent with empty marked cut, state `(0,0,0,1)`. The operation theorem supplies a valid actual marked tournament at every iteration. Here I=-9, so (2) gives the all-k family

    F_k = (52*80^k+27)/79,
    C_k = 18*(80^k-1)/79,
    L_k = g_k = 0,
    t_k = (16*80^k+63)/79,                  k>=0. (3)

The formula is integral because 80 is 1 modulo 79; all counts are positive odd. Each iteration adds seven vertices, so the full tournament has `7k+2` vertices. For k=0 the full graph is the transitive two-vertex tournament.

The first output is 53 and the next is 4,213. For the first iteration, `verify.py` reconstructs every adjacency matrix and independently counts Hamiltonian paths and the ordered covers at every intermediate stage. The final full tournament has nine vertices, native state `(h,F,q,g)=(17,53,36,0)`, and positive state `(18,0,0,17)`; the complete parent masks and cut are retained in the receipt. Equation (2), not extrapolation from these samples, proves all of (3).

As an affine profile this orbit lies on

    a=(9/26,0,0,4/13), b=(-9/26,0,0,9/13),
    v(F)=aF+b,             F -> 80F-27.

Its integer profile requires `F=1 mod 26`. The actual W-orbit is a specific sparse subset of this congruence class; the statement does not make every integral point on the ray realizable. The denominator 26 is a phase/lattice sidecar any attempted extension must retain.

## 5. Failure anatomy and next boundary

The first failed implication was “fresh migration gives more scalar values, hence it may give interchangeable scalar digits on one invariant profile.” The finite certificate shows that the short-word integer spectral rays are rigid: at a fixed B, each allowed ray identifies one matrix. However migration **can** preserve a nontrivial affine ray, as (2) demonstrates; an all-length claim that productive migration destroys eigenvalue 1 is false.

The strongest survivor is an actual nonbinary infinite family with a signed affine offset and an exact modular phase. A repaired positive target would use a finite collection of lawful profiles with residue-dependent transitions, or an authenticated higher-dimensional region, rather than assume one common ray. Neither target is settled here. No depth-five computation or finite target-prefix search is being promoted to an all-length assertion.

## 6. Reproduction and controls

With the files together, from any working directory:

    PYTHONDONTWRITEBYTECODE=1 python3 /path/to/verify.py
    PYTHONDONTWRITEBYTECODE=1 python3 -O /path/to/verify.py

The verifier uses only the Python standard library, exact integers, and exact `Fraction` arithmetic. It uses no asserts, floating point, random seed, external data, or assumed target-state realizability. The characteristic-polynomial calculation and determinant calculation are separate exact algorithms, and eigenvectors are checked by multiplication. The identity matrix is a positive control for repeated eigenspaces; a polynomial with a repeated root and integer root 80 checks root extraction. The actual W witness is counted independently of the transfer formulas.

Normal, optimized, relocated-normal, and relocated-optimized receipts are byte-identical. The relocated runs copy only `verify.py`, `migration.py`, and `bank.json` into a fresh directory and run with working directory `/`; no exploratory file or source checkout is needed.

The supporting migration compiler is reproduced verbatim solely for portable bank regeneration; its general correctness theorem and independent audit are in the sibling constructor packet. The finite computation above consumes that theorem without claiming to reprove it from finitely many graphs.
