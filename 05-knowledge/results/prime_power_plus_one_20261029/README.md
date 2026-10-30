# Prime powers plus one: recursive tournaments, norm clocks, paired lattices

**PROVED scoped deductions / FINITE-EXACT controls / CITED classical
inputs. Same-agent proof review; no independent referee or claim of
literature priority.** Session continued on2026-10-30 from the existing
isolated worktree; the directory retains the preceding session's date.

The organizing distinction is between three constructions of the same
integer p^k+1: adding a vertex to a recursive tournament, parametrizing
a finite-field norm-one group, and taking four times a tournament
lattice norm. Their operations are explicit, but their symmetry groups
are different. A common size does not identify the objects.

## Recursive families and complete extension recurrences

For the carousel tournament R_p on any odd prime p, the parity cut
on R_p^[k] produces a rigid, unmarked tournament U_k with

    |U_k|=p^k+1,
    cut weight=[p^k-(p-2)^k]/2,
    largest transitive set=((p+1)/2)^k+1.

| base | size of U_k | cut imbalance | largest transitive set |
|---|---|---|---|
| R3 | 3^k+1 | 1 | 2^k+1 |
| R5 | 5^k+1 | 3^k | 3^k+1 |
| R7 | 7^k+1 | 5^k | 4^k+1 |
| local11 inside Paley23 | 11^k+1 | 9^k | 4^k+1 |
| local13 inside Paley27 | 13^k+1 | 9^k | 4^k+1 |
| local29 inside Paley59 | 29^k+1 | 27^k | 8^k+1 |

The last three bases have full groups of orders11,39,29, respectively;
their seed cuts and transitive maxima are FINITE-EXACT. The all-depth
conclusions follow from the proved substitution and attachment rules.
R7 is the carousel, not the different Paley7 tournament.

For prime carousel bases, all binary-coloring orbits F_k and
trivial-stabilizer orbits R_k satisfy the complete recurrences

    F_0=R_0=2,
    F_(k+1)=[F_k^p+(p-1)F_k]/p,
    R_(k+1)=[R_k^p-R_k]/p.

There are F_k-1 unmarked one-vertex children. At depth k>=2 exactly
R_k are rigid. The ternary case is the unique prime base for which the
binary regular palette stays fixed at2. At parent orders9,25,49 the
rigid child counts are2,1554,87460002. These are extension-family counts,
not the whole A000568 sequence. Totient/Möbius necklace formulas extend
the construction to composite odd carousel bases.

The mixed15-vertex cores R3[R5] and R5[R3] have the same parity weight
and transitive maximum, yet group orders375 and1215. Their rigid
16-vertex children remain nonisomorphic. The ordered substitution tree
must therefore survive any proposed scalar compression.

## A sharp guard for several attachments

For every doubly regular tournament H on q>=7 vertices:

    minimum deletions making H decomposable = (q-1)/2,
    number of minimum deletion sets = q(q+3)/2.

Every induced subtournament larger than(q+1)/2 is prime. At the boundary
there is exactly one proper nontrivial module, either a pair or a whole
regular neighborhood under a source/over a sink. The proof uses the exact
Gram identity and parity of the row sums; it applies beyond Paley fields.

Retaining the actual prime obstruction types sharpens the size bound:
every R5^[k] tolerates any two new vertices without creating Paley7,
while three can suffice. Every Paley7^[k] tolerates any five without
creating Paley11, while six can suffice. Both thresholds are exact for
every k>=1. These give concrete robust forbidden-pattern families, not
a proof of the general tournament Erdős--Hajnal conjecture.

## Two arithmetic structures attached to q+1

For q=p^k odd, the norm-one group in F_(q^2) has order q+1 and is
parametrized by P1(F_q). Its multiplication has the exact chart

    t star s=(ts+delta)/(t+s),

with a retained nonsquare delta and infinity as identity. A primitive
clock has Frobenius degree2k; its reciprocal trace has degree k.
Trace loses a branch that future multiplication can detect. The full
clock has negative determinant character, while its square preserves
the switching class. The quadratic relation is a tournament only for
q=3 mod4; it is symmetric for every5^k.

For any doubly regular q-tournament, the inherited integral operator
theta=(I+S)/2 on the root lattice has the complete Smith form

    coker(theta) = (Z/C)^r,  C=(q+1)/4, r=(q-1)/2.

The complementary cokernels of theta and1-theta have an explicit
perfect pairing; their direct sum has a perfect alternating pairing
and square order C^(q-1). Examples include (Z/4)^7 at the non-field
order15 and (Z/7)^13 at the field order27. This is an actual paired
finite group. No identification with a Selmer or Sha group has been
proved, and Sha finiteness is not a consequence.

For C=p prime, its Lagrangian-subspace counts satisfy the classical
recurrence L_j(p)=(p^j+1)L_(j-1)(p). The tournament lattices at
q=11,19,27 realize p=3,5,7 respectively. Thus all three requested
families occur as literal growth factors in one counting problem.

For p=3 mod4 and odd k>=3, a primitive divisor ell of p^k+1 forces
a new paired torsion factor absent at every smaller exponent. It is
also absent from the supplied affine--Frobenius tournament group.
Examples are ell=7 at3^3,61 at3^5, and43 at7^3. This is a precise
connection between arithmetic growth and symmetry, with the group
actions kept distinct.

## Evidence, sources, and reproduction

- [GROWTH.md](GROWTH.md): all-size proofs T1--T3 and D1--D3, including the sharp deletion classification.
- [ARITHMETIC.md](ARITHMETIC.md): integer formulas, norm clocks, exact phase, lattice and pairing proofs A1--A4 and L1--L3.
- [RECOVERY.md](RECOVERY.md): inheritance, niche sources, literature boundary, typed connections.
- [certificate.json](certificate.json): field charts, all finite results, prime obstruction banks and embedding maps.
- [HANDOFF.md](HANDOFF.md): remaining proof obligations and the next decisive experiments.

From the repository root:

    python -X utf8 -B 05-knowledge/results/prime_power_plus_one_20261029/verify.py
    python -X utf8 -B 05-knowledge/results/prime_power_plus_one_20261029/seal.py

The verifier requires Python and SymPy, and imports the pinned
regular-parent verifier for native isomorphism/module routines. It runs
1,077,498 explicit local truth gates, retained under python -O, plus the
inherited hostile verifier's gates. Normal and optimized outputs agree.
The certificate SHA256 is
ddd23b5271c72ec5a314795345bb5aad11fe0ee8df8de1e18952d969075f67c9.

The complete candidate-module universes at q=7,11,15,19,23,27 give
35,77,135,209,299,405 minimum deletion sets. Independent induced-subset
module enumeration checks q=7,11. Other controls include16 complete
norm-clock charts,56 factorizations, exact lattice Smith forms through
q=59, complete symplectic-plane enumerations over F2,F3,F5,F7,
full native groups of the named local bases, and the inherited
81/27 pair with1620/1728 regular binary masks.

All-size theorems rest on the written proofs; finite computations do
not replace their quantifiers. LRC(14), general tournament
Erdős--Hajnal, universal Collatz convergence, Sha finiteness, and
Fermat-prime completeness remain OPEN. No new full A000568 term or
enumeration speedup is claimed.
