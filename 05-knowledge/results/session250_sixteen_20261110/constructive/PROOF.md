# Rooted ordinal gluing and the strong-core normal form for unit-offset seeds

Status: PROVED IN THIS NOTE; no independent canon audit. Chan–Pak remains OPEN.
No assertion about what the historical `1-5-4-5-1` decomposition represents
is needed. This is a complementary constructive boundary-state mechanism.

**Main result:** every nonconstant unit-offset rooted seed reduces to a
strongly connected deletion core, with a transitive prefix and suffix.
Removing that padding preserves the odd part of its ray slope. The exact
seven-coordinate gluing law proves this reduction; see Theorem5.

## Inheritance and scope

THM-5126, `01-canon/theorems/THM-5126-unit-offset-rays-two-adic-lines-isolation-of-7-and-21-extension-mean-law.md`,
gives the inherited ray `H(L o 0^t)=c+2^t M`, where the root has out-set O
in its deletion B, `c=H(B)-T`, and T counts ordered two-path covers whose
first end and second start both belong to O. The closest mechanism is this
cut-defect law and THM-3121's path-cover substitution kernel.
Session249 `sequence_gf_20261110/grammar/PROOF.md` gives positive support
grammars; it does not justify forgetting the endpoint data below.
THM-5128, 5129 and 5130 are RESERVED empty stubs in the inspected HEAD,
and are not dependencies. The present results are proved directly.

Hostile: scalar H does not determine gluing at a root. Corrected near miss:
an empty inserted block is not the r=0 specialization of the nonempty-gap
formula. Least-used sidecar: the two opposite ordered-cover boundary sums.
Concept board: actual path support; unit-offset defect; ordinal gluing;
endpoint ownership; positive dyadic gap equations. No external literature
priority claims are made.

## 1. A seven-coordinate exact gluing law

Let B be a nonempty tournament and O a selected subset. All paths below are
directed and nonempty unless expressly allowed empty. Define

* d = H(B).
* q = number of ordered pairs of paths covering B, allowing either empty.
* S, E = numbers of Hamiltonian paths starting, respectively ending, in O.
* M = number of ordered covers (P,Q) with P possibly empty and Q starting in O.
* N = number of ordered covers (P,Q) with Q possibly empty and P ending in O.
* T = number of ordered covers (P,Q) with P ending in O and Q starting in O.

Thus q=2d+F_B(2)=Q_B(2). The letters S,E are endpoint sums on B,
not the similarly named endpoint count on the rooted extension in THM-5126.

**Theorem 1.** For nonempty coloured blocks A,C, their ordinal join
B=A▷C, with O=O_A union O_C, has signature

    d = d_A d_C,                    q = q_A q_C,
    S = S_A d_C,                    E = d_A E_C,
    M = M_A q_C + d_A M_C,
    N = N_A d_C + q_A N_C,
    T = T_A d_C + d_A T_C + E_A S_C + M_A N_C.       (1)

Proof. Every path in an ordinal join has its A segment before its C
segment. Restricting two covering paths to the factors is a bijection to
independent pairs of factor covers, allowing empty paths. This proves q;
the single-path formulas follow the same way. For M, if Q starts in A,
choose a cover counted by M_A and any two-slot cover of C. If Q starts
in C, all A is in P, giving d_A M_C. Reverse the endpoint reasoning for N.

For T, classify the ordered pair (end(P),start(Q)) by its factor owners.
Owners (A,A) force all C into Q, giving T_A d_C. Owners (C,C) force all
A into P, giving d_A T_C. Owners (A,C) force P to be all of A and Q all
of C, giving E_A S_C. Owners (C,A) allow both slots in both factors,
giving M_A N_C. These cases are disjoint and exhaustive. QED.

The law is associative because it computes a literal combinatorial
signature. It constructs the actual joined tournament and its rooted ray,
not merely a formal sequence of candidate values. All seven coordinates
are nonnegative and the rule is polynomial with nonnegative coefficients.
For a fixed right factor, appending a constant coordinate makes it a
nonnegative integer linear transfer, directly compatible with Session249.

**Corollary 2 (unit seeds do not freely concatenate).** Write c_i=d_i-T_i.
Then

    c_(A▷C)-1 = (c_A-1)d_C+(c_C-1)d_A
                -(d_A-1)(d_C-1)-E_A S_C-M_A N_C.     (2)

If both factors have unit offset, their join has unit offset iff
`(d_A-1)(d_C-1)+E_A S_C+M_A N_C=0`. In particular two unit-offset
factors with d_A,d_C>1 never give a unit-offset joined root. This is an
exact obstruction to a tempting multiplicative seed-generation shortcut.

## 2. A dyadic seam with an exact stopping rule

Insert between A and C an unselected transitive tournament P_r of r>=1
vertices: B_r=A▷P_r▷C, O=O_A union O_C. Its signature is
`(d,q,S,E,M,N,T)=(1,2^r,0,0,0,0,0)`.

**Theorem 3 (seam equation).** Put

    alpha=d_A d_C-T_A d_C-d_A T_C,     beta=M_A N_C.

Then the rooted extension with this deletion and out-set has

    c_r = alpha-beta 2^r,
    M_r = M_A q_C 2^r+d_A M_C.                       (3)

If both selected sets are nonempty, beta>0 and there is at most one
positive gap length producing a unit-offset seed. It exists exactly when
`(alpha-1)/beta` is an integral power of two at least two. When it exists,
construct B_r and its root explicitly; the inherited tail law certifies
every `1+2^t M_r`, t>=1, as an actual Hamiltonian-path value.

Proof. Apply (1) twice. Each y in O contributes at least H(B-y)>0
to M by taking Q=(y), and at least H(B-y)>0 to N by taking P=(y).
Thus beta is positive precisely when both O_A,O_C are nonempty.
Solving c_r=1 gives the stated test. QED.

If beta=0, the defect is constant and the criterion is alpha=1.
For the empty gap, instead

    c_empty=alpha-beta-E_A S_C.                      (4)

The E_A S_C contribution disappears as soon as a nonempty unselected gap
is inserted. Equation (4) is essential, not an optional convention.
For A=C=K1 with their vertices selected, alpha=beta=1: empty gap has
c=-1, whereas r=1 has c=-1 and r=2 has c=-3. Naively substituting r=0
in (3) would incorrectly give c=0 (also violating oddness).

## 3. Strong-core reduction: the seam obstruction is universal

**Lemma 4 (path-cut injection).** For every vertex y of a nonempty B,
`U_B(y)>=H(B)` and `V_B(y)>=H(B)`, where U and V are the individual
summands of M and N. Cut a Hamiltonian path immediately before y to get
the U-cover, or immediately after y to get the V-cover. Concatenation
recovers the path, so both maps are injective. Consequently
`M>=|O|d` and `N>=|O|d`.

**Theorem 5 (exact unit-offset core normal form).** If O_A and O_C are
both nonempty, then the ordinal join has

    c_(A▷C)<=d_A d_C-M_A N_C<=0.                    (7)

Consequently the nontrivial unit-offset rooted tournaments are EXACTLY
those whose deletion has the form

    B=TT_a▷C▷TT_b,       O nonempty and contained in C,           (8)

where C is one strongly connected component (a singleton is allowed),
the root attached to C with out-set O has unit offset, and a,b>=0.
The signature transports as

    c_B=1,                M_B=2^b M_C.             (9)

The remaining case O empty has unit offset exactly when B is transitive;
then M=0 and the ray is constant1.

Singleton boundary: when C=K1 and O is its vertex, its signature is
`(1,2,1,1,1,1,0)`, so (8) gives M=2^b and the line of odd slope1.
For a non-singleton example, take C the directed triangle and O any two
of its vertices. Its signature is `(3,12,2,2,6,6,2)`, c=1 and M=6;
its rays begin13,25,49. A prefix of length2 and suffix of length3 gives
`(3,384,0,0,48,24,2)`: still c=1 and odd slope3, with first ray point97.
Both signatures are recounted by the companion, without using (1).

Proof. Inequality (7) follows from (1), nonnegativity and Lemma4.
The SCCs of a tournament are totally ordered by their common direction:
two distinct SCCs cannot have arcs in both directions (they would merge),
and their acyclic quotient tournament is transitive. If O meets two SCCs,
cut between them and apply (7), ruling out c=1. If O lies in one SCC C,
all other factors are unselected. Formula (1) gives
`c_B=c_C product_(other SCCs D) H(D)`. For this integer to equal1, every
factor must equal1 and c_C=1. A tournament has H=1 iff it is transitive:
a directed triangle has three paths, which inject into paths of any
containing tournament by deterministic insertion of the other vertices;
nontransitive tournaments contain a directed triangle. A non-singleton
strong component is nontransitive, hence has H>=3. Thus all other SCCs
are singletons. Applying (1) to the transitive prefix and suffix proves
(9), and proves the converse as well. QED.

The insertion used above is elementary: for a path v_1,...,v_n and new
x, prepend x if x->v_1; otherwise insert it just before the first v_j
that x beats, or append it if there is no such v_j. Deleting x recovers
the old path, proving injectivity.

This is the constructive payoff: removing all external SCC padding
preserves the odd part of the ray slope. For a fixed core C,O, every
padding variant is explicitly reconstructed by (8); trailing padding
only starts the same odd line later. Thus arbitrary ordinal gluing cannot
create a new odd line from separate marked components. It suffices to
search strongly connected deletion cores for all nonconstant unit-offset
ray seeds. This reduces that mechanism, not the general Chan–Pak problem.

The apparent positive-gap search of Theorem3 actually has no solution
when both marked sides are nonempty: alpha<=d_A d_C<=beta, so
alpha-beta*2^r<=-d_A d_C for r>=1. The power-of-two criterion remains
exact, but Lemma4 explains structurally why its test always fails here.

## 4. Many gaps: a positive defect equation

Take nonempty coloured factors B_1,...,B_k separated by unselected
transitive gaps of positive lengths r_1,...,r_(k-1). Write D=product d_i.
Then the same induction gives

    T = sum_i T_i product_(h!=i) d_h
        + sum_(i<j) M_i N_j
          product_(h<i or h>j) d_h
          product_(i<h<j) q_h
          2^(r_i+...+r_(j-1)).                       (5)

Proof. A cover counted by T has first-end owner j and second-start owner
i. The case i=j gives the first sum. The case i<j gives the second sum:
outside the owner interval there is only one occupied path slot; inside
it both slots are allowed, contributing q_h and 2^r for gaps. The case
j<i is impossible: an intervening nonempty unselected gap could belong
to neither path without changing one marked endpoint. QED.

Thus unit offset is exactly the positive exponential Diophantine equation

    sum_(i<j) C_ij 2^(r_i+...+r_(j-1))
       = D-1-sum_i T_i product_(h!=i) d_h = K,        (6)

with explicit nonnegative integer C_ij from (5). Every positive coefficient
implies `r_i+...+r_(j-1)<=floor(log2(K/C_ij))`. If K<0 no solution exists;
if K=0 every positive coefficient rules out solutions. Any gap crossed by
a positive coefficient is therefore bounded in every solution. Uncrossed
gaps are genuinely invisible to the defect and can vary freely. This is
an exact finite constructive search after removing those free coordinates,
not an unbounded brute-force size search. For unit offset, Theorem5 is
stronger: multiple selected factors are impossible. For other prescribed
offsets, (5) remains a usable exact seam equation. Recover M by (1), add
the root, and use its certified ray. It does not prove that every required odd line
has a seed, nor that arbitrary rooted tournaments admit these decompositions.

## Reproduction and connection ledger

`python3 examples.py` in this directory tests (1) against literal path
enumeration, tests seam and empty-gap formulas, and emits example signatures.
Tests are proof-derived examples, not a census or independent canon audit.
Source: coloured tournaments with actual endpoint owners. Target: seven
integers and polynomial gluing. Preserved: H, both rooted-cover boundaries,
unit-offset condition and ray slope. Lost: isomorphism and internal paths.
Required sidecar: the explicit factor tournaments and selected sets, retained
by the constructor. The positive gap equation is a concrete link between
rooted seed construction and support grammars; no formal root state is
declared realizable without its tournament construction.
