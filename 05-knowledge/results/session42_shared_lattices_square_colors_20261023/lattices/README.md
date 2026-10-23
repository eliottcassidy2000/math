# Session42: cyclic attachment lattices, larger kernels, and shared alignment

## Status and inherited boundary

PROVED for explicitly supplied actions and exact comparison blocks. This packet
extends the fixed-prime alignment quotient in the incoming
`odd_cell_completion_20261023/NETWORKS.md` at commit
`3f1a9a2c09d230bd77b141aabefd9458b031f707`. Its nine-cell kernel is retained,
not discarded. The all-size theorem below classifies a supplied chart problem;
it does not manufacture a compact complete observer for arbitrary native ports.

## 1. Every cyclic block transporter has two congruences

Let row charts be powers of a specified permutation r of order m and column
charts powers of c of order n. They must be genuine allowed internal charts
(for example automorphisms of the source cell). They need not be transitive.
Let (a,b) act on a binary block B by

    ((a,b)B)(i,j) = B(r^(-a)i,c^(-b)j).

For fixed B and B', define R={ (a,b) in C_m x C_n : (a,b)B=B' }.
It is empty unless the blocks belong to the same declared orbit. If nonempty,
choose t=(a0,b0) in R. Then R=t+H, with H=Stab(B). This follows immediately
by composing a transporter with the inverse of another transporter.

Lift H to the full-rank lattice L in Z^2. It contains (m,0) and (0,n).
There is a unique basis

    (h,0), (s,v), with h,v>0 and 0<=s<h.

Proof: L meets the horizontal axis in hZ; its second-coordinate projection
is vZ. Choose a vector of second coordinate v and reduce its first coordinate
modulo h to obtain s. Subtracting its multiples from every vector of L leaves
a horizontal multiple of h. Uniqueness follows from these intrinsic definitions.

Thus an attachment is equivalent to the following TWO integer congruences:

    b-b0 = 0 mod v,
    v(a-a0)-s(b-b0) = 0 mod hv.                         (1)

Equivalently, after the first congruence, a-a0=s((b-b0)/v) mod h.
The divisions are exact; replacing this by a naive modular division is unsafe.
The period conditions are h|m, v|n and h | s*n/v.
The relation has exactly mn/(hv) members, not necessarily one choice per row.

Horizontal one-sided kernel size is m/h. The minimum positive vertical
period is vh/gcd(h,s), so the vertical kernel has size n*gcd(h,s)/(hv).
These formulas expose rather than erase nonfunctional attachment choices.
Uniform blocks give h=v=1,s=0; rigid cyclic blocks give h=m,v=n,s=0.

## 2. A genuine nine-cell kernel fits the same lattice

Use the actual regular tournament R9=C3[C3,C3,C3], label vertices by
x=a+3u with a,u in C3, and let r:x->x+1 mod9. This is an automorphism:
it rotates the outer triple and applies an inner rotation when carrying.
The actual equitable attachment B(x,y)=1 iff x=y mod3 has row/column sum3.
For the declared C9 charts, H={ (a,b): a=b mod3 }, so

    h=3, v=1, s=1; |H|=27; both one-sided kernels have size3.

The full internal automorphism group of R9 has order81, not9. Its full
attachment relation has 2187 chart pairs and one-sided kernel27, as proved
and tested in the inherited packet. The cyclic calculation does NOT claim
those additional automorphisms disappear.

This gives a meaningful pair of descriptions: a compact cyclic subproblem
with two congruences, and the larger full-group problem with hidden modules.

## 3. Why an arbitrary chart list is not a cyclic coordinate

Even an abelian internal chart group can fail to be cyclic. The translations
of the two C3 coordinates of R9 form C3^2. Enumerate them as integer labels
code(u,v)=u+3v in {0,...,8}. Take the identity 9x9 cross block, and a target
obtained by applying the row translation (1,0). Its comparison relation is

    u_row=u_col+1 mod3, v_row=v_col.

In the artificial integer labels it is NOT a coset subgroup of C9 x C9.
Indeed the graph function at labels 0,1,2 is respectively 1,2,0; any affine
cyclic group graph with f(0)=1,f(1)=2 would have f(2)=3, contradiction.
The relation is perfectly affine in the correct four F3 coordinates.

Therefore 'every nonempty block relation has two congruences' requires two
CYCLIC actions, not merely two finite sets of charts. For abelian chart groups
of ranks d_i,d_j, the correct lift has d_i+d_j coordinates. Nonabelian chart
groups require the subgroup-coset formulation in Section6.

## 4. Global simultaneous congruences, never independent edge witnesses

Let a named cell i have declared cyclic group C_(m_i). Every attachment e=ij
supplies (1) in the SAME variables x_i,x_j. Unary ownership marks are retained
as their exact allowed chart cosets, or as explicit filters if not cosets.
Collect congruences as A x = b mod d, one modulus per row. Since they arise
from genuine actions, changing x_i by m_i leaves every row unchanged modulo d.

For the pure coset system, define

    Lambda={ z in Z^V : Az=0 mod d },
    M=diag(m_i) Z^V.

If x0 is a simultaneous solution, every solution is x0+Lambda modulo M.
The count is exactly

    product_i m_i / [Z^V:Lambda] = product_i m_i / det(Lambda).

An integer HNF/SNF computation supplies the lattice and its index. Existence
is the integer solvability of [A | -diag(d)](x,t)=b. Its Smith divisibility
conditions are exact. Any integer x can be reduced into its chart ranges,
with the auxiliary t adjusted because M is contained in Lambda.

This is an all-size exact decision/counting rule for supplied cyclic chart
relations, including composite orders and nontrivial one-sided kernels.
It is not a polynomial-size representation theorem for all native wiring.

For a well-defined additive chart observation Lx mod q (q divides every
L_i*m_i), constancy is equivalent to L annihilating a basis of Lambda mod q.
An arbitrary real wiring statistic or nonlinear electrical invariant is NOT
such a character and must be evaluated on the retained compatible states.

## 5. Exact composite-order joint alignment quotient

Fix one source orbit for each local block. Put G=product_i C_(m_i),
W_e=(C_(m_i) x C_(m_j))/H_e, W=product_e W_e.
These local block choices are torsors for W_e. The shared vertex charts act
by the homomorphism D:G->W given by endpoint restriction. Therefore

    W / image(D)

is the exact finite abelian group of joint attachment classes above the
chosen local types. Smith normal form, not a field dimension, is the natural
all-composite-order answer. In particular

    number of classes = (product_e m_i*m_j/|H_e|) * |ker D| / product_i m_i.

The inherited beta+r-1+a formula is its equal-prime, zero-or-slope-line
specialization. Composite groups can carry several distinct invariant factors.

For any connected network with at least two cells of the nine-cell cyclic attachments in Section2,
all x_i agree modulo3, while the lifts by3 remain independent. Hence
|ker D|=3^(|V|+1), W=C3^E, and the joint class group is C3^beta,
beta=|E|-|V|+1. Local kernels change the automorphism fibre dramatically but
do not create extra cycle coordinates in this particular family.

## 6. Full finite groups: quotient-functional does not mean functional

For arbitrary finite internal groups G_i, each block transporter is empty or
a coset of H_e <= G_i x G_j. Let P_i,P_j be the projections of H_e and let
K_i={g:(g,1) in H_e}, K_j={h:(1,h) in H_e}. Then K_i is normal in P_i,
K_j normal in P_j, and

    P_i/K_i is isomorphic to P_j/K_j.

Proof: send gK_i to hK_j whenever (g,h) belongs to H_e. Two choices of h
differ by K_j; two choices of g by K_i; multiplication is preserved, and both
projections are onto. Thus every relation is functional AFTER the exact
kernel quotient, but the kernel fibres must still be retained for shared ports.
This is the elementary subgroup argument often called Goursat's lemma.

After fixing reference internal isomorphisms between source and target cells,
the global transporter inside product_i G_i is empty or x0 K, where K is the
intersection of the pulled-back SOURCE stabilizers for all blocks and marks.
The proof is simply that the ratio of two global solutions fixes every
retained source object. This gives an exact all-size semantic specification,
not an assertion that K can always be cheaply computed or compressed.

A fixed literal indexed Hamiltonian path may be retained as an additional
marked object. Forgetting its positions, allowing multiple unspecified paths,
or imposing unrelated numerical filters can change the carrier or introduce
a UNION of cosets. Such unions must not be silently replaced by one lattice.
Converse is likewise a separately retained branch unless an actual specified
action encodes it. The native electrical mod8 condition remains an exact
filter, not a cyclic relation merely because it contains a modulus.

## 7. An all-size visible-kernel lower bound for ownership

For a fixed actual named-cell source, let K_i^vis be the intersection of the
one-sided stabilizers of ALL actual blocks incident to cell i, together with
its retained unary decorations. The direct product of the K_i^vis embeds in
the full cell-preserving automorphism group: each factor acts only in its own
cell and fixes every incident block. Choices in different cells are independent.

Let b(K_i^vis) be the minimum number of actual vertex owners whose pointwise
stabilizer in that permutation group is trivial. Any pointwise-recovering set
of owner marks must contain at least b(K_i^vis) marks IN cell i. Therefore

    minimum total owners >= sum_i b(K_i^vis).                 (2)

Proof: if the local marks fail to form a base, a nonidentity element of that
local kernel moves the source while fixing every mark everywhere. This is
an obstruction regardless of how connected the coarse support is. The bound
need not be sufficient: synchronized quotient symmetries can survive after
all independent local kernels have been killed.

For the full R9 equality-of-outer-coordinate family, at every cell incident
to a nonuniform attachment, K_i^vis=C3^3 and its
base size is3: one mark is required in each inner triple. Thus (2) recovers
the inherited3m lower bound. For the declared C9 subproblem, K_i^vis=<3>
has base size1 at a nonisolated cell; this difference records the carrier,
not a contradiction. An isolated cell instead retains its full internal
group (Aut(R9) or the declared C9); the same respective base sizes3 and1
still hold, but the asserted kernel subgroups are different.

Every new attachment shrinks a visible kernel by exact subgroup intersection.
A proposed cheap owner signature should therefore be tested against this
intersection before an expensive complete-chart search. Adding an attachment
with individually nontrivial kernel may still eliminate the old kernel if
the two kernels intersect trivially. Keeping each kernel size without its
embedded subgroup loses this information, just as marginal block orbits do.

## 8. Coupled integer moves after chart authentication

One may lift (1) to integer equalities with bounded chart variables and bounded
auxiliary quotient variables. Combine them with the actual bounded wiring
model only when their variables name the same owners and ports. The feasible
set is then an integer fibre C z=d, l<=z<=u. Session41's conformally primitive
kernel-move criterion applies to this FULL matrix: a projected linear output
is unique exactly when every primitive move feasible at one base point has
zero projection. Separate edge cycles are not enough after coupling.

The conceptual pipeline is therefore exact and ordered:
1. declare genuine groups, reference charts, marks and path clocks;
2. compute local transporter cosets without deleting kernels;
3. solve their shared-variable global compatibility;
4. retain the exact joint alignment class;
5. apply coupled move tests inside authenticated wiring states, preserving
   converse and actual electrical guards, and compare baselines across branches.

No step proves that a proposed compressed local inventory exhausts all native
tournament ports. That completeness remains the all-size research obligation.

## Reproduction

`python3 check.py` and `python3 -O check.py` use explicit failures rather than
assert statements. They exhaust small binary block actions, test HNF membership
and kernel counts, inspect actual R9 automorphisms and the noncyclic hostile,
and compare shared congruences against direct finite chart enumeration.
