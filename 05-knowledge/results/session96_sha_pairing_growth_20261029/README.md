# Sha, doubling, and the pairing data that F1 shadows forget

2026-10-29. Source pin `d92d6c0d742c7bab872d0e62769e70e9f0595d2c`.
This session follows the requested connection between the Tate–Shafarevich
finiteness conjecture, square/twice-square orders, field-with-one-element
ideas, and tournament growth modes. It yields exact bridges, obstructions,
and a finite-prime stopping criterion. **It does not prove general Sha
finiteness or construct an F1 model of Sha.**

## 1. The arithmetic theorem is stronger than a numerical square test

[Poonen–Stoll](https://math.mit.edu/~poonen/papers/sha.pdf) show that a finite
Sha of a principally polarized abelian variety has underlying group
T×T or C2×T×T. The pairing is antisymmetric, not necessarily alternating.
Its characteristic element c satisfies b(x,x)=b(x,c),2c=0. The value b(c,c)
determines the square/twice-square alternative, while the entire element c
is needed to decide whether the given pairing is alternating. A square-order
example can have c nonzero but b(c,c)=0. The finite nondivisible quotient has
the same structure; this does not remove a possible divisible radical.

The principal-polarization hypothesis is essential. [Konstantinou,
Theorem1.1](https://arxiv.org/html/2404.16785v1) proves that every positive
squarefree integer occurs as the squarefree part of a finite Sha order for
some abelian variety over Q. This cited result concerns the larger class,
not the principally polarized setting of the requested square/twice-square law.

[Marked finite pairings](pairings/PAIRING_DEFECT_RECURRENCES.md) gives an
independent finite-algebra proof of the characteristic correction and
hyperbolic peeling. The arithmetic inputs and classical algebra are credited;
the session contribution is their precise connection to the other carriers.

Writing |G|=2^epsilon R², orthogonal products obey the exact carry law

    (epsilon,R)*(eta,S)=(epsilon XOR eta, 2^(epsilon eta)RS).

Two odd defects make a square cardinality, but the original doubled pairing
can remain nonalternating: the marked element becomes(c,c), not zero.
For exponent-two pairing objects, adjoining a hyperbolic plane changes
dimension by+2, whereas duplicating the object multiplies dimension by2 and
squares its cardinality. These are different operations.

## 2. An actual tournament-to-pairing map, with a sharp obstruction

[Tournament pairing growth](tournaments/PAIRING_GROWTH_MODES.md) recovers
the repository's precise growth constructions. For an even-order tournament
with integral skew dominance matrix M, use

    G_M=Z^n/M Z^n,
    b_M([x],[y])=x^T M^(-1)y modulo Z.

This is a well-defined perfect alternating pairing. Moreover M modulo2 is
J-I and is invertible for even n. Thus its determinant is an **odd square**.
This inherited determinant fact has a decisive consequence for the proposed
bridge: these tournament groups have no2-primary part at all. They cannot
by themselves represent Sha's exceptional characteristic defect.

The growth modes retain different amounts of information:

* Endpoint padding source ->T->sink, n->n+2, is integrally congruent to
  M direct-sum J2. It preserves the entire paired discriminant group.
* Self order join T->T, n->2n, squares the group order. It need not produce
  a direct square of the group. The explicit order-four witness has
  G_M=C3×C3, but its self join has C9×C9, not C3^4.
* Skew-Sylvester doubling has an explicit integral quotient map to
  coker(I-2M²), together with its transported pairing. The diagonal twin
  corrections and the named integer coordinates are retained.

The second example explains why doubling/squaring is not a cardinality-only
duality: the boundary attachment can change extension data. The note gives
a sufficient marked-port condition for direct-sum preservation and an exact
exponent bound. No arithmetic realization of these groups as Sha is asserted.

### A fresh adjacent cokernel supplies a useful negative control

The incoming [golden circulant cokernel tower](../stabilizer_necklace_kernels_20261028/COKERNELS.md),
published at05c5d4362e4, retains a cyclic clock and has depth-two order4864,
with Smith factors1, seven copies of2, and38. Its19-primary valuation is
one:4864=2^8*19. Therefore that particular finite group cannot carry any
perfect antisymmetric Q/Z pairing. Its order is neither a square nor twice
a square. This does not contradict the tournament construction above: the
clocked sign-word circulant is a different matrix, not an even tournament
skew dominance matrix. An exact additive cokernel and a preserved clock
still do not supply the required paired arithmetic structure. The incoming
source's general tower proof is credited; this packet uses only its explicit
depth-two group and the proved finite-pairing obstruction.

## 3. F1 gives a paired skeleton, but hides coupled sectors

[Symplectic skeleton and gluing](f1/SYMPLECTIC_SKELETON_AND_GLUE.md) uses
classical finite symplectic geometry and the Weyl-group interpretation of
q=1. The complete isotropic flag polynomial specializes to2^n n!, while the
Lagrangian count specializes to2^n. These are signed paired-coordinate
choices, not literal matrices over a one-element field.

For two symplectic summands of dimensions2m and2n, the exact Lagrangian
composition recurrence is

    L_(m+n)(q)=sum_r I_(m,m-r)(q) I_(n,n-r)(q) |Sp_(2r)(F_q)|.

The index r measures the rank of the hidden coupling between the summands.
After choosing their isotropic intersections, the remainder is the graph of
an anti-symplectic isomorphism. This supplies an actual geometric map behind
the recurrence.

At q=1 the r-sector vanishes to exactly order r, with leading coefficient

    2^(m+n-r) r! binom(m,r) binom(n,r).

Every mixed sector disappears at the constant term. Thus thin Lagrangian
counts multiply, even though finite-field Lagrangians do not simply split.
Over F2 in the smallest case, there are **15 Lagrangians:9 separable and6
mixed**. At q=1 the total and separable counts both become4.

This refines the earlier F1 root-algebra warning: multiplicative structure
and a simple q=1 number can survive while additive gluing data disappear.
Here the loss is measured by a hierarchy of vanishing orders. It does not
provide a q-variable for Sha, prove torification, or recover its characteristic
2 data without additional structure.

## 4. What would actually certify a piece of Sha finiteness?

[Selmer tower certificates](selmer/SELMER_TOWER_FINITE_CERTIFICATES.md)
separates finite order shape from stopping. For any p-primary torsion group
with finite p-torsion, let b_n=log_p|G[p^n]|. Its increments are nonnegative
and nonincreasing. A single plateau b_(n+1)=b_n forces the entire p-primary
group to equal G[p^n]. The proof does not require a conjectural cutoff.

The usual Selmer exact sequence turns this into a usable certified test.
Suppose r0 is a proved Mordell–Weil rank lower bound, the torsion ratio

    tau_n=|A(K)_tors[p^(n+1)]|/|A(K)_tors[p^n]|

is certified, and L_n<=|Sel_(p^n)|, U_(n+1)>=|Sel_(p^(n+1))| are certified
bounds. Then

    U_(n+1) <= L_n p^r0 tau_n

forces the actual rank to be r0 and Sha[p^infinity]=Sha[p^n] to be finite.
The proof is an exact sandwich of group orders. Both bound directions and
the torsion correction matter. This is a criterion for verified arithmetic
inputs, not a theorem that every descent search reaches it.

A strong hostile explains the remaining difficulty. For any fixed depth N,
a sufficiently large finite hyperbolic group (Z/p^E)^2 and the divisible
group(Q_p/Z_p)^2 have identical first N torsion layers, inclusion maps and
restricted pairings. The finite ambient pairing can be perfect while its
observed low layers are zero. No bounded growing prefix forces eventual
stopping. Global finiteness additionally requires control of prime support.

## 5. An actual arithmetic calibration, not only abstract models

[Exact two-descent control](example/EXACT_TWO_DESCENT_CONTROL.md) reproduces
a classical calculation for E:y²=x³-x. Full2descent reduces possible Selmer
classes to sixteen signed2-supported square-class triples. Real signs leave
eight; primitive congruence equations modulo8 remove four. The survivors
are precisely the four rational2-torsion Kummer images.

This proves |Sel_2|=4, rank(E(Q))=0 and **Sha(E/Q)[2^infinity]=0** using the
standard descent interface. The verifier checks all32,768 residue vectors
and supplies explicit rational points on the four surviving coverings.
The additive equations and projective primitivity are essential. Modular
survival is not promoted to solubility without those rational witnesses.

The example is known mathematics reproduced as a calibration of the
certificate, not a new elliptic-curve result. It gives no conclusion about
odd-primary Sha from this calculation alone.

## 6. The proposed unifying framework and its obligations

The appropriate common object is a **marked paired carrier with explicit
attachment and level maps**, rather than its order or a thin matching
skeleton. This phrase is a research organization, not a claim of a new
universal equivalence. The packet supplies three concrete comparisons:

1. Integral tournament maps preserve a genuine alternating discriminant
   pairing, but exclude all2-primary arithmetic defects.
2. Symplectic direct-sum gluing has exact coupled sectors invisible at q=1.
3. Arithmetic descent needs localization conditions, actual image bounds,
   Mordell–Weil contributions and the p-power stopping clock.

A future arithmetic construction must preserve these data before its
recursive simplicity can constrain Sha. The finite-pairing recursion alone
cannot supply missing height or prime-support bounds. The most actionable
next target is to implement the certified ratio test for a chosen nontrivial
higher descent, with independently checked local images; that implementation
is not claimed in this packet.

## Verification

Every producer script and the [independent audit](audit/AUDIT.md) uses exact
arithmetic with explicit checks active under normal and optimized Python.
The committed packet is replayed from a clean checkout, and file bytes and
relative links are checked after publication. Enumeration scopes are stated
in each note. Existing classical theorems, new repository syntheses, finite
controls, and open arithmetic obligations are kept separate throughout.
