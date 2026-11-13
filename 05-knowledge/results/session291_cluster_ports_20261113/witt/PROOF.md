# Grading before collapse: a Witt-inspired tournament bridge

Status: elementary statements proved below; independently reviewed (see ../audit/REVIEW.md). This is a
restricted source-inspired bridge, not an application of a de Rham–Witt theorem
to tournaments and not a proof of Chan–Pak. Classical necklace and determinant
identities are credited as such.

## 1. The marked necklace carry

Let P be a finite set of distinguished paths, each with a monomial weight w(P)
in any chosen nonnegative endpoint/port variables. Distinguishable paths may
have equal weights. Put F=sum_P w(P), and let psi_p substitute x_i^p for x_i.
For prime p define

    C_p(F) = (F^p - psi_p(F))/p.

This has nonnegative integral coefficients. Indeed cyclic rotation acts on
p-tuples of paths. Its fixed tuples are precisely the constant tuples, whose
weight sum is psi_p(F). Every other orbit has size p and constant monomial
weight. Thus C_p counts those free orbits, with their marks intact.

Direct expansion gives the exact carry and product rules

    C_p(F+G) = C_p(F)+C_p(G)
                + sum_(1<=j<p) binom(p,j)/p F^j G^(p-j),
    C_p(FG) = psi_p(F) C_p(G)+psi_p(G) C_p(F)+p C_p(F)C_p(G).

This C_p is NOT a named Witt norm N, and has the opposite sign from the
conventional delta associated to the Frobenius lift psi_p. The useful feature
is the required mixed correction, not an identification of these operations.

More generally the primitive cyclic words of length p^r have weight inventory

    L_(p^r)(F) = (F^(p^r) - psi_p(F)^(p^(r-1)))/p^r.

The subtracted words have period dividing p^(r-1); all remaining words have
full period. This proves integrality and positivity without a congruence-only
argument. It is the classical primitive-necklace mechanism with path marks.

### Physical-realization obstruction

If F is the inventory of all Hamiltonian paths of a tournament and h=F(1),
then h is odd. For every odd prime p and r>=1,

    L_(p^r)(F)(1) = (h^(p^r)-h^(p^(r-1)))/p^r

is even: its numerator is even and its divisor odd. It therefore cannot itself
be the Hamiltonian-path count of a tournament. For h=1 it is zero, also
impossible for a nonempty tournament. Even an exact positive combinatorial
orbit construction is not automatically a tournament construction.

The ordered sum of p copies of T physically realizes h^p. But rotation of its
path tuples is not in general an automorphism of that tournament: the ordered
components distinguish their positions. No deletion, quotient or division by p
with the required count follows. This is precisely the additional realization
obligation any norm-inspired Chan–Pak argument must meet.

## 2. Recovering what ordinary spectral ghosts forget

There is nevertheless a useful exact carrier enlargement. Work in

    R = Q[z_v : v in V]/(z_v^2 : v in V),  B=diag(z_v) A(T).

The variables remember vertex resources. Any repeated-vertex walk has zero
weight. All identities below are finite because the positive-degree ideal is
nilpotent. Define

    Q_T(z) = det(I+B)/det(I-B).

The formal log-determinant identity gives

    log Q_T = 2 sum_(k odd) tr(B^k)/k
            = 2 sum_(simple directed odd cycles C) z_(V(C)).

The second equality follows because repeated vertices vanish, and each simple
cycle of length k contributes k based walks. Exponentiating shows

    Q_T = sum_(vertex-disjoint odd-cycle collections K) 2^|K| z_(V(K)).

Overlapping cycles multiply to zero. Distinct cycles with the same vertex set
are counted distinctly. These are classical determinant/cycle-gas mechanisms;
the comparison with the cited Witt paper is about retaining the carrier.

By the repository's inherited odd-cycle formula (OCF), for every S subset V,

    H(T[S]) = sum_(U subset S) [z_U] Q_T.

Set H(empty)=1. In particular H(T) is the sum of all squarefree coefficients.
This summation is a LINEAR readout, NOT the ring substitution z_v=1: such a
substitution cannot be a homomorphism because z_v^2=0.

Ordinary tr(A^k) ghosts permit repeated vertices; the square-zero lift kills
those repetitions before the scalar readout. This explains structurally the
repository's overlap corrections in THM-505. It does not make the scalar
spectrum alone sufficient: the carrier still has 2^|V| basis elements.

## 3. A productive research obligation, not a claimed coverage lemma

The two constructions expose complementary requirements. Cyclic path orbits
supply integral arithmetic but fail the tournament parity/realization test.
Square-zero grading supplies exactly the disjoint-resource information needed
for OCF, but does not yet retain named endpoint feasibility. A constructive
next step must keep both marked endpoint data and disjoint-resource data until
composition is finished, and prove a physical tournament realization for the
resulting scalar. A useful quotient must be a congruence for those composition
operations, not merely equality of H or ghost values.

This suggests a testable search program: compare candidate port profiles under
all allowed attachments; refine profiles on the first differing continuation;
seek a finite or recursively described stable quotient. Neither finite success
nor analogy with Witt coordinates supplies an all-order stable quotient.

## 4. Endpoint-preserving carrier and exact seam elimination

The same resource ring has a more directly constructive companion. Let

    R_T = (I-AZ)^(-1), Z=diag(z_v),
    E_T(s,t,z) = sum_(i,j) s_i z_i (R_T)_(i,j) t_j.

The coefficient of s_i t_j z_S counts Hamiltonian paths of T[S] from i to j.
Proof: the finite geometric expansion of R_T lists directed walks; the initial
z_i and each arrival variable record every visited vertex exactly once, or kill
the walk if a vertex repeats. This also handles singleton paths. Named ports
must remain in s,t until the required feasibility conditions are evaluated.

For a split V=U disjoint-union W, write A in blocks A_U,C,D,A_W and put
R_W=(I-A_W Z_W)^(-1). Exact block inversion yields

    (R_T)_(U,U) =
      (I-A_U Z_U-C Z_W R_W D Z_U)^(-1).

All inverses exist: their constant term is I and every other term has positive
nilpotent degree. The effective term C Z_W R_W D Z_U records excursions from
U through W and back. Its iterates allow arbitrarily many excursions subject
to disjoint vertex resources. Thus one cannot simply delete this term when
contracting a block on the grounds that its total H or its two port totals
match another block. A substitution is guaranteed safe if the needed effective
operator and boundary entrance/exit readouts agree in the retained carrier.

This is ordinary resolvent/Schur-complement algebra, not a newly claimed
matrix identity. It provides a concrete sufficient interface equality test
and an exact recursively composable carrier for any finite tournament. Its
cost remains exponential in the original vertex-resource variables, including
eliminated interior vertices. Finding a smaller family-stable
representation is the genuine unresolved compression problem; finite matrix
size at a fixed order is not a uniform finite-state theorem.
