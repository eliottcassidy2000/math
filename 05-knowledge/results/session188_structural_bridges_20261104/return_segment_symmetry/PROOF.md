# Factorial and power-of-two symmetry of tournament run coefficients

## 1. Setting and distinction from a separator theorem

Fix a tournament Q on named child labels, two distinct physical ports x,y with x->y, and arbitrary block-uniform incidences to the ports. Let D be the resulting augmented quotient. For a nonnegative child multiplicity vector k, write K_D(k) for the number of directed words that use child i exactly k_i times and each port once. Every adjacent pair must be an arc of D. Thus equal consecutive child symbols are forbidden. Child multiplicity zero is allowed in this coefficient notation.

These are the exact run coefficients in two-port Hamilton responses. This packet gives divisibility and parity constraints on the coefficients, not a theorem that an actual quotient-converse response is nonzero. The two-return swap used for modulo-eight profile reduction is already in `../coarse_mod8_affine/PROOF.md`; the new statements concern the full return-segment permutation action and its power-of-two parity contraction.

## 2. Return-segment action

Fix a label i with k_i=r>=2. Every word has a unique decomposition

    P i L1 i L2 i ... i L_(r-1) i S,

where P,S contain no i and each Lj is a nonempty word containing no i. Each Lj begins at an out-neighbor of i and ends at an in-neighbor of i. Permuting the r-1 return words therefore preserves every arc, all multiplicities, and the fixed prefix and suffix. It defines an action of the symmetric group on the full word set.

Every return word uses at least two distinct non-i symbols: a single distinct symbol would force either adjacent equal symbols or a directed two-cycle through i. If all labels except i and possibly one other child j occur at most once, no two return words can be equal. Equality would require at least two distinct other symbols to occur twice. Ports occur once and cannot supply either repeated symbol.

**Factorial divisibility.** If at most one child other than i has multiplicity greater than one, then

    (r-1)! divides K_D(k).

Indeed, every return permutation orbit has exactly (r-1)! elements. This holds for every r, every multiplicity of the one permitted repeated child, and arbitrary actual port incidences. In particular it strengthens the three-return parity controls; it is not restricted to r=3 or r=5.

The divisor is sharp. Put R=r-1 and take Q=C3[{i},{j},TT_R], with directions i->j->TT_R->i. Set k_i=R+1, k_j=R, and give each of the R remaining labels multiplicity one. Take x as a source and y as a sink. Every word is

    x i j v_(pi(1)) i j v_(pi(2)) i ... j v_(pi(R)) i y,

for a unique permutation pi. There are exactly R! words. This also covers R=1.

## 3. Exact parity contraction at a power-of-two return count

Set R=2^t with t>=1 and k_i=R+1; other child multiplicities are unrestricted. In a return permutation orbit, let m1,...,ms be the multiplicities of its distinct return words. Its size is

    R!/(m1! ... ms!).

This multinomial is odd precisely when adding the m_j in binary entails no carries. Since R has only one nonzero binary digit, that occurs precisely when s=1 and m1=R. Equivalently, use v2(n!)=n-s2(n): the orbit is odd iff sum_j s2(m_j)=s2(R)=1. Thus all odd orbits are singleton orbits whose return segments are identical. Even orbits cancel modulo two.

Let W run over all nonempty directed return words from i to i: W contains only child labels other than i, its successive arcs are directed, i points to its first symbol, and its last symbol points to i. Require R times the multiplicity of each label in W to be at most that label's k-coordinate. W contains no ports, because each port occurs only once whereas R>=2. Put

    k'_i=1,     k'_j=k_j-R times multiplicity_j(W) for j!=i.

Then the exact parity identity is

    K_D(k) = sum_W K_D(k')    (mod 2).                 (1)

The sum is finite. To prove the identity, every surviving word has the unique form P i (W i)^R S. Contracting its repeated return to P i S gives an arbitrary directed word counted by K_D(k'). Conversely, expand the unique i in any such word using W. Both directions preserve all exterior joins and both port names. Distinct W contribute distinct original words.

## 4. A quotient-triangle obstruction

Let H={j!=i : k_j>=R}. If K_D(k) is odd, at least one W in (1) exists and every label of W belongs to H. Its first label is an out-neighbor of i and its last label an in-neighbor. Along W take the first transition from an out-neighbor u to an in-neighbor v. Then

    i->u->v->i,

so Q restricted to {i} union H contains a directed triangle through i. In particular at least two other child labels must have multiplicity at least R. Failure of this triangle condition proves K_D(k) even for every port query.

This is a necessary condition, not sufficient for oddness. Even with such a triangle, port incidences or cancellation in the sum can give an even coefficient. Nor does a zero coefficient imply a collision of complete named responses.

The heavy-label threshold is sharp: for Q=C3 on i,j,k, multiplicities (R+1,R,R), and source/sink ports, there is exactly one word x i (j k i)^R y. Hence K=1. This also shows why unrestricted factorial divisibility would be false once two other labels repeat.

## 5. Reproduction and scope

The standard-library verifier counts run words by a remaining-multiplicity recursion, independently partitions explicit small word sets into return-permutation orbits, verifies the factorial controls, checks every summand in (1), and retains sharp and non-sufficient controls. Normal and optimized runs must agree. Finite controls support the proofs; they are not used in place of the all-order orbit argument.

The result concerns Hamilton run words, not ordered global path covers: cuts can create empty return pieces, so the same action is not asserted for cut words. The original child profiles and their odd scalar factors remain necessary when using these coefficients in a physical response. No physical attachment names are identified or changed.
