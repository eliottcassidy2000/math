# Fresh census audit and the exact mixed-prime descent target

**Status:** PROVED lemmas and reduction; FINITE-EXACT independent audit of the
ten classes in `README.md`. 20 October 2026. No claim beyond degree ten is made
by the census, and no new canon theorem identifier is reserved here.

## 1. Findings that can now be promoted

The completeness argument was checked independently at every reduction:

1. An odd group has no unordered-pair reverser and hence admits an invariant
   tournament. A transitive such tournament is regular.
2. A prime cycle of length strictly between half the degree and the degree
   conflicts with regularity at a fixed vertex. This rules out both 5 and 7
   in an odd transitive degree-nine group, leaving a 3-group.
3. At degrees 5 and 7 the four and eight circulant tournaments give exactly
   the stated transitive images. Their vertex-stabilizer arguments do not
   require a classification of finite simple groups or prime-degree groups.
4. The remaining actions are subdirect products on their intrinsic vertex
   orbits. The only nontrivial equal-order coupling is the diagonal F21 on
   7+3; its two quotient maps are conjugate by a root transposition.
5. The inherited all-p-group theorem covers the omitted p-groups, including
   the trivial group. Consequently the complete table proves that degree ten
   is the least failure degree, with precisely one failing ambient-conjugacy
   class there: C15 on 5+3+1+1.

Fresh runs reproduce the original JSON byte for byte. The original orbit
enumeration has 306 checks and its direct-cube/Burnside replay has 621.
The newly written `independent_automorphism_audit.py` shares no implementation
with those programs. It independently recovers 69 tournament types and 388
graph types, counted with repetition across the ten symmetry predicates, and
tests 10,212 full-automorphism stabilizer-chain generators. Its 14,122 checks
pass, identically in normal and optimized Python. In particular, the graph
sign audit does not reuse the original twin-quotient shortcut.

## 2. The exact copy-orbit coordinate

Fix odd K <= G=S_n. For a full graph type X write A=Aut(X), let chi be its
edge-reversal character, and put

    C_K(A) = {J <= A : J is G-conjugate to K}.

Suppose this set is nonempty. Let a_X be its number of A-conjugacy orbits.
Then a_X is exactly the number of N_G(K)-orbits of K-fixed presentations of X.

**Proof.** A presentation g.X is K-fixed exactly when g^(-1)Kg <= A. Replacing
g by n g a with n in N_G(K), a in A changes this subgroup by A-conjugacy, and
conversely A-conjugacy of the subgroups gives such a double-coset relation.
Thus the two orbit sets are naturally in bijection. This counts embedded
subgroups, not labelled generators or abstract embeddings modulo Aut(K).

If X is globally sign-odd, put A^+=ker(chi) and let a_X^+ count the A^+-orbits
on the same copy set. Then its relatively-even contribution is exactly

    z_X = a_X^+ - a_X.                                      (1)

Indeed an A-orbit splits into two A^+-orbits exactly when its stabilizer
N_A(J) lies inside A^+. That is precisely the condition that the associated
K-fixed presentation is N_G(K)-even. An unsplit orbit contributes zero to
both sides of (1); a split orbit contributes one. No division by an abstract
automorphism count or by |K| is missing.

For a globally even graph every presentation is relatively even and its
contribution is a_X, rather than the right side of (1). This distinction is
essential: (1) is only for globally odd graphs.

Consequently the general defect identity becomes

    tau_T(K)-tau_E(K)
      = sum_(even X admitting K)(a_X-1)
        + sum_(odd X admitting K)(a_X^+-a_X)
        - sum_(tournament T admitting K)(a_T-1).             (2)

This is the existing relative Walsh identity expressed entirely in terms of
automorphism-group actions on copies of K. It explains both degree-ten rows:
C15 has uncompensated splitting, whereas diagonal F21 has two tournament
fusion losses exactly compensating two odd-graph splitting contributions.

## 3. Necessary and sufficient conditions for lossless descent

For fixed K the following are equivalent:

1. Each full tournament type meeting Fix(K) is one normalizer orbit; each
   graph type meeting Fix(K) is one normalizer orbit; and no globally odd
   graph has a relatively-even presentation.
2. A_T acts transitively on C_K(A_T) for every such tournament; and A_X^+
   acts transitively on C_K(A_X) for every such graph, with A_X^+=A_X when X
   is globally even.

The tournament and globally-even graph equivalences are the double-coset
description. For an odd graph, A^+-transitivity implies A-transitivity and
no splitting. Conversely, A-transitivity plus no splitting gives exactly one
A^+-orbit. This proves the equivalence.

Odd ambient pronormality of K is a sufficient group-theoretic condition:
any two copies lie inside the relevant A, or inside A^+ on the graph side,
and the pronormal conjugator belongs to their generated subgroup. The
criterion above isolates precisely the weaker object-relative facts actually
needed. It is not asserted here that it supplies a strictly larger known
class of K, nor that it is necessary for equality alone. Diagonal F21 is an
explicit counterexample to that latter necessity.

Under these lossless conditions, the inherited q-local relative Walsh theorem
descends for every odd prime q: one bijection on full types admitting K
preserves simultaneously every extension predicate

    Aut(x) contains a G-conjugate of H,
    K <= H <= N_G(K),  H/K a q-group.

To check the predicate, if gHg^(-1) lies in Aut(x), both K and gKg^(-1) are
copies in Aut(x). Copy-transitivity supplies an automorphism a aligning them,
so ag normalizes K. Thus ambient containment equals normalizer-quotient
containment. Apply the inherited relative theorem. This proof keeps the
quotient prime, distinguished base, and embedded action; it does not claim
joint control at distinct quotient primes.

## 4. The precise next research target

The strongest currently established mixed-prime mechanism is a pronormal
odd base followed by one odd-prime quotient localization, including its
normalizer-tower variants. The census shows both how it succeeds and why
testing H itself for pronormality is too restrictive.

The next target is therefore not another unrestricted odd-group conjecture.
It is to characterize the copy-transitivity conditions in section 3 for
nonpronormal bases, or to replace lossless descent by an explicit cancellation
of the three terms of (2). Any proposed criterion must distinguish:

- C15 on 5+3+1+1: nonzero defect one;
- diagonal F21 on 7+3: nonpronormal H and sign escape, yet exact equality;
- C3 on 3+1+1+1: already equal by the p-group theorem despite fusion and
  sign escape, so it is a hostile to necessity of lossless descent.

A concrete computational extension is a complete degree-eleven subgroup
census with independent prime-degree classification, rather than appending
selected examples to this table. It is not performed or counted here.

Connection contract: K-fixed presentations -> embedded copies in Aut(x);
preserved datum is ambient K-presence; lost datum is the copy orbit and, on
odd graphs, its A^+-splitting. Equations (1) and (2) restore exactly those
sidecars. The branch-integrality and twin-packing results share the lesson
that weighted quotient information is not an unweighted pairing, but their
carriers and weights are different. No numerical or proof transfer is inferred
without an additional map.
