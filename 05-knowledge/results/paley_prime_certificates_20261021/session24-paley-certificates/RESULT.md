# Explicit minimum signed-regular Paley certificates at prime order

Status: PROVED elementary corollary, with a CITED classical group theorem and FINITE-EXACT controls. No claim of literature novelty. Prime fields only.

## Inheritance and boundary

The closest mechanism is THM-4743, `01-canon/theorems/THM-4743-converse-extended-gluck-for-tournaments-primitive-classification-paley-tournaments-and-c3-towers.md`: the Paley exception classification already says precisely q=3,7, even for prime powers. This note does not claim that classification as new. Its refinement is an explicit fixed three-point certificate, its minimal size, and an exact count at prime order.

The canonical hostile is P7. The relevant repaired near miss is THM-309, `01-canon/theorems/THM-309-five-cycle-design-paley.md`, section "Proof Sketch": full affine maps are NOT all ordinary automorphisms; nonsquare slopes reverse orientation. We use that repaired distinction, not the unrelated cycle-count assertions in that file. The small-sidecar control is the full permutation enumeration at p=3,7 below. A targeted search of these sources and the current theorem files did not locate the explicit all-prime {0,1,3} certificate or the triple count. This is a bounded inheritance check, not an exhaustive novelty search.

Concept board: setwise stabilizers; affine normalization; centroid multiplier groups; ordinary versus signed symmetry; the characteristic-seven quadratic obstruction. No implication about Hamiltonian-path gaps 7 or 21 follows.

## 1. Affine triple classification

Let p>3 be prime and A=AGL(1,p). For a three-set S, its centroid c=(sum S)/3 is fixed by every setwise stabilizer element. Translation by -c therefore embeds A_S in F_p^*. Its action on S is faithful (an affine map fixing two distinct points is identity), so A_S is also a subgroup of S3. Being cyclic, it has order 1, 2 or 3.

Normalize S to {0,1,t}, t not 0,1. An involution induces a transposition, and its unique fixed point is the midpoint of the transposed pair. Thus order 2 occurs exactly for

    t in {-1, 2, 1/2}.

These three parameters are distinct for p>3. An order-three element permutes the points in a 3-cycle. Taking the cycle 0 -> 1 -> t -> 0 forces x -> (t-1)x+1 and hence

    t^2-t+1=0.

The inverse cycle yields the same condition. The quadratic has exactly two roots when p=1 mod 3, and none when p=2 mod 3: its roots are primitive sixth roots for p>3. These roots do not belong to the midpoint set (substitution gives 3 or 3/4). Thus these two cases are disjoint and exhaustive.

For t=3, midpoint failure means p divides 4, 1, or 5, hence only p=5 among p>3. The three-cycle failure means p divides 7, hence p=7. Therefore

    Stab_AGL(1,p)({0,1,3}) = 1 for every prime p>=11.

At p=3 this is not a three-set; that case is excluded rather than substituted into the formula.

## 2. Exact counting

Each unordered triple is represented exactly six times by an ordered choice of its first two distinct elements u,v and the normalized third parameter t=(w-u)/(v-u). For every permissible t there are p(p-1) choices of u,v. Write I=1 if p=1 mod 3, otherwise I=0. The numbers of triples with stabilizer orders 2,3,1, respectively, are therefore

    N2 = p(p-1)/2,
    N3 = I p(p-1)/3,
    N1 = p(p-1)(p-5-2I)/6.

The number of free affine orbits is (p-5-2I)/6. In particular it is 0 at p=5,7 and 1,2,4 at p=11,19,31. The formula counts labelled subsets; division by p(p-1) is justified only for the regular triples.

## 3. Transfer to the FULL signed Paley group

Classical input: Peter Müller, "Permutation groups of prime degree, a quick proof of Burnside's theorem", arXiv:math/0310200v1 (2003), Proposition 1, https://arxiv.org/html/math/0310200, journal DOI https://doi.org/10.1007/s00013-005-1421-z. The proposition states that a permutation preserving a nonempty proper subset of nonzero field differences is affine. Its following proof of Burnside's theorem also explicitly obtains the affine subgroup in the original translation coordinates. This avoids the potential gap of merely saying "abstractly isomorphic to an affine group".

Apply it with the quadratic residues Q in F_p, for p=3 mod 4. Every ordinary automorphism of P_p preserves Q-differences, so is affine; its multiplier must be a square. Conversely every square-affine map preserves P_p. Thus Aut(P_p) is exactly the square-affine group. Negation reverses orientation, since -1 is a nonsquare. Composing any anti-automorphism with negation gives an ordinary automorphism. Consequently

    Aut±(P_p) = AGL(1,p).

Equivalently one can apply Burnside's prime-degree theorem: Aut(P_p) contains translations and is not 2-transitive because ordered arcs and reversed arcs are distinct orbits. The ordinary group's odd order is true but unnecessary for this step.

For p>=11, {0,1,3} is therefore a signed-regular subset of P_p. No subset of size below three is signed-regular: the empty set is fixed by everything; a singleton {u} is fixed by x -> 2u-x; a pair {u,v} is fixed by x -> u+v-x. These are nonidentity anti-automorphisms. The minimum size is exactly three.

At p=3 or 7, no regular subset exists: by complementation it suffices to consider sizes at most floor(p/2). Sizes 0,1,2 fail as above; at p=7 every triple fails by N1=0. Thus the prime-order exception classification also follows directly here without Gluck or Feit–Thompson.

Do not transfer this fixed triple to prime powers. The Frobenius map fixes every prime-subfield point, and for a nonprime Paley field it is a nonidentity ordinary automorphism fixing {0,1,3} pointwise (when that is a three-set). The established prime-power existence theorem remains inherited from THM-4743, not proved by this certificate.

## Reproduction and scope

Run from the workspace:

    python3 session24-paley-certificates/verify.py
    python3 -O session24-paley-certificates/verify.py

The script uses explicit RuntimeError checks, not removable assertions. It exhaustively tests all triples against all affine maps for primes 5<=p<=31; independently tests every normalized parameter via its six possible ordered image pairs for primes 5<=p<=1009; checks {0,1,3} at each such prime; tests all subsets at p=3,7; and independently enumerates all permutations at p=3,7 to recover the full signed Paley group. Normal and optimized outputs must byte-match. These finite controls validate implementation and hostiles, not the universal quantifiers of the proofs.

Artifacts and source-file hashes are recorded in receipt.json. No repository file was modified by this task.
