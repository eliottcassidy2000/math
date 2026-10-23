# Session37 independent statistics and parity audit

Status: PASS after the two precision corrections below. All new audit programs use only the standard library or C++17 and import no Session37 worker code. No Git operations or repository edits were made.

## Accepted conclusions

1. The unsigned path carrier is the commuting square-zero algebra R_n. Its paired Grassmann realization x_i -> theta_i eta_i is faithful, lies in the even sector, and preserves unsigned vertex relabelling. One odd generator per vertex instead gives an alternating path sum.
2. The squarefree determinant/permanent identity, mixed loop weight b+(-1)^(length-1)f, even-cycle response hierarchy, and determinant/permanent convolution are correct. The nilpotent quotient is essential; they must not be treated as ordinary univariate determinant/permanent identities.
3. The optimal exterior-generator counts are 2n for an even-sector embedding and 2n-1 for an unrestricted algebra embedding, over Q and n>=1. The lower-bound arguments are valid even for nonhomogeneous embeddings.
4. The nine-vertex determinant/permanent collision is real: both complete polynomials agree but H is 2843 versus 2841. The corrected no-collision statement through order eight is supported by the primary complete fixed-Hamiltonian-path census and the separate inheritance auditor's independent cycle census. This audit inspected those receipts and algorithm/coverage arguments but did not rerun their entire 2^21 universe.
5. The full exterior algebra Lambda(Lambda^2 Q^n) has the same S_n character as labelled tournaments. Its graph stabilizer sign is edge-reversal parity TIMES edge-permutation parity. It is not the Royle/Walsh graph predicate.
6. The surviving exterior graph algebra is genuinely supercommutative: internally admissible connected types with even edge number can repeat, while an odd-edge connected type can occur only once. The bigraded Euler product, connected balance, dyadic recurrence, odd-order complement symmetry, and homology determinant formula are correct.
7. Vertex-parity free-Lie-super PBW and the signed curl quotient are a separate valid construction. The extra even generators formed from odd primitive squares are indispensable. The signed cyclic quotient kills precisely even repetitions of odd-weight primitive roots; H does not descend to it as a supertrace.
8. The one-reversed-arc family U_m has the claimed exact odd-cycle polynomial, Hamiltonian specialization, recurrences, and non-real-rootedness for every m>=3. The previously known order-nine threshold is correctly credited rather than claimed as new.

## Precision corrections made before acceptance

### Linear species versus extra operations

Equality of the S_n characters for every n permits a choice of linear-species isomorphism under bijections. It is too strong to deny such a species isomorphism altogether. What is not supplied is a canonical intertwiner or compatibility with join, deletion, converse, or other operations. The source note now makes this distinction.

More strongly, tournament converse cannot be identified with exterior parity under an S_3 x C_2 intertwiner. On S_3-invariants, tournament converse has trace 2, whereas exterior parity has trace 0. The actual invariant degree polynomials at n=3 are 1+q^2 in the Walsh/tournament model and 1+q^3 in the exterior model.

### Empty-order forest convention

The linear-forest identity with prefactor 2^(n-1) needs n>=1. At n=0 it conflicts with the empty Hamiltonian convention. This qualification has been added; all other empty-object conventions audited here are consistent.

## Independent exact checks

### Characters and graph signs

`check_characters_and_cycles.py` evaluates signed pair-permutation orbit factors on all 272 conjugacy types for n=0,...,12. It checks total exterior trace equals the fixed-tournament count, all invariant degree multiplicities are nonnegative integers, and the exterior invariant superdimension vanishes for n>=3.

`check_graph_signs.cpp` separately enumerates all 1,253 graph isomorphism types through n=7 directly from labelled masks and relabellings. It evaluates reversal and selected-edge permutation signs independently. No graph atlas or graph package is used.

`check_components.py` compares the resulting degree polynomials with the character averages, expands the bigraded component product, and checks the connected balance and doubling correction. Connected even/odd counts for n=1,...,7 are

- b_n: 1, 0, 0, 1, 4, 22, 199.
- f_n: 0, 0, 1, 1, 4, 22, 199.

The difference from Royle-even connected counts already matters at n=6: b_6+f_6=44, but e_6=43 because the repeated triangle is excluded in the exterior model. The homology identity

det(g|C_1(G)) = det(g|H_1(G)) sign(g|V(G)) sign(g|components(G))

follows from the graph chain exact sequence. It confirms the infinite admissible families: P_n iff n=0 or1 modulo4 (n>=1), C_n iff n=3 modulo4 (n>=3), and K_n iff n is odd.

### Square-zero algebra

`check_squarezero.py` covers every labelled tournament through n=4, including empty order: 76 cases. It checks every subset coefficient of the squarefree MacMahon identity and the Hamiltonian determinant quotient. Seven (f,b) choices give 532 mixed-statistics checks. A separate signed weighted matrix with nonzero diagonal tests the MacMahon identity outside tournaments. All 1,365 paired exterior basis products through n=5 are checked.

The four-vertex hostile with outmasks [8,1,3,6] has five Hamiltonian paths but alternating path sum -1. Thus even taking the absolute value does not repair naive one-generator fermionization.

For the unrestricted minimality proof, the degree-two component of a commutator forces all degree-one parts of the commuting images into one line. Every product term containing two such parts vanishes, even with intervening higher-degree factors. A nonzero product of n images must therefore have exterior degree at least 2n-1. The explicit one-odd-plus-(n-1)-even-pairs construction attains it. The even-sector argument analogously forces degree at least 2n.

### Determinant/permanent collision and root hostile

The character/cycle program directly enumerates all 9! permutations of each of three matrices. It counts Hamiltonian paths, determinant/permanent coefficients, and odd-cycle packings independently of the workers' rooted-cycle and subset routines.

For masks 29759988540 and 18875161993, the common polynomials in ascending coefficients are

- per(xI+A): [384,692,629,380,154,68,28,0,0,1].
- det(xI+A): [20,-82,167,-192,154,-68,28,0,0,1].

Their odd-cycle packing vectors are [1,771,303,11] and [1,772,304,10], giving H=2843 and 2841. The full realized difference

C_T(z;u,v)-C_U(z;u,v) = z^9(u^3-u)+z^8(v^2-u^2)

is consistent with the independently checked permutations and the workers' full length-profile tables. Mixed (3,6) and (4,5) changes cancel because both map to the same z^9 uv term. The first even response agrees at 1567; the second divided by 2 is 96 versus 95.

The root fixture has I(z)=1+49z+7z^2+z^3, H=135, and discriminant -348172. Its single odd-cycle counts by length 3,5,7,9 are 8,20,18,3.

### Infinite reversed-arc family

`check_triangle_family.py` uses an independent minimum-root subset DP for simple-cycle counts, a disjoint-packing DP, and a separate unrestricted Hamiltonian-path DP. It directly checks m=2,...,5 (n=6,...,15):

- m=2: [1,6,1], H=17.
- m=3: [1,49,7,1], H=135.
- m=4: [1,452,98,8,1], H=1377.
- m=5: [1,4509,1354,148,9,1], H=15795.

The all-m proof is independent of these finite controls. With P(t)=t+t^2+t^3, the new-cycle contribution is

z Odd_t[P(t)^2(1+z+3P(t))^(m-2)],

which gives the stated formula by evaluation at t=+1 and -1. The normalized top-coefficient Newton difference is

-4(21m^2-48m+4)/(m^2(m-1)) < 0 for every m>=3.

This rules out all-real-rootedness, and hence a finite Hermitian positive-semidefinite free-mode determinant. It does not exclude signed, interacting, non-Hermitian, or nilpotent determinant-ratio realizations. H_2 is exceptional, so the two-term Hamiltonian recurrence correctly starts at m=5.

The minimal root-failure order and the minimal determinant/permanent collision order are distinct conclusions. Their shared value nine does not establish an implication between the two properties.

### Signed necklaces and PBW

`check_signed_necklaces.py` derives the strong alphabet and ordinary Euler exponents from this audit's independent tournament character counts, then enumerates all 7,412 nonempty weighted SCC words through weight eight. It checks the primitive-root sign criterion, explicit signed rotation averages, signed logarithmic Burnside series, and the PBW Hilbert product. The signed necklace counts are

1, 0, 2, 2, 9, 43, 399, 6415.

The free-Lie-super dimension sequence is

1, 1, 1, 2, 8, 44, 398, 6413.

The general super-Lyndon and PBW inputs remain explicitly cited results, not new proofs in this audit. Their stated locations were checked directly in the [primary Bokut–Kang–Lee–Malcolmson paper](https://arxiv.org/pdf/math/9809024): the section-2 definition includes odd Lyndon squares, Theorem 2.2 supplies the basis, and Proposition 2.11 supplies PBW. The inherited curl-stabilizer classification is used as inherited input; this audit checked the resulting determinant character and index-two orbit argument rather than repeating the earlier classification proof.

## Reproduction and provenance

From this directory:

    g++ -O3 -std=c++17 check_graph_signs.cpp -o /tmp/session37_graph_signs
    /tmp/session37_graph_signs > graph_signs.json
    python3 check_characters_and_cycles.py
    python3 check_components.py
    python3 check_squarezero.py
    python3 check_triangle_family.py
    python3 check_signed_necklaces.py

Every Python assertion is an explicit exception check and remains active with -O. Normal and optimized outputs match. `INPUT_SHA256SUMS` records audited source-note versions; `SHA256SUMS` records this audit's deliverables. The report does not claim a physical spin-statistics theorem, quantum dynamics, arithmetic transcendence consequence, or operation-preserving equivalence between the different algebras.
