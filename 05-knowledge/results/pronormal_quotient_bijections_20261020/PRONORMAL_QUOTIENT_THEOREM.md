# Prime-local symmetry beyond p-groups: pronormal bases and normalizer quotients

**Status: PROVED below; independent review requested before promotion.** Exact controls are supplied separately. These results concern isomorphism types of tournaments and **sign-even** graphs: every graph automorphism must have trivial edge-reversal sign. This is not the Eulerian-degree condition.

## 1. Current starting point and scope

The repaired HYP-9340 already proves equality for every odd p-group, including every C3 action on six vertices. Its universal odd-order version is false: the C15 action with orbit sizes (3,5,1,1) at n=10 has 32 tournament types and 31 sign-even graph types.

Dependencies, in truth-precedence order:

- `05-knowledge/results/c15_symmetry_counterexample_20261020/C15_ORDER_TEN_COUNTEREXAMPLE.md` and its independent review;
- `05-knowledge/results/prime_local_walsh_bijection_20261020/SYLOW_ACTION_ADDENDUM.md`;
- the weighted Walsh matching theorem in that packet's `PRIME_VALUATION_THEOREM.md`;
- THM-4721, `01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md`, Theorem R and Proposition Q.

The new point is a second localization: first fix an odd subgroup K, then work at an odd prime q in its normalizer quotient. A group-theoretic descent condition makes the resulting relative statement a statement about ordinary isomorphism types.

Notation: G=S_n, N=N_G(K), Gamma=N/K. For an object x fixed by K, let A_x=Aut(x), and let Gamma_x=(N intersect A_x)/K. A K-fixed graph is N-even when its sign character is trivial on N intersect A_x.

## 2. General relative fixed-cube theorem

**Theorem 1.** Let K<=S_n have odd order, and let q be any odd prime. There is an explicitly computable bijection between

1. N-orbits of K-fixed tournaments; and
2. N-orbits of N-even K-fixed graphs,

preserving the Gamma-conjugacy class of a Sylow q-subgroup of Gamma_x. One such bijection preserves simultaneously every predicate

    Gamma_x contains a Gamma-conjugate of U,

where U ranges over all q-subgroups of Gamma, including every Boolean combination of these predicates.

### Proof: the correct Fourier dual of the fixed cube

Let O_1,...,O_d be the K-orbits on unordered vertex pairs. Each |O_i| is odd. No element of K reverses an edge, so choosing one orientation per orbit gives all 2^d K-fixed tournaments. A K-fixed graph is the union of any subset of these orbits, also giving 2^d objects.

For a K-invariant graph X use the restricted character

    phi_X(T)=product over e in X of T(e),

with T(e)=+1 on the increasing-label orientation and -1 on the reverse. On the fixed tournament cube, reversing orbit O_i changes this product by (-1)^|O_i|=-1 exactly when O_i belongs to X. Thus these 2^d restricted characters are all the characters of the d-cube, up to individual signs. In particular their matrix has inverse 2^(-d) times its transpose.

The odd orbit-size hypothesis is essential here. Simply restricting arbitrary characters to fixed points can identify characters; restricting to K-invariant graphs works because all orbit sizes are odd.

N permutes the characters with the usual edge-reversal sign cocycle. Elements of K act trivially on the fixed cube and have trivial graph sign, so the signed action factors through Gamma. Taking Gamma-invariants over R=Z[1/2] gives two R-bases: tournament-orbit indicators and signed graph-orbit sums for N-even graph orbits. For one representative X of such an orbit set

    Phi_X^N(T)=sum over Y in N.X of epsilon(Y) phi_Y(T).

The transporter sign epsilon is well-defined precisely because X is N-even. The square integer matrix M_K=(Phi_X^N(T)) is invertible over R. Character orthogonality grouped into Gamma-orbits gives

    M_K diag(|Gamma_X|) M_K^t = 2^d diag(|Gamma_T|).       (1)

The inherited weighted-pairing theorem therefore supplies a perfect matching supported on q-unit entries and preserving v_q(|Gamma_x|).

Choose a Sylow q-subgroup Q of Gamma_T. The signed summands in Phi_X^N(T) are constant on Q-orbits of graph copies, by the same two-cocycle cancellation as in the inherited Sylow-action lemma. Non-singleton Q-orbits have q-divisible size. A q-unit coefficient therefore forces a Q-fixed copy Y of X. Hence Q<=Gamma_Y. Equality of q-valuations makes Q a Sylow q-subgroup of Gamma_Y. Transporting the copy back to X proves conjugacy in Gamma. The subgroup-presence and Boolean-combination claims follow from Sylow containment. QED.

This proof uses the original restricted signed coefficients. It does not infer a set bijection from a linear isomorphism without the q-unit matching and localization steps.

## 3. Pronormality is a sufficient descent condition

A subgroup K is **pronormal in G** if K and gKg^(-1) are conjugate inside the group they generate, for every g in G. This is a standard notion; for terminology see Revin and Vdovin, *The existence of pronormal pi-Hall subgroups in E_pi-groups*, [arXiv:1504.02834v2](https://arxiv.org/abs/1504.02834v2). No classification theorem or Hall-subgroup existence theorem from that paper is used below.

**Lemma 2.** If K is odd and pronormal in S_n, then:

(a) Every full isomorphism type meeting Fix(K) is one N-orbit.

(b) A K-fixed graph is N-even if and only if it is globally sign-even.

**Proof of (a).** Suppose x'=g.x and both objects are K-fixed. Then K and gKg^(-1) lie in A=Aut(x'). Pronormality gives a conjugator a inside their generated subgroup, hence inside A, with a g K g^(-1) a^(-1)=K. Thus ag lies in N and (ag).x=x'.

**Proof of (b).** Only the nontrivial direction needs proof. If a K-fixed graph X is odd, set A=Aut(X), A^+=ker(sign_X), and choose g in A with sign_X(g)=-1. The odd group K lies in A^+, as does gKg^(-1). A pronormal conjugator a as above belongs to their generated subgroup, hence to A^+. Therefore ag lies in N intersect A and has sign -1. The graph is not N-even. QED.

In particular the relative-even theorem immediately yields

    tau_T(K)=tau_E(K) for every odd pronormal K<=S_n.     (2)

This sufficient class is different from the inherited class of all odd p-groups: a p-subgroup need not be pronormal, and the p-group theorem did not require it.

## 4. The mixed-prime survivor

**Theorem 3.** Let K be an odd pronormal subgroup of S_n and q an odd prime. On the isomorphism types admitting K, there is one bijection from tournaments to sign-even graphs preserving every predicate

    Aut(x) contains an S_n-conjugate of H,

simultaneously for every subgroup H with

    K <= H <= N_(S_n)(K), and H/K a q-group.             (3)

Consequently tau_T(H)=tau_E(H) for each H in (3). No requirement that H itself be pronormal is made.

**Proof.** Lemma 2 identifies both relative object sets in Theorem 1 with ordinary types. It remains to verify the predicate, rather than silently replace ambient conjugacy by normalizer conjugacy. Fix a K-invariant representative x. If gHg^(-1)<=A_x, then both K and gKg^(-1) lie in A_x. Pronormality supplies a in A_x for which ag normalizes K. Thus (ag)H(ag)^(-1)<=N intersect A_x. Conversely such a normalizer containment is plainly an ambient containment. Since subgroups containing K correspond bijectively to subgroups of Gamma, the desired predicate is exactly Gamma-conjugate containment of H/K in Gamma_x. Apply Theorem 1. QED.

The preserved datum is the conjugacy class of a Sylow q-subgroup of the **normalizer quotient stabilizer**, not the entire automorphism group or all of its Sylow subgroups. The datum is well-defined on an ordinary type: Lemma 2 makes all choices of its K-copy conjugate under its automorphism group.

### Ambient-Sylow specialization

Every Sylow subgroup P of S_n is pronormal: P and gPg^(-1) are Sylow subgroups of their generated subgroup, so Sylow conjugacy supplies the required conjugator. Hence Theorem 3 applies whenever P is an ambient Sylow p-subgroup for an odd prime p and H/P is a q-group inside N(P)/P. This proves a family of genuinely mixed-prime equalities, even when p!=q and H is not an ambient Sylow subgroup.

Taking K=1 recovers the inherited all-q-group theorem. Taking H=K recovers the odd-pronormal statement. The useful new part is a nontrivial base and a nontrivial quotient.

### Iteration through pronormal normalizer quotients

**Lemma 4.** Suppose K is pronormal in G, K<=H<=N_G(K), and H/K is pronormal in N_G(K)/K. Then H is pronormal in G.

**Proof.** Fix g in G and L=<H,gHg^(-1)>. Pronormality of K supplies a in <K,gKg^(-1)><=L such that ag normalizes K. Put J=(ag)H(ag)^(-1), so H,J<=N_G(K), K<=H intersect J, and J<=L. In the quotient N_G(K)/K, pronormality of H/K supplies a conjugator in <H/K,J/K>. Since <H,J> contains K, this conjugator lifts to b in <H,J><=L. Choosing its direction correctly gives bJb^(-1)=H: equality follows because both subgroups contain K and have equal quotient images. Then ba in L conjugates gHg^(-1) to H. QED.

Thus one can start with K=1 and repeatedly replace K by the full preimage of an odd-order pronormal subgroup of N_G(K)/K. Each result remains odd and pronormal. In particular, taking a Sylow subgroup for an odd prime at each stage creates an explicit family of mixed-prime normalizer towers to which (2) applies, and Theorem 3 allows a final arbitrary q-subgroup extension. The construction requires an ambient normalizer quotient at each step; an arbitrary abstract Sylow tower inside H is not enough.

## 5. Why C15 succeeds at n=8 and n=9 but fails at n=10

Let H act as a 5-cycle and a disjoint 3-cycle, with all other vertices fixed.

- At n=8 and n=9 the subgroup P=C5 is a Sylow 5-subgroup of S_n. Its normalizer quotient is C4 x S_(n-5), and H/P is a C3 subgroup. Theorem 3 therefore proves the inherited counts agree: 2=2 at n=8 and 8=8 at n=9.
- At n=10 a Sylow 5-subgroup of S_10 has order 25. The displayed C5 is no longer an ambient Sylow subgroup. It is also not pronormal: inside C5 x C5 on two disjoint 5-sets, a permutation swapping the sets conjugates the first factor to the second, while that abelian generated group cannot conjugate one factor into the other.
- Neither choosing C3 as the base nor using normality of C5 in C15 repairs this: the single C3 is not an ambient Sylow 3-subgroup, and abstract normality inside H is weaker than the required ambient descent condition.

The known K_(5,5) obstruction makes this loss concrete. Its two parts admit two distinct C5 factors; the H-fixed presentation has a 3-cycle plus two fixed vertices on one part and the original 5-cycle on the other. The new odd part swap escapes the chosen normalizer. Thus the boundary is tied to an actual conjugacy and sign obstruction, not just a failed divisibility computation.

The theorem does not claim the ambient-Sylow condition is necessary; arbitrary p-groups can satisfy equality by the separate inherited theorem, and other nonpronormal mixed groups may satisfy it for other reasons.

## 6. A prime-power base with nontrivial mixed-prime semidirect action

At n=15 take P=C5 x C5 x C5 on three disjoint five-vertex blocks. Since v_5(15!)=3, P is an ambient Sylow 5-subgroup. Its normalizer quotient is C4 wr S3. Let H=P semidirect C3, where C3 cyclically permutes the three blocks. Thus H is the imprimitive C5 wr C3, of order 375, and H/P=C3. Theorem 3 proves tau_T(H)=tau_E(H).

Here the value is exactly one on each side. P-invariant tournaments have uniform cross-block arcs, and each five-block is a circulant tournament, with all four choices isomorphic under multipliers. Invariance under the block 3-cycle forces the quotient to be a directed triangle, giving the unique type C3[C5,C5,C5]. For graphs, any nonempty internal circulant graph on five vertices has an odd automorphism extending by the identity (the options are C5, its complement, and K5); all blocks must therefore be independent. H invariance leaves only an empty or complete three-block quotient. The empty graph is sign-even. The complete quotient gives K_(5,5,5), whose transposition of two parts reverses 25 cross edges and is odd. Thus only the empty graph survives.

The base P sector itself contains two tournament types and two sign-even graph types, corresponding to the two three-vertex quotient types on either side. This is a concrete nonabelian mixed-prime consequence with a prime-power-order normal base.

## 7. Exact controls and source boundaries

### The user's six-vertex C3 case

`n6_c3_certificate.py` independently constructs all fixed labelled objects, tests full S6 isomorphism and full automorphism signs, and includes all 720 permutations. It uses no inherited atlas and no graph-isomorphism library.

| embedded subgroup | fixed tournaments/graphs each | tournament types | sign-even graph types | normalizer-orbit count each |
|---|---:|---:|---:|---:|
| one 3-cycle, three fixed vertices | 128 | 11 | 11 | 12 |
| two synchronized 3-cycles | 32 | 3 | 3 | 4 |
| independent C3 x C3 | 8 | 1 | 1 | 1 |

The first two predicates intersect in exactly the C3 x C3 type on each side. The exact Sylow-C3 strata have sizes ten and two; one further type has Sylow C3 x C3. For the single-cycle normalizer census, the extra globally odd graph is K_(3,3). For the double-cycle census, it is the disjoint union of two triangles. On each tournament side the C3 x C3 type contributes two normalizer orbits, exactly compensating the extra odd graph. Thus n=6 illustrates cancellation beyond the elementary ambient-Sylow proof, rather than a still-open case.

### A genuinely mixed-prime quotient census at n=10

`c7_normalizer_quotient.py` takes K=C7 with three fixed roots. This is an ambient Sylow 7-subgroup; Gamma=C6 x S3 has order 36. It enumerates all 512 objects on the 9-dimensional fixed cube, builds the entire 24-by-24 signed Walsh matrix, checks weighted orthogonality exactly, computes determinant +/-2^88, and finds a q=3-unit matching. It checks the whole quotient Sylow action for every matched pair, not only the subgroup order.

Let the multiplier a:x->2x on F7 and root cycle b=(7 8 9) generate the two C3 factors in the quotient. The same matching simultaneously gives:

| group above K | quotient subgroup | type count each |
|---|---|---:|
| C7 | 1 | 24 |
| F21 fixing all three roots | <a> | 12 |
| C7 x C3 | <b> | 4 |
| diagonal F21 on vertex orbits 7+3 | <ab> | 2 |
| F21 x C3 | <a,b> | 2 |

Here F21=C7 semidirect C3 with nontrivial multiplier action. Its diagonal action on the three roots is faithful only after quotienting the C7 kernel. The two diagonal lines are Gamma-conjugate, because a transposition of roots inverts b. The table covers all quotient-C3 subgroup conjugacy classes and their full elementary-abelian join.

The diagonal F21 is itself **not pronormal in S10**, so this example genuinely needs the extension theorem rather than just (2). Conjugating by a transposition of roots sends its quotient line <ab> to <ab^(-1)>. The two groups generate L=F21 x C3. Modulo the normal C7 kernel, L/C7=C3 x C3 is abelian, and its two distinct diagonal subgroups cannot be conjugated by L/C7. Hence the two original H groups cannot be conjugated inside their generated group L. Nevertheless their common ambient symmetry predicate has the proved count 2=2.

The complete joint histogram for (<a>,<b>,<ab>) presence is (false,false,false):10; (false,true,false):2; (true,false,false):10; (true,true,true):2 on both sides. These are all three line predicates for one quotient prime, not simultaneous independent Sylow predicates for the full automorphism group.

Both scripts use explicit runtime checks that remain active with python -O. Standard and optimized JSON outputs agree byte for byte. The n=10 finite computation is a normalizer-quotient census; its interpretation as global types and global evenness uses the proved ambient-Sylow descent. It is not advertised as an independent full S10 census.

## 8. Mechanism, boundary, and next question

Map: Fourier character transport on a K-fixed Boolean cube, followed by a q-unit matching on normalizer orbits. Preserved: quotient Sylow-q action, hence all extension predicates (3). Lost: joint information from multiple quotient primes, whole automorphism groups, labelled naturality, and arbitrary embeddings of a mixed-prime H. Sidecar: the base subgroup K and its ambient conjugacy-control property. Cheapest hostile: the C15 action at n=10, where the base C5 fails pronormality and the K_(5,5) sign escapes.

This does not revive all odd-order H, and it does not give one global bijection preserving all odd-prime predicates. The inherited five-vertex joint-prime obstruction remains decisive. A useful next question is to characterize the weakest object-relative fusion/sign-control conditions that permit this descent, or to obtain a comparably explicit mixed-prime statement without a distinguished pronormal base.

No literature-priority claim is made. The standard definition of pronormality is cited; all descent and fixed-cube arguments used here are proved above from the inherited Walsh mechanism and elementary finite-group facts.
