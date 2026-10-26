# The complete quotient response is exactly unsigned triangle curl

Status: PROVED elementary walk identities and full-support detector; CITED classification of labelled cyclic-triangle structures; FINITE-EXACT checks and independent audit. This completes, rather than conjecturally strengthens, the full-support/full-series issue left cautious in the first draft of RESULT.md.

## 1. Main equivalence

For tournaments `T,U` on the SAME labelled vertex set `[n]`, the following are equivalent:

1. Their labelled sets of cyclic triples agree.
2. Their absolute signed triangle-curl fields agree: `|c_T(i,j,k)|=|c_U(i,j,k)|`, where `c_T(i,j,k)=S_ij+S_jk+S_ki` and `S=A-A^T`.
3. A sequence of homogeneous-interval reversals transforms `T` into `U`.
4. Their full ordinary commutative walk series `W_T(x)` and `W_U(x)` agree.
5. Their full-support parts `W_T^+(x)` and `W_U^+(x)` agree.
6. `H(T[S_1,...,S_n])=H(U[S_1,...,S_n])` for every tuple of nonempty tournament factors.
7. The same equality holds when every factor is transitive.

A homogeneous interval is a vertex subset whose vertices have identical adjacency to every exterior vertex. Reversing it means reversing all internal arcs and leaving all other arcs unchanged. The whole vertex set is allowed, giving global converse.

The equivalence of1 and3 is CITED, not new: Boussaïri, Ille, Lopez and Thomassé, “The C3-structure of the tournaments,” author manuscript dated June14,2005, Theorem2 (PDF pages2–3; proof pages8–9). Its Corollary1 (PDF page7) also says that if one tournament is indecomposable, equality of cyclic-triple sets forces equality or global converse. Primary source inspected and independently checked: https://perso.ens-lyon.fr/stephan.thomasse/liste/bilt.pdf .

The connection of this structural theorem to complete Hamiltonian substitution response is proved below. No new priority claim is made for the cited classification or the older substitution kernel.

## 2. Ordinary walk series composes exactly

For any nonempty tournament `T`, let

    W_T(x)=1^T X(I-A_T X)^(-1)1,

where `X=diag(x_v)` and ordinary commuting variables may repeat. This sums directed walks by vertex content, including singleton walks.

**Lemma1 (walk substitution).** For tournament substitution `T=Q[S_1,...,S_q]`,

    W_T(x)=W_Q(W_(S_1)(x|S_1),...,W_(S_q)(x|S_q)).       (1)

Proof. Split any walk in `T` at every change of module. Each resulting maximal run is a nonempty directed walk in one factor. The module word is a directed quotient walk, because its consecutive letters are different and connected by a quotient arc. Conversely, independently replacing every occurrence in a quotient walk by a factor walk produces a unique walk in `T`. Repeated occurrences of a vertex are allowed on both sides. Constant terms are zero, so formal composition is well-defined. QED.

Reversing every directed walk gives the content-preserving identity

    W_(S^rev)=W_S.                                      (2)

Contracting a homogeneous interval makes it one substitution block. Equations (1) and (2) therefore show that interval reversal preserves the ENTIRE `W`, not just the Hamiltonian coefficient or one fixed blowup.

Combined with the cited theorem, this proves1=>3=>4. Conversely,

    [x_i x_j x_k] W_T = H(T[{i,j,k}])
                       =3 if the triple is cyclic, else1, (3)

so4=>1. Also1<=>2 since absolute triangle curl is3 precisely on cyclic triples and1 on transitive triples.

The equivalence of5,6,7 is the nonempty-factor Stirling inversion in RESULT.md, TheoremC. The missing implication is5=>1. It does NOT require deleting slots or assuming anything about rational denominator factorization.

## 3. A finite full-support coefficient detects every triangle

Fix a triple `R={i,j,k}` in an `n>=3` vertex tournament. Define the content vector

    b_v=n-1 for v in R,
    b_v=1   for v outside R.

Its total degree is `4n-6`.

**Lemma2 (positive detector).**

    [x^b] W_T^+ > 0  iff  T[R] is cyclic.                (4)

Proof, cyclic direction. Start with any Hamiltonian path of `T`, whose existence is Rédei's theorem. At its occurrence of one triangle vertex, insert `n-2` complete directed turns around the triangle, returning each time to that vertex before continuing. Each turn adds one visit to each of the three triangle vertices and no visit outside the triple. The resulting walk has exactly the prescribed content.

Proof, transitive direction. Any walk with that content visits the exterior vertices exactly `n-3` times. Hence its visits inside `R` split into at most `n-2` maximal consecutive runs. Inside a transitive triple each run is strictly increasing in its total order, and visits any vertex at most once. Thus each triangle vertex appears at most `n-2` times in the entire walk, contradicting the prescribed count `n-1`. QED.

Therefore `W^+` determines every cyclic triple. This closes all implications in Section1, and shows that the cautious first-draft distinction between full `W` and `W^+` disappears on the restricted class of tournaments.

The coefficient can itself be recovered from finite Hamiltonian observations by the explicit multivariate signed-Stirling inverse. For this detector, exterior block sizes are fixed at1 and the three marked block sizes range from1 through `n-1`. Thus at most `(n-1)^3` transitive-factor probes per triple suffice. They involve at most `4n-6` vertices. This is an existence/reconstruction bound, not an optimized measurement scheme or polynomial-time counting claim.

## 4. SCC factorization and the exact information loss

If `T=S_1▷...▷S_r` is its ordered chain of strong components, a directed walk may visit an increasing subsequence of SCCs, with one nonempty run in each visited SCC. Thus

    1+W_T = product_(i=1)^r (1+W_(S_i)),                 (5)

and on disjoint labelled vertex sets,

    W_T^+ = product_(i=1)^r W_(S_i)^+.                   (6)

Both expressions are commutative: every permutation of the SCC blocks preserves the entire walk series and every vertexwise-substitution Hamiltonian count. This is larger than preserving signed curl under the specific cyclic SCC rotations in the root's curl classification.

These factorization identities are direct consequences of the walk decomposition. They do not assume uniqueness of rational factors or reconstruct the SCC partition from an arbitrary formal series. The cyclic-triangle classification supplies the necessary structural information separately.

The path-colour polynomial likewise satisfies

    Q_T(t)=product_i Q_(S_i)(t),

by inherited THM-3166. The scalar-profile carrier forgets even more, because two strong tournaments can already have identical `Q` and different cyclic-triangle structures (THM-3235 twins).

## 5. Sharp labelled hostile and a nontrivial strong hostile

### 5.1 Smallest labelled failure of signed curl up to converse

On four fixed labels, compare the transitive orders

    0<1<2<3,               0<1<3<2.

Their entire walk series is identical:

    W=product_(v=0)^3(1+x_v)-1.

For the increasing triples `012,013,023,123`, their signed curl vectors are

    (1,1,1,1),             (1,1,-1,-1).

These are neither equal nor global negatives. Thus even complete quotient response does not determine signed curl up to converse. Order4 is minimal for this labelled statement: at order3 there is only one triple, and equality of its cyclic/transitive type fixes the curl up to a global sign; orders1,2 have no triangle coordinates.

These transitive tournaments are, of course, isomorphic after forgetting labels. This is not an unlabelled-isomorphism counterexample.

### 5.2 A smallest-order strong labelled mechanism

Take `C3[P_2,1,1]`, with the two-vertex module on labels0,1, and reverse just that internal arc. Both tournaments are strong. Their curl vectors are

    (1,1,3,3),             (-1,-1,3,3),

again neither equal nor global negatives, while their entire `W` is identical by the interval-reversal lemma. This also occurs at order4, and the tournaments are isomorphic by exchanging the two module labels.

### 5.3 Strong, nonisomorphic and nonconverse example at order7

Let `A=1▷C3` (order4), `B=P_2`, and form

    T=C3[A,B,1],            U=C3[A^rev,B,1].              (7)

Interval reversal proves `W_T=W_U`, so all vertexwise-substitution Hamiltonian counts agree.

The maximal proper strong modules in this displayed cyclic substitution have distinct sizes4,2,1. Any isomorphism must preserve their sizes. Since `A` has a source but no sink and `A^rev` has a sink but no source, `T` and `U` are not isomorphic. A converse-isomorphism would also have to preserve these three distinct block sizes, but reversing the outer triangle changes all its block arcs, so it is impossible. Both tournaments are strong because the quotient is a directed triangle.

The verifier independently checks nonisomorphism and nonconverse by all vertex permutations with score pruning, and checks concrete equal-profile and equal-blowup counts. No claim is made that7 is the minimum unlabelled nonconverse order.

### Strongest reconstruction survivor

If `T` is indecomposable, the common response state determines it exactly up to global converse, by the cited Corollary1. Strong connectivity alone is insufficient, as (7) demonstrates. The missing hypothesis is absence of nontrivial homogeneous intervals, not an additional numerical moment.

## 6. Relationship with the inherited forest polynomial

The exact normalization remains

    E_T(u)=sum_even L 2^(-L)(n-L)! N_L(T) u^(n-1-L),

where `N_L` is the signed directed-linear-forest sum of THM-4757. It is the same module-profile information in another basis. The cyclic-triangle hypergraph is a finer labelled quotient-response state. Forgetting labels and collapsing to the path-colour/forest coefficients can merge strongly connected, nonconverse quotients, even though those coefficients are complete for placing the tournament inside one uniform module slot.

The genuinely recovered hierarchy is therefore:

    signed curl / labelled orientation information
      -> unsigned curl = cyclic-triangle hypergraph
      = full walk series = full vertexwise-substitution response
      -> complete path-cover / path-colour / forest profile
      -> one Hamiltonian count.

The first arrow loses interval-orientation and SCC-order information; the second loses vertex-addressed threading information; the last loses run counts. None of these losses is repaired by calling the resulting object a “complete invariant” without naming its interface.

## 7. Verification and source discipline

`verify_clans.py` checks every labelled tournament through order5. It checks full walk coefficients through total degree7 within every cyclic-triangle class, and10,504 instances of the exact full-support detector (4). It also checks the labelled hostile examples, the order7 nonisomorphism/converse obstruction, equal profiles, and selected blowups (including the uniform twofold, order14).

The independent auditor used a separate implementation:

- every1,099 labelled tournament through order5;
-10,504 detector instances;
-3,830 homogeneous-interval inversions;
- all full walk coefficients with every exponent in0..2 within each C3 class;
-560 further detector checks on112 deterministic order6..8 samples.

Independent paths: `session36-audit/independent_c3_walks.py`, `session36-audit/c3_walk_receipt.json`.

The all-order equivalence uses the cited structural theorem and the elementary proofs here. Finite checks audit its application and do not replace that theorem. Exact operation scope, label dependence and source credit are essential.
