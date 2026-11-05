# An exact self-converse path-orbit mass and a constructive search universe

Separate Session197 theorem. This refines the [S196 path-orbit mass identity](https://github.com/eliottcassidy2000/math/blob/cf3b44b9a442c3cff0f28cab924e73c4c00bb00c/session196/symmetry/PROOF.md) to self-converse tournaments. It provides an exact enumeration universe and fibre weight, not a spectrum coverage claim. No literature-priority claim is made.

Let T be a nonempty tournament and G=Aut(T). A permutation tau is an anti-involution if tau squared is the identity and it reverses every arc. For such tau let f_tau(T) count Hamilton words fixed by

    J_tau(v1,...,vn)=(tau(vn),...,tau(v1)),

and let C_G(tau)={g in G:g tau=tau g} be its centralizer.

## 1. Exact mass theorem

**Theorem.** Every self-converse tournament admits an anti-involution. The quantity

    c(T)=f_tau(T)/|C_G(tau)|                            (1)

is a positive odd integer, independent of the chosen anti-involution. For every n>=1,

    sum_[T] self-converse, |T|=n c(T)
       = 2^floor((n-1)^2/4).                           (2)

In particular, the number of ordinary self-converse tournament classes is even for every n>=3.

For the optional empty-tournament convention H(empty)=1, there is one empty word and a trivial automorphism group. Taking tau to be the identity gives c=1, and equation(2) also holds at n=0. This case does not use an order-two argument.

**Anti-involutions and conjugacy.** The n=1 case is immediate with tau the identity. At n>=2 an anti-automorphism has even order, since an odd power still reverses arcs. Its square lies in G, whose order is odd: G acts freely on the odd set of Hamilton words, as proved in S196. Thus the anti-automorphism has order2m with m odd, and its mth power is an anti-involution.

Any two anti-involutions tau,sigma are conjugate by an element of G. Put g=tau sigma in G, of odd order m. Then tau g tau=g^(-1), and

    g^((m-1)/2) tau g^(-(m-1)/2)=sigma.

This proves conjugacy without needing a classification of tournament automorphism groups. Consequently their fixed-word counts and centralizer orders agree.

**Integrality and oddness.** J_tau is an involution of the Hamilton words, so f_tau is odd and positive because H is odd. The centralizer acts freely on that fixed-word set: commuting with tau preserves J_tau-fixedness, and an automorphism fixing a full vertex word is the identity. Therefore |C_G(tau)| divides f_tau. Its order is odd, so the quotient in (1) is a positive odd integer.

## 2. Direct normalization proof and explicit construction

Normalize a distinguished Hamilton word by naming its entries0,1,...,n-1. The consecutive arcs are forced forward. Reversal of this named path is the fixed permutation r(i)=n-1-i.

To make r an anti-automorphism, each remaining forward-orientation bit on a pair i<j must equal the bit on its reflected pair

    (i,j) <-> (n-1-j,n-1-i).

There is no bit complementation: applying both converse and reversal cancels it. Consecutive pairs reflect to consecutive pairs, so their forced arcs are consistent. All other pair-orbits can be assigned independently.

There are d=binom(n-1,2) unconstrained pairs and floor((n-1)/2) reflected fixed pairs among them. Hence their number of orbits is

    [d+floor((n-1)/2)]/2 = floor((n-1)^2/4).

Choosing one free bit per orbit produces exactly the right side of (2) distinct normalized tournaments. Every output has an actual Hamilton word and an explicit anti-involution. This is a complete independent-bit generation universe for self-converse path-pointed objects, not a uniform sampler of ordinary classes.

Fix an ordinary self-converse class T and one anti-involution tau. Each tau-symmetric Hamilton word normalizes to one of these objects. Two such words give the same normalized tournament precisely when an automorphism carries one to the other. If w'=g w and both are tau-symmetric, then tau g and g tau agree on every named vertex of w, so g commutes with tau. Thus the fibre inside T consists of centralizer-orbits and has size c(T).

Conversely any normalized object isomorphic to T pulls its fixed reversal r back to an anti-involution tau' of T. The conjugacy proved above lets the isomorphism be adjusted by an automorphism so that it pulls r back to tau. Its Hamilton word is then tau-symmetric. Thus no normalized objects of this class are missed. Summing the exact fibres proves (2).

Since every summand is odd, (2) also gives the asserted parity of the number of self-converse classes. At n=1,2 the exponent is zero, as it must be.

## 3. Equivalent orbit interpretation and concrete examples

The usual path-orbit quotient H(T)/|G| is odd. Converse, followed by any anti-isomorphism back to T, induces an involution on these G-orbits. Its fixed-orbit count is exactly c(T): on every invariant free G-orbit there are |G| pairs consisting of a word and its unique reversing anti-involution; all anti-involutions form one G-conjugacy class. This agrees with the direct normalization proof above.

For the self-converse eleven-vertex7p example in `SELF_CONVERSE.md`, G=C7 and its anti-involution inverts G, so its centralizer is trivial. The fixed-word count is111, and that class contributes111 to (2), while contributing8069 to the unrestricted path-orbit mass.

For the module-indecomposable fourteen-vertex7p example, the centralizer is again trivial and the contribution is6111, versus2163467 in the unrestricted mass. These different weights are not ordinary-class multiplicities and do not imply every prime occurs as a quotient.

The construction can therefore be used as a complete self-converse search universe: every such class appears, its exact multiplicity is known, and there are only2^floor((n-1)^2/4) normalized presentations. This count alone gives no arithmetic coverage guarantee. Conversely, self-converse candidates for prime or7p targets need not be discarded merely because the first symmetric4151 cores in S196 were non-self-converse.

## 4. Exact controls

`verify_converse_mass.py` constructs the complete normalized universe through order6, canonicalizes every output under all vertex permutations, and checks every individual class fibre against f_tau/centralizer. It independently enumerates the ordinary self-converse classes among all labelled tournaments through order5, so the constructed universe's completeness is tested without assuming it from the proof. All Hamilton words, automorphisms, anti-involutions and centralizers in that control universe are explicitly enumerated.

Run normally and under `python3 -O`, comparing receipts. Finite checks control the all-order proof; they do not replace its normalization and group-action arguments.
