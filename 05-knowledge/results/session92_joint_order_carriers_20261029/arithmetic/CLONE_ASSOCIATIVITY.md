# Source-preserving divisor clones: the LCA correction that survives rebracketing

**Status:** PROVED elementary finite-measure statements; FINITE-EXACT controls. No novelty-priority claim. No theorem about aliquot orbit termination, prime untouchables, or all tournament carriers.

## Result in brief

A quotient can preserve every scalar divisor sum and still fail to compose. For three copies of a prime p, the native divisor p has three named sources `(p,1,1)`, `(1,p,1)`, `(1,1,p)`. Locally uniform cloning through `((pp)p)` assigns probabilities `(1/4,1/4,1/2)`; through `(p(pp))` it assigns `(1/2,1/4,1/4)`. Neither is the globally uniform `(1/3,1/3,1/3)`.

The exact repair is to retain the **number of descendant realizations of each intermediate divisor**, not merely its native divisor weight. These multiplicities give an associative, source-preserving message-passing rule. They are the arithmetic counterpart of keeping joint port compatibility rather than independently resampling each cursor.

## 1. Provenance and recovered correction

Source snapshot: `b07ec51b204e5c0e09e52eed206cb2bd3cc29ebc`.

- `05-knowledge/results/session88_symmetry_certificates_20261028/operators/TRANSFER.md`, Theorem 5.3 and “Coprimality hostile and clone weights”: for coprime A,B, native divisor weights form a product measure. After deleting the top pair, the flag angle has squared cosine `AB/sigma(AB)`. For `6×6`, naive mass is144 rather than 91, proper mass 108 rather than 55. Divisor multiplication produces clone fibres, and silently resetting mass invalidates the product proof.
- `05-knowledge/results/session85_arithmetic_codes_free_bits_20261028/LINEAR_GLUING_INFORMATION.md`, sections2–5: exact shared-source joins, missing separator constraints, affine messages, and LCA rank ownership. Tree structure alone cannot recreate correlations that were thrown away.
- Current user’s strong-child tournament witnesses: separate cursor responses agree, but assembled Hamilton counts differ. Those witnesses motivate the question; this note does not re-prove their numerical counts or claim the divisor model is isomorphic to them.

Closest inherited mechanism: corrected clone pushforward. Corrected near miss: mass correction restores native sums, not product independence. Underused sidecar: an integer extension-count table at every internal node. Wildcard made exact: changing tree shape exposes hidden nonuniform source weights even when every total stays correct.

## 2. Two-factor clone repair and its exact boundary

Let D(A) be the positive divisors of A, including 1 and A. Map

    m:D(A)×D(B)→D(AB),  m(a,b)=ab.

For d|AB put c(d)=#{(a,b):ab=d}. Writing the p-adic exponents of A,B,d as r,s,t,

    c(d)= product_p [min(r,t)−max(0,t−s)+1].                 (1)

The brackets are positive for d|AB. This follows by independently choosing the exponent i of a in the interval `max(0,t−s)≤i≤min(r,t)`.

**Theorem A (native measure repair).** Among measures constant on every multiplication fibre, the unique one pushing forward to native divisor-value weight d is

    w(a,b)=ab/c(ab).                                        (2)

Indeed the c(d) equal summands must add to d. Thus every native function f satisfies

    Σ_(a|A,b|B) w(a,b) f(ab)=Σ_(d|AB) d f(d).

The unique pair mapping to AB is(A,B), so deleting it gives total mass s(AB). This statement preserves arbitrary native divisor observables, not just sigma.

**Theorem B (independence iff coprime).** The strictly positive matrix w has rank one over Q (equivalently, its normalized coordinate law is independent) iff gcd(A,B)=1.

If coprime, each c(d)=1 and w(a,b)=ab. Conversely take a common prime p. Restrict rows and columns to1,p. The matrix is

    [[1, p/2], [p/2, p²/c(p²)]].

Here c(p²)∈{1,2,3}, depending on the two exponent bounds. Its determinant is `p²(1/c(p²)−1/4)>0`. Hence w cannot have rank one. This is an obstruction to a product probability measure, not to every possible compression algorithm.

**Hostile 6×6, repaired.** After removing(6,6), the flags “neither top, left-only top, right-only top” have weights

    (22,33/2,33/2),  total 55.

The corresponding three-state flag angle is

    yz/((x+y)(x+z))=9/49,

not36/91. The native mass has been repaired, but the unitary-split spectrum has not. Even correct pushforward weights do not permit transplanting the old independence proof.

## 3. Associative source-preserving repair for arbitrary factor trees

Let A_1,...,A_k≥1 be fixed named factors. The source carrier consists of tuples `(d_1,...,d_k)` with d_i|A_i, mapped to their product d. Let C(d) count this fibre. Declare its desired measure to assign weight d/C(d) to each tuple. This is the unique fibre-uniform lift of native divisor weight d.

Choose any binary tree with these leaves, preserving their identities. At a leaf i set C_i(a)=1 for a|A_i and0 otherwise. At every internal node v with children L,R set

    C_v(d)=Σ_(ab=d) C_L(a)C_R(b).                           (3)

For a fixed root value d, choose a child split with probability

    P_v(a,b | d)=C_L(a)C_R(b)/C_v(d),  ab=d.                 (4)

Zero-count states are omitted, so every denominator is positive.

**Theorem C (shape-independent source law).** Recursing with(4) produces the uniform law on all leaf tuples of total d, for every binary association tree. Multiplying by root weight d gives the exact declared source measure. Every joint named leaf guard, not merely every product observable, has the same expectation under every tree.

**Proof.** Equation(3) counts extensions because child leaf sets are disjoint and each pair of child tuples combines uniquely. In the product of all conditional probabilities(4), every internal count cancels once in numerator and denominator. The remaining numerator is the product of leaf counts, all1. The denominator is C_root(d), itself independent of the tree because it counts the same source tuples. Thus every leaf tuple has probability1/C_root(d). This also proves normalization and exact equality for arbitrary functions of all named leaves. ∎

In p-adic coordinates, C_v is the coefficient table of

    product_(i below v) (1+x+...+x^v_p(A_i)),

and tables for different primes multiply. This gives an exact convolution implementation without enumerating all full source tuples. It is not a polynomial-time claim in the binary lengths of arbitrary unfactored integers: factorization, divisor counts and message support sizes still cost resources.

**Minimal hostile.** For A_1=A_2=A_3=p, condition on total p. At the root of `((pp)p)`, the two product splits `(1,p)` and `(p,1)` are not equiprobable: their descendant multiplicities are1 and2. Correct probabilities are1/3 and2/3. Resetting them to1/2, then splitting the left p uniformly, produces `(1/4,1/4,1/2)`. The other bracketing gives `(1/2,1/4,1/4)`. All three laws push forward to the same total divisor p; the first-leaf guard distinguishes them.

Deleting the top root value simply restricts the same proof to proper divisors. Changing the root weighting from d to any prescribed nonnegative μ(d) also leaves the conditional proof intact. Uniform fibre weights are a declaration, not a claim that an arbitrary previously supplied source measure should be overwritten. For nonuniform original sources, carry their extension **masses**, not counts, and use the same product/sum argument when the declared source law factors over leaves.

## 4. Marginal answers are not reusable Boolean guards

A function on sources descends to the native product quotient iff it is constant on every product fibre. For two factors sharing p, the guard “left divisor is p” separates `(p,1)` from`(1,p)` and therefore does not descend. Its conditional mean on total p is1/2. Replacing the Boolean guard by this mean changes its algebra:

    E[g²|p]=1/2,  (E[g|p])²=1/4.                           (5)

More generally the defect `E[fg|d]−E[f|d]E[g|d]` is the conditional covariance. It is zero for every pair of source functions iff every fibre has size1 (under strictly positive fibre weights). For a non-singleton fibre choose its indicator g; its conditional variance is positive.

For k copies of p at total p, the k named leaf guards each have mean1/k, but distinct guards have product identically0, not1/k². Exactly one named leaf owns p. Recovering all leaf guards recovers the source tuple; storing only the common marginal1/k does not. Under the declared uniform law the remaining uncertainty is log2(k) bits, not automatically an integer number of binary guards. This distinguishes information entropy from the number of named tests.

Equation(4) preserves the full law and all these products without pretending the joint law is a product of individual responses. It gives a concrete arithmetic model of why the tournament cursor interface needs its joint order sidecar.

## 5. Connection contract and limits

| Source | Target/map | Preserved | Lost if forgotten | Repair |
|---|---|---|---|---|
| Divisor tuples | Product divisor | Every native sum after(2) | Named source coordinates, coordinate products | Full fibre law or sufficient source-specific data |
| Binary factor tree | Same named leaf tuple space under rebracketing | All joint guard expectations via(3)-(4) | Extension multiplicity if internal nodes are reset | Integer multiplicity/mass messages |
| LCA/shared ports | Arithmetic extension-count messages | Each physical source counted once; compatibility conditioned jointly | Separately uniform child answers | Exact separator product state and count |
| Tournament run-order witness | Conditional analogy only | Mechanism “same local responses, different joint assembly” | No actual tournament-to-divisor map provided | A separate proved map is required before numeric transfer |

No free reconstruction of arbitrary hidden relations is claimed. The arithmetic relation here is explicitly multiplication of independent named divisor choices; this is what makes count messages sufficient. For an additional cross-leaf guard, its dependence must be carried at a scope where it can be tested. In particular replacing count messages by their support alone reintroduces the hostile.

All divisions are over Q or R. Reduction modulo a prime dividing a count denominator is not defined by these formulas; an integer-lift or different carrier is required. No Frobenius cancellation is inferred from a rational conditional law.

## 6. Reproducibility

Run:

    python3 check_clone_associativity.py
    python3 -O check_clone_associativity.py

Both modes return **49,074 active exact checks** over 1,600 factor pairs and 79 association trees, with 9 named factor families. The checker verifies native pushforwards, proper-divisor masses, independence iff coprime, exact source probabilities under every tested tree, both rebracketing hostiles, and the 6×6 angle correction. Proofs above establish general statements; the finite tests are regression controls.

**Independent audit:** sibling audit accepted the proofs of A–C and the guard/scope boundaries, and independently checked 256 four-factor rows over five shapes, totaling 17,925 fibre checks under normal and `-O`. See the session audit packet; these are additional controls, not included in the 49,074 count.
