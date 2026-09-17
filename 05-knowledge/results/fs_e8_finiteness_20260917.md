# Fargues-Scholze fibers for E8: the exact remaining depth obligation

**Date:** 2026-09-17. **Status:** PROVED reduction from CITED inputs;
CONDITIONAL affirmative theorem; unrestricted target remains OPEN in this
research attempt. No claim of a new unconditional solution or priority.

## Target and outcome

Let F be a finite extension of Q_p, let G be a group of type E8,
fix a prime ell different from p, and put k = Qbar_ell. Interpret the
question as finite fibers of the representation-level map

```text
FS_G : Irr_k(G(F)) -> {semisimple Weil L-parameters}/Ghat(k)-conjugacy.
```

The domain consists of isomorphism classes of smooth irreducible
representations. This is not a claim that a map of parameter stacks is a
finite morphism, or that the inertia-only parameter map has finite fibers.

Working with the split form does not restrict the question. For the split
E8 group G_0, the center and diagram automorphism group are trivial, so
Aut(G_0)=G_0. Kottwitz, Proposition 6.4, identifies H^1(F,G_0) with the
character group of pi_0(Z(Ghat_0)^Gamma), which is trivial. Hence every
F-form of E8 is split. This uses the p-adic theorem, not the paper's
differently scoped global statements;
see [Kottwitz, printed p. 627](https://math.uchicago.edu/~drinfeld/langlands/Kottwitz.pdf).

The strongest unconditional conclusion established here is the equivalence

```text
For each fixed phi:
  FS_G^{-1}(phi) is finite
    iff
  {depth(pi) : FS_G(pi) = phi} is bounded above.
```

The proof works for connected reductive groups more generally and does not
require p>7. Thus a counterexample to finite fibers would necessarily contain
representations of arbitrarily large depth with the same FS parameter.

The requested E8 statement also follows from a precise additional
modular-functoriality hypothesis by Cotner-Feng, Theorem 10.2.1. That
hypothesis is still explicit in the version retrieved on this date. The
[primary-source audit](../reference/CORE-PAPERS-FS-E8-FINITENESS-2026-09-17.md)
records versions, theorem numbers, and the distinction between residue
characteristic and auxiliary primes.

## Inheritance, portfolio, and connection contract

No earlier repository theorem about this FS map was located in the bounded
startup search. The closest proved mechanism is DHKM's finiteness of the
excursion action at fixed compact-open level. The canonical hostile is an
infinite union of finite sets: N -> {point}, filtered by {0,...,r}. The
corrected near miss is confusing the residue-prime restriction p>7 with
the auxiliary-prime hypotheses in modular base change. The least-used
relevant sidecar is representation depth.

Anchor: finite E8 fibers. Niche: Hecke-algebra specialization at a fixed
excursion character. Wildcard: whether forgetting monodromy itself could
explain an infinite fiber. The five live concepts are level, depth, wild
inertia, auxiliary primes, and semisimplification. Their interaction is
recorded below rather than inferred from the common name E8 in the earlier
octonion/lattice session.

The operative connection is

```text
source: irreducible representations with nonzero K-invariants
map:    pi |-> pi^K, then specialize the central excursion action at phi
target: simple modules over a finite-dimensional algebra H_{K,phi}
preserved predicate: distinct irreducible representations remain distinct
lost information: representations with no K-invariants are absent
sidecar: one K detecting every representation in the fiber
test:   a depth bound supplies such a K; fixed parameter alone is not yet
        shown to supply it.
```

## 1. Unconditional finiteness at a fixed compact-open level

Fix a compact open pro-p subgroup K of G(F). Choose a square root of the
residue-field cardinality q in k, and let R = Z_ell[sqrt(q)] inside k.
Write

```text
P_R = c-Ind_K^{G(F)} R,
A_R = Exc(W_F,Ghat)_R,
H_R = R[K\G(F)/K].
```

Normalize Haar measure by vol(K)=1. P_R is generated as an R[G(F)]-module
by the characteristic function of K. Its K-invariants identify with H_R
(up to the harmless opposite-algebra convention for endomorphisms).

By DHKM, *Finiteness for Hecke algebras of p-adic groups*, Corollary 3.5(1),
P_R is admissible over A_R through the FS action. In particular H_R is a
finite A_R-module. That action is central. The double-coset basis gives

```text
H_R tensor_R k = H_k = k[K\G(F)/K],
```

so H_k is finite over A_k := A_R tensor_R k. For a semisimple parameter phi,
let chi_phi : A_k -> k be its excursion character. Then

```text
H_{K,phi} = H_k tensor_{A_k,chi_phi} k
```

is a finite-dimensional k-algebra.

If FS_G(pi)=phi and pi^K is nonzero, the simple Hecke module pi^K factors
through H_{K,phi}. The standard idempotent-corner correspondence is injective
on irreducible representations with nonzero K-invariants. A finite-dimensional
algebra has only finitely many isomorphism classes of simple modules, since
its semisimple quotient is a finite product of matrix algebras. Therefore

```text
{pi : FS_G(pi)=phi and pi^K != 0}
```

is finite. No integral lattice in pi or phi, and no commutation assertion
for Bernstein centers under arbitrary coefficient base change, is needed.
We also do not identify A_R tensor_R k with a new inverse-limit excursion
algebra; its displayed scalar extension is all the proof uses.

Injectivity here does not assert a Morita equivalence for the entire
category of K-generated representations. Given a simple H_k-module S,
c-Ind_K^{G(F)}(k) tensor_{H_k} S is generated by its K-invariants S.
Every proper submodule has zero K-invariants. Exact averaging over K
implies their sum also has zero K-invariants, so there is a unique
irreducible quotient. This recovers pi from pi^K.

There is also a uniform quantitative consequence, though no effective
numerical bound is extracted here. If H_R is generated by m_K elements
over A_R, then dim_k H_{K,phi}<=m_K for every phi. Its semisimple quotient
is a product of matrix algebras, so

```text
sum_{FS_G(pi)=phi, pi^K!=0} (dim_k pi^K)^2 <= m_K.
```

In particular the number of such irreducibles is at most m_K, uniformly
in phi. This uses only an injection into the simple modules of the
specialized algebra; surjectivity onto all those modules is unnecessary.

**Dependencies:** the FS central action; DHKM Corollary 3.5(1); the usual
Hecke/idempotent equivalence for irreducibles detected by K; elementary
finite-dimensional algebra. See the
[DHKM primary PDF, printed p. 11](https://arxiv.org/pdf/2203.04929).

## 2. A depth bound supplies a common compact-open level

For every fixed real r>=0, the category of representations of depth at most
r has a finitely generated projective generator over R. This bounded-depth
input is explicitly used immediately before DHKM Corollary 3.5.

Choose finitely many G(F)-module generators for this projective generator
and a compact open pro-p K fixing all of them. Every object in the category
is a quotient of a sum of copies of that generator, hence is generated by
K-fixed vectors. This holds over k as well: regard a k-representation as an
R-representation, take an R-linear generating surjection, and extend it
to a k-linear surjection by p tensor a |-> a f(p). Fixed-vector depth
conditions do not change under restriction of scalars. Thus every
irreducible of depth at most r has nonzero
K-invariants. Section 1 now proves

```text
{pi : FS_G(pi)=phi and depth(pi)<=r} is finite.
```

Choosing m_K for this common K even gives a finite N_r independent of phi
such that the cardinality of that set is at most N_r. The depth bound r
is fixed before choosing K and N_r.

If the entire phi-fiber has bounded depth, choose one such r and apply this
statement. Conversely, a finite set of smooth irreducibles has bounded
depth, because each smooth irreducible has finite depth: a nonzero smooth
vector is fixed by an open subgroup, which contains G_{x,s+} for sufficiently
large s at a fixed building point x. The empty fiber
is included: it is finite and vacuously bounded. This proves both directions
of the displayed equivalence in the target statement.

The conclusion is finite-to-one in the set-theoretic sense; it does not
assert one uniform cardinality bound for all phi.

## 3. The E8 arithmetic and the conditional closing step

Fintzen's *Types for tame p-adic groups*, Table 1, gives

```text
|W(E8)| = 2^14 * 3^5 * 5^2 * 7 = 696729600.
```

For a prime p, p>7 is exactly p not dividing this order. This is an
elementary check on the factors, independently consistent with multiplying
the invariant degrees 2,8,12,14,18,20,24,30. It is not a computational test
of the FS correspondence. Split E8 already satisfies tame splitting.

[Cotner-Feng II, Theorem 10.2.1](https://arxiv.org/html/2609.16387v1#S10.SS2)
proves finite fibers and depth preservation for tame-split groups with
p not dividing the Weyl-group order, assuming its Hypothesis 4.1.4 for
cyclic base change at all auxiliary primes different from p.

In the requested range the implication is therefore rigorous:

```text
cyclic-base-change instances of Hypothesis 4.1.4
    => depth(pi)=depth(FS_G(pi))
    => bounded depth inside each fixed phi-fiber
    => finite fibers.
```

One way to locate the essential mathematical input is the wild-inertia
comparison for Yu representations. Exhaustion by Yu's construction in this
prime range supplies a description of supercuspidals, but it does not by
itself identify their FS parameters. The comparison makes the last depth
in a normalized Yu datum visible in the fixed parameter. Only then is
the formerly unbounded depth coordinate controlled. Compatibility with
parabolic induction and finiteness of parameter lifts from Levi subgroups
extend the result from supercuspidals to all irreducibles.

The three prime roles must stay separate:

| Role | E8 condition | Effect |
|---|---|---|
| Residue characteristic p | p>7 | Prime to the Weyl-group order |
| Fixed coefficient prime ell | ell!=p | Defines the chosen FS map |
| Auxiliary prime used in modular base change | varies over primes !=p | Current comparison retains hypotheses at 2,3,5 |

The retrieved paper establishes the relevant functoriality for auxiliary
primes greater than five and explicitly expects the remaining cases in a
future revision. Choosing a fixed coefficient ell>5 does not remove the
other auxiliary primes from this proof. E8's bad primes are 2,3,5; its Weyl
order additionally contains 7.

## 4. Attempted bypasses and the first failed implications

### Keep only wild inertia

Finiteness needs less than a complete parameter calculation. However, the
available proof of the wild-inertia comparison already uses the missing
base-change compatibility when reducing the tame splitting degree by its
prime divisors. The exact invocation is in
[Cotner-Feng I, Section 10.3.4, equations (10.3.2)-(10.3.3)](https://arxiv.org/html/2609.16381v1#S10.SS3.SSS4).
Omitting the later full-inertia calculation does not remove
that earlier invocation. A tame extension can have degree divisible by
2,3,5 even when p>7. Thus residue characteristic alone does not discharge
the auxiliary assumptions.

The unconditional comparison in Cotner-Feng II, Theorem 1.1.1(1), for
supercuspidals with an attached maximally unramified torus is a valid
restricted positive result. A split ambient group does not imply that all
its relevant tori are maximally unramified. This subclass cannot simply be
declared exhaustive.

### Reverse the bounded-depth continuity statement

DHKM Corollary 3.5(2) controls the wild-inertia cutoff needed for a chosen
representation-depth bound. The required implication runs the other way:
a fixed parameter, or a fixed wild cutoff, would need to control
representation depth. Reversing these quantifiers is not justified.

The hostile filtered map N -> {point} has finite fibers on every bounded
piece and an infinite total fiber. It isolates exactly the missing
coordinate; it is a counterexample to the inference, not to FS finiteness.

### Infer a global packet from finite Bernstein pieces

Finiteness inside each Bernstein component leaves open how many components
meet a fixed parameter fiber. Likewise a product of semisimple Artinian
categories need not have finitely many factors. Compact objects can have
finite support while the entire product remains infinite. The depth
sidecar addresses this missing global control.

### Use a theorem that assumes finite fibers

DHKM, *Local Langlands in families: The banal case*, Theorem 8.2 is not an
independent closure: Definition 6.19(C1) already assumes finite fibers for
the supplied correspondence. The repair is to verify that assumption
first, not to feed the desired conclusion back in as a premise.

## 5. A second route: semisimplification does not itself cause infinitude

There is a useful independent check on the target's forgotten information.
For a fixed semisimple Weil parameter lambda, compatible monodromy is
classified by orbits of H_lambda = Z_Ghat(lambda) on

```text
V_lambda = {N in Lie(Ghat)^{lambda(I_F)} :
            Ad(lambda(Fr))N = |Fr| N}.
```

Use a convention for Fr consistently; |Fr| is q or q^{-1} and is not a
root of unity. This is the Weil-Deligne convention
Ad(lambda(w))N=|w|N. Vogan, *The local Langlands conjecture*, Proposition
4.5, Corollary 4.6, and Lemma 4.8 imply that there are finitely many such
orbits; see the [primary text, pp. 32-34](https://math.mit.edu/~dav/md.pdf).

Here is an added algebraic proof sketch of the mechanism (Vogan states the
key lemma without proof). Put L=Z_Ghat(lambda(I_F))^0, g=Lie(L),
theta=Ad(lambda(Fr)), and t=|Fr|. The central torus has only root-of-unity
theta-eigenvalues, so g_t has zero central component. For N in g_t,
ad(N) sends a theta-eigenspace with eigenvalue a into the eigenspace with
eigenvalue ta. As there are only finitely many eigenvalues and t is not a
root of unity, ad(N) is nilpotent. Hence N is nilpotent in the semisimple
part. There are finitely many nilpotent
L-orbits. For one meeting g_t, its intersection with g_t is the fixed
locus of the diagonalizable Zariski closure of t^{-1}theta, hence smooth
in characteristic zero. Its tangent space at N is
[g,N] intersect g_t = [g_1,N], exactly the centralizer-orbit tangent
space. These orbits are open in the intersection, so quasi-compactness
gives finitely many. The argument is algebraic over any algebraically
closed characteristic-zero field, so Vogan's complex formulation applies
to k here as well. Finite components of the centralizers do not spoil
finiteness. An arbitrary smaller conjugating group would not supply this
eigenspace/tangent-space identity.

Consequently a compatible classical LLC with finite packets would imply
finite FS fibers. But the required full E8 correspondence and its
compatibility cannot be introduced as already-proved inputs here. This
route identifies a sufficient theorem, not a new proof of it.

## 6. The precise weaker lemma worth pursuing

Exact depth preservation is more than this question requires. It suffices
to establish the following statement, with no effective formula for B_phi:

```text
For every semisimple parameter phi, there exists B_phi<infinity such that
FS_G(pi)=phi implies depth(pi)<=B_phi.
```

A stronger but still one-sided target is a function B_{F,G}(r) bounding
representation depth whenever parameter depth is at most r. An estimate
depth(pi)<=C_{F,G}(1+depth(FS_G(pi))) would be enough. These are OPEN
targets here, not consequences of the fixed-level Hecke argument.

This gives a concrete stopping reason for the attempt: the unconditional
inputs control each level and the conditional comparison controls the
missing depth coordinate; no unconditional step connecting those two
controls was established. No counterexample to the requested statement
was found, and no missing hypothesis is being silently promoted to a
theorem.

A more local proposed replacement for full small-prime functoriality is
boundedness across the explicit tame toral lifts: for each relevant cyclic
extension E/F of degree 2,3,5, bound the FS depth of pi(Psi_E) in terms of
the FS depth of pi(Psi), uniformly over the Yu data under consideration.
In the existing induction, bounded-length tame splitting towers and
descent through proper twisted Levis lead to the unramified endpoint.
Such one-sided estimates, together with the established comparisons at
the other steps, would retain the needed depth control without requiring
equality of parameters. This is a research proposal, not an established
bound. Valuations must be normalized consistently: the tame lifts in
Part I use the extension of the original field's valuation.

## Audit and reproducibility boundary

The result is a proof reduction and a dated primary-literature audit, not
a numerical packet enumeration. The load-bearing proof was checked
independently for coefficient extension, excursion specialization,
injectivity on irreducibles, and the bounded-depth generator. No theorem
identifier was reserved or promoted in canon. The unrestricted E8 claim
retains OPEN status in this note.
