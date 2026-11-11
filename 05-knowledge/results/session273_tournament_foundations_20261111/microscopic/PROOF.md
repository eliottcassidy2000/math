# A whole entropy spectrum above one macroscopic tournament

**PROVED in this packet; no Chan–Pak closure.** New contribution: a quantitative
sampling obstruction and an arbitrary-entropy lifting theorem, not the inherited
ordered-join product or the deletion transform.

## 1. Source and target

A tournament has one directed arc between every distinct pair. H counts directed
Hamiltonian paths. Ordered join puts all cross-arcs from earlier to later blocks;
its inherited law is H(A→B)=H(A)H(B). A singleton has H=1.

The macroscopic observer is normalized arc-edit distance, hence also the weaker
cut-distance and all fixed-size induced sampling distributions. The microscopic
observer is the ordered strong-component list. We use no asserted reduction of
Yang–Mills or any continuum conjecture to tournaments.

## 2. Arbitrary entropy above the transitive limit

**Theorem 1.** For each real λ≥0 there is a sequence T_n of n-vertex tournaments
such that

1. T_n differs from a transitive tournament in O_λ(n) arcs;
2. log H(T_n)/n tends to λ;
3. every fixed-size unlabelled induced-subtournament distribution tends to the
   point mass at the transitive tournament.

Thus even the exponential rate of H is not determined by the limiting
tournament kernel. This statement is about existence for each fixed λ; the
constant implicit in O_λ(n) need not be uniform in λ.

**Proof.** For a uniformly random labelled s-tournament, every permutation is a
Hamiltonian path with probability 2^{-(s-1)}, so some s-tournament S has
h=H(S)≥s!/2^{s-1}. In particular log(h)/s is unbounded as s grows. For λ>0 fix
one S with log(h)/s>λ. Let r_n=floor(λn/log h). For large n, sr_n≤n; let T_n be
r_n ordered copies of S followed by n−sr_n singletons. Then
H(T_n)=h^{r_n}, and log H(T_n)/n→λ. Fix a linear order inside S. Only internal
block arcs disagree with the global transitive order, at most
r_n binom(s,2)=O_λ(n). A uniform k-set sees a discrepant arc with probability
at most O_λ(n) (k)_2/(n)_2→0, proving the sampling assertion. The λ=0 case is
TT_n. ∎

This is not merely a fixed exceptional core hidden by padding: a positive
fraction of vertices may belong to nontrivial components. For the explicit
triangle chain T_{3r}=C3→...→C3, H=3^r, while just r arc reversals make it
transitive. Its entropy is log(3)/3, not zero.

A growing-block variant gives log H(T_n)/n→∞ with o(n²) discrepant arcs: choose
s_n→∞ with s_n=o(n), use floor(n/s_n) copies of a maximizing s_n-tournament,
and fill the remaining vertices with singletons. The lower bound above gives
entropy at least (1−s_n/n)(log(s_n!)/s_n−(1−1/s_n)log2), which diverges.

## 3. A sharp cubic sampling scale for a fixed core

Let T_n=TT_{n−s}→S, where S is fixed and nontransitive. Let c_3(S) be the number
of directed cyclic triples. Let μ_{n,k} be the unlabelled law of a uniform
k-vertex induced subtournament, for fixed k≥3. Then

    ||μ_{n,k}−δ_{TT_k}||_TV
      = Pr(the sample contains a cyclic triple)
      = c_3(S) (k)_3/(n)_3 + O_{S,k}(n^{-4}).             (1)

The upper bound without an error term is c_3(S)(k)_3/(n)_3.

**Proof.** A tournament is transitive iff it has no directed triangle. Every
triangle of T_n lies in S. Each particular cyclic triple is included with
probability (k)_3/(n)_3. The union bound gives the upper bound. Distinct triples
have a union of at least four vertices, so Bonferroni's pair error is
O_{S,k}(n^{-4}). The identity with TV follows since the comparison law is a
point mass. ∎

For q independent samples the TV distance from q samples of TT_n is at most
q times (1). Hence any test restricted to independent unlabelled k-samples
needs q=Ω_{S,k}(n³) for constant distinguishing advantage. Access to marked
core vertices, adaptive neighbourhood queries, or a full adjacency matrix is
NOT included. This is an information-access theorem, not a runtime lower bound
for arbitrary tournament algorithms.

For r disjoint ordered triangle blocks on n vertices the exact sampling law is

    Pr(transitive sample)
     = [z^k](1+3z+3z²)^r(1+z)^{n−3r} / binom(n,k).       (2)

At n=3r the nontransitive probability is O_k(n^{-2}), whereas H is exponential.
This strengthens the crude discrepant-edge bound of Theorem 1.

## 4. What restores the lost information?

Let ν_T(S)=(number of strong components of T isomorphic to S)/|T|, including
singletons. Then the exact microscopic identity is

    log H(T)/|T| = sum_S ν_T(S) log H(S).                 (3)

This is an ordinary finite sum; no limiting exchange is implicit. Component
order is irrelevant to this scalar but indispensable for the graph class.
Convergence of ν_T alone does not pass (3) to the limit without integrability:
unbounded component order permits escape of log H-mass. For a fixed finite
component alphabet, convergence does suffice.

A higher-resolution motif version is also exact. If F is strong of order k≥3,
all induced copies of F in an ordered join lie in one component. Thus for a
chain of r copies of an s-vertex S, with n=rs,

    n^{k−1} Pr(T_n[k]≅F) → (k!/s) ind_F(S),              (4)

where ind_F counts unlabelled vertex subsets inducing F. Ordinary limiting
motif densities are zero; their renormalized coefficients retain the hidden
component information. Formula (4) by itself does not reconstruct arbitrary
varying cores.

## 5. Connection and research consequence

- **Tournaments → infinite-limit methods:** convergence of macroscopic local
  observables can coexist with every prescribed finite H-entropy and even
  divergent entropy. A limiting classification needs a microscopic component
  measure and a uniform-integrability sidecar before claiming exact arithmetic.
- **Continuum/contact-data idea → tournaments:** the Session272 observation that
  separated correlators can forget contact terms suggests retaining the scaled
  motif coefficients (4), not just the limiting kernel. This is a structural
  analogy implemented by (1)–(4), not a theorem about Yang–Mills fields.
- **Chan–Pak:** neutral padding already hides each attained exact h in the same
  transitive limit. The theorem shows the loss persists at the extensive scale.
  Proving the right graph limit or a smooth entropy law cannot certify the
  missing integer targets. A constructive source interface must retain exact
  finite components and endpoint data.
- **A000568:** an ordered word over any fixed collection of q distinct strong
  s-vertex classes gives q^r distinct ordinary classes on sr vertices, all
  approaching the same transitive limit. Their H values may coincide or differ;
  the class entropy log(q)/s and Hamiltonian entropy are distinct observers.
  This is a lower subfamily, not the superexponential full A000568 census.

## Provenance

The ordered-join product is standard and inherited throughout the repository,
including THM-1370 (`h-spectrum-omits-7-21-all-n`). Session268's interacting
CA family supplies an earlier concrete example of exponentially many classes
with far fewer H-values; we generalize the macroscopic loss, not its cell
construction. THM-3372 (`multiaffine-deletion-transform-variance-and-skew-join-current`)
already retains full deletion responses; (4) is a scale-separated reading, not
a new deletion polynomial. The old GMC factorial-Hankel faithfulness theorem
THM-2810 inspired checking which unbounded information survives a finite/limit
carrier; its proof is not used here. All claims above have elementary proofs
independent of the repo's external major-conjecture claims.
