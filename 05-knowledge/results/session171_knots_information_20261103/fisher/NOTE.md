# Holonomy information: a probe bank preserves a subspace of questions

**Status: PROVED finite information identities and explicit applications; FINITE-EXACT controls.** The Shannon, Fourier, and Fisher ingredients are standard finite-dimensional mathematics. The contribution is their precise composition with S170's signed-unit holonomy codes and probe addresses, together with exact examples separating identity certification, partial evidence, and noisy readout. No literature-priority claim is made.

The key change of question is from “how many sheets were tested?” to **“which questions about the return remain answerable?”** A bank can preserve one question perfectly and erase another completely, even when both are single-bit questions. Probe rank determines total information under a uniform prior; a dual subspace determines its meaning.

## 1. Inherited carrier and a new observer vocabulary

[S170 holonomy tomography](https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/05-knowledge/results/session170_holonomy_homotopy_20261103/holonomy_code/NOTE.md) supplies the following explicit interface. At Cayley–Dickson level n≥3, write d=n−1. A return promised to lie in the rth lower-central subgroup, 2≤r≤d, acts on signed units by

    D_h(s,y,t)=(s+h(y),y,t),
    h∈C=RM(d−r,d),
    K=dim C=Σ_(j=0)^(d−r) binom(d,j).

Here C is the binary vector space of truth tables of square-free Boolean polynomials of degree at most d−r. The code's minimum distance is 2^r. An exact sign readout at addressed initial sheet y measures the linear functional h↦h(y). Initial sign and the newest bit t do not change that functional. Balanced fixed-label closed binary words satisfy the r=2 promise; arbitrary transport words need not.

For a list of sheet addresses S=(y₁,…,y_q), including possible repetitions, write

    E_S:C→F₂^q,       E_S(h)=(h(y₁),…,h(y_q)),
    a=rank E_S,       N=ker E_S.

Call N the **blind defect space**. Define the **answerable character space**

    A_S = {ℓ∈C*: ℓ vanishes on N} = im(E_S*).

These are operational definitions, not new names for a new code family. The first space describes indistinguishable defects; the second describes exactly which linear binary questions about them remain answerable.

The uniform law used below is a declared experiment on C. It is not asserted to be the law of naturally occurring runner movies, random tournaments, or uncontrolled circuits. S170 realizes all r=2 codewords by suitable scheduled binary circuits; that construction does not give their occurrence probabilities in a separate application.

## 2. Exact Shannon budget, with a distribution-free boundary

All entropies use logarithms base two. Let X be a uniform random element of C and let T=E_S(X).

**Theorem 1 (linear observation budget).**

    H(X)=K,    H(T)=I(X;T)=a,    H(X|T)=K−a.

Every feasible transcript has exactly 2^(K−a) compatible defects, uniformly distributed conditionally. A function V:C→R is exactly recoverable from T if and only if it is constant on every coset of N. In particular, a linear question ℓ is recoverable precisely when ℓ∈A_S.

**Proof.** Rank–nullity gives |N|=2^(K−a) and |im E_S|=2^a. Every nonempty fiber is a translate of N. Uniformity makes those fibers equiprobable and uniform internally, giving all entropy identities. A function factors through E_S exactly when it is constant on its fibers. For linear ℓ, this is precisely ℓ|_N=0. ∎

For an arbitrary law p on C, the fibers are unchanged but their probabilities need not be uniform. The general statements are instead

    H_p(X|T)≤K−a,        I_p(X;T)=H_p(T)≤a.

The first bound is attained exactly when every positive-probability conditional fiber law is uniform on the entire fiber. The second requires uniformity on the whole image. If p has full support, exact task recovery still requires constancy on all fibers; otherwise only its support matters.

Thus “K−a bits are lost” is an exact uniform-prior assertion, not an unconditional count of uncertainty for every source law. The worst-case indistinguishability fiber size is distribution-free.

## 3. Fourier decomposition gives an exact Fisher loss ledger

For ℓ∈C*, set χ_ℓ(x)=(-1)^ℓ(x). Under the uniform law the 2^K characters form an orthonormal basis of real functions on C. Write

    f(x)=Σ_(ℓ∈C*) f̂(ℓ)χ_ℓ(x),
    f̂(ℓ)=2^(−K)Σ_(x∈C) f(x)χ_ℓ(x).

**Theorem 2 (answerable Fourier projection).** Conditional expectation onto the probe transcript deletes exactly the characters outside A_S:

    E[f(X)|T] = Σ_(ℓ∈A_S) f̂(ℓ)χ_ℓ(X).

For the finite exponential-tilt experiment

    p_θ(x)=exp(θ f(x))/Σ_z exp(θ f(z)),

the full and observed Fisher informations at θ=0, with natural logarithm in the score, are

    I_X(0)=Σ_(ℓ≠0) f̂(ℓ)²,
    I_T(0)=Σ_(ℓ∈A_S, ℓ≠0) f̂(ℓ)²,
    I_X(0)−I_T(0)=Σ_(ℓ∉A_S) f̂(ℓ)².

**Proof.** Conditional on T=E_S(x), X is uniform on x+N. Its character average is

    χ_ℓ(x)|N|^(−1)Σ_(z∈N)χ_ℓ(z).

If ℓ vanishes on N, this is χ_ℓ(x). Otherwise choose z₀∈N with ℓ(z₀)=1; pairing z with z+z₀ makes the average zero. Linearity proves the projection formula. Differentiation of the finite positive exponential family gives score f−E f. The observed score is its conditional expectation, by differentiating each fiber sum. Orthogonality then gives the three sums. ∎

This explicitly instantiates the earlier repository identity

    I_X−I_T=E Var(score|T)

from [Information-preserving observers still need their source and action markers](https://github.com/eliottcassidy2000/math/blob/2d137c00eaf36a6fbe9f06044ac56d2c3bf7c7d6/05-knowledge/results/information_observers_collatz_20261014a.md), §2. We do not claim that identity as new.

**Corollary 2.1 (one-bit questions can be wholly visible or wholly invisible).** For a nonzero ℓ, the tilt by χ_ℓ has I_X(0)=1 and observed information either 1, if ℓ∈A_S, or 0 otherwise. Equivalently the Boolean tilt f=ℓ has full information 1/4 and observed information either 1/4 or zero.

Under the uniform source, a nonzero ℓ outside A_S is a fair bit independent of T. If ℓ∈A_S, it is a deterministic function of T. This is exact, not an asymptotic assertion about small noise.

For several real statistics f₁,…,f_b, the same calculation gives the matrix identity

    J_full−J_observed = Σ_(ℓ∉A_S) v_ℓ v_ℓᵀ,
    v_ℓ=(f̂₁(ℓ),…,f̂_b(ℓ)).

Hence the missing Fisher matrix is positive semidefinite and its null directions are precisely linear combinations of the supplied statistics whose nonconstant Fourier support lies in A_S. This is local information at the uniform parameter, not a claim about every nonlinear parametric family or every parameter value.

## 4. Identity testing: partial evidence is not a certificate

Set V(x)=1 if x=0, and zero otherwise. Let ε=2^(−K). Its Boolean exponential tilt asks whether observations distinguish increased mass on the identity defect.

**Corollary 3 (exact identity information).**

    I_full=ε(1−ε),
    I_observed=2^(a−2K)−2^(−2K),
    I_lost=2^(−K)−2^(a−2K).

If a<K, the all-zero readout still leaves 2^(K−a) possible defects and cannot certify identity. Nevertheless, for a>0, that readout has positive statistical information about the identity tilt.

**Proof.** The conditional probability of identity is 2^(−(K−a)) on T=0 and zero otherwise. The zero transcript has probability 2^(−a). Taking the variance of this conditional probability gives I_observed; Var(V)=ε(1−ε) gives I_full. ∎

The all-zero transcript is compatible with every x∈N, including nonzero defects when a<K. This is the precise distinction between probabilistic evidence under a supplied prior and an all-sheet identity certificate.

## 5. Concrete S170 examples

### Level four: one sheet retains one bit and loses three

Here n=4, d=3, r=2, so h(y)=c₀+c₁y₁+c₂y₂+c₃y₃ and K=4. Testing y=0 reads c₀. Its rank is one and its blind space has dimension three. Thus it loses exactly three bits under the uniform law.

The question h(0) is completely observed. The question h(e₁)=c₀+c₁ is independent of the observed bit and loses its full Fisher information 1/4 under its Boolean tilt. The identity tilt retains Fisher information 1/256, compared with full information 15/256, a retained fraction of 1/15; the lost amount is 14/256=7/128.

Testing 0,e₁,e₂,e₃ has rank four and recovers every coefficient. Repeating the test at 0 three times still has rank one in the exact noiseless model. Repetition can help under a separate independent-noise model, but it creates no additional noiseless source information.

### Level five: first-order probes miss exactly six quadratic coefficients

Here n=5, d=4, r=2 and K=1+4+6=11. Probe the five points of weight at most one. Their evaluations recover c₀ and the four linear coefficients. All six quadratic coefficients remain free, so the residual entropy is exactly six bits and every transcript has 64 compatible defects.

In particular y↦y₁y₂ is a nonzero balanced defect missed by this entire bank. Evaluating h at e₁+e₂ supplies a question not determined by the five observations. Its Boolean tilt loses all 1/4 Fisher information.

The complete weight-at-most-two bank has eleven points and rank eleven, by Boolean Möbius inversion, exactly the S170 optimal bank. The full sixteen-point cube has the same eleven bits of noiseless information. Its five redundant readout bits are not wasted: the observed full truth table has minimum distance four and hence corrects one arbitrary readout error, whereas an eleven-point information set has image F₂^11 and no uniform one-error protection. Information about the source and protection against corruption are separate resources.

### Adaptive addresses have a transcript-dependent rank budget

Suppose a deterministic procedure chooses each next linear probe using previous replies and eventually stops, with all branching rules known. For any realized complete transcript t, the compatible sources form the affine solution set of its queried equations. All sources satisfying those equations follow that same branch, so there are no additional hidden restrictions. If their rank is a(t), then under the uniform law

    H(X|T=t)=K−a(t),
    H(X|T)=K−E a(T),
    I(X;T)=E a(T).

Thus adaptivity can change expected rank or cost, but it cannot reduce the worst-case all-zero identity-certification requirement below K independent binary probes. On that branch, fewer than K equations leave a nonzero common kernel. This recovers the S170 adaptive lower bound with its information content made explicit.

Unlike a fixed linear bank, an adaptive transcript need not have one global answerable character subspace. The Fourier projection theorem above is stated for the fixed map E_S, not silently extended to this case. The checker includes a three-bit adaptive tree with leaf sizes 4,2,2, expected rank 3/2, and residual entropy 3/2 bits.

## 6. Physical connection, scope, and a bandwidth warning

S170 supplies scheduled constant-profile tournament implementations of the signed-unit probes. The present theorem applies once their returned bit is authenticated as the evaluation h(y) and the subgroup promise is established. It does not establish either premise from a Hamilton count or an unmarked tournament.

The companion S171 contextual coding construction uses a different physical interface: a chosen context can reveal h(y) through an exact Hamilton-count formula. An individually addressed binary response can therefore instantiate this theorem, once that formula is proved. But a context touching several blocks may return a count encoding a subset weight or still more information. Such a measurement is not a single binary linear functional. The K-probe lower bound applies to addressed binary sign/local-bit probes, not to arbitrary real-, integer-, or vector-valued contexts.

The source and map ledger is:

- Source: a signed-unit defect h with an authenticated code-family promise.
- Map: a specified list of exact linear evaluations, with initial sheet addresses retained.
- Preserved: the quotient C/N and exactly the functions constant on its fibers.
- Destroyed: the residual defect inside a fiber; also path words, schedules, geometric winding, and knot embeddings, which were never encoded by h alone.
- Distributional assumption: uniform on C for the exact bit budget and Fourier Fisher formulas; the fiber criterion itself does not need uniformity.
- Decisive hostile: h(y)=y₁y₂ is invisible at every weight-zero or weight-one address, yet changes at e₁+e₂.
- Boundary: source entropy is not a proof certificate, Fisher sufficiency is relative to a specified statistical experiment, and a noisy sensor requires a separate noise model.

No torus-knot type, global runner theorem, or tournament classification follows from a sufficient observation of this finite defect code.

## 7. Reproduction

Run with Python's standard library:

    python3 check.py > results.json
    python3 -O check.py > results.optimized.json
    cmp results.json results.optimized.json

The checker uses exact integers and fractions and explicit exceptions, not removable assertions. It enumerates the complete 16-word and 2,048-word codes, verifies every fiber size, ranks, all nonzero characters in the smaller examples and the declared character subset in the larger examples, the Fourier Fisher identity for explicit multi-character scores, identity-tilt fractions, the two minimum distances, and the adaptive tree. Shannon entropies here are exact integer or rational expectations of log₂ power-of-two fiber sizes, so no floating-point entropy approximation is used.

All-family statements follow from the displayed proofs. The finite certificates verify the examples and implementation, not an extrapolation to larger n.
