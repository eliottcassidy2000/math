# Completing observables through products and transitions

2026-10-20. **PROVED elementary finite statements; FINITE-EXACT controls on eight specified tiny systems.** This separate addendum leaves the earlier simplex/sphere/cube packet unchanged. It connects its consumer-minimal linear completion theorem to deterministic dynamics and binary registers. No historical novelty is claimed for observable quotients, Koopman pullbacks, finite function algebras, or linear-system observability. There is no Collatz convergence or planar Jacobian consequence.

## 1. The exact question

A list of observables can do three different jobs:

1. distinguish deterministic states;
2. evolve its own expected values without retaining all states;
3. identify an arbitrary probability law on the distinguishable states.

These jobs need not require the same number of coordinates. The distinction is already sharp for a three-bit cyclic register: observing one bit for three times identifies the initial state, but the three separate output means do not determine an arbitrary law on the eight initial states. Products retain the missing joint information.

### Inheritance and source boundary

- `05-knowledge/results/observable_hankel_clocks_20261020/OBSERVABLE_HANKEL_CLOCKS.md`, §2, proves a shift-stable observable quotient and stresses that its module size is not the period of one sequence. This addendum changes the requested carrier from a linear span to a pointwise algebra.
- **THM-848**, `01-canon/theorems/THM-848-hamiltonian-path-walsh-needle-stalk.md`, §3 and its operation-boundary discussion, already distinguishes a minimal Krylov observable from a state operation-congruence.
- **THM-3088**, `01-canon/theorems/THM-3088-punctured-projective-direction-algebra-and-exceptional-parity-saturation.md`, identifies a finite quotient by its pointwise idempotent algebra. **THM-3072**, `01-canon/theorems/THM-3072-a4-flag-three-c2-tomography-and-edge-cycle-cospan.md`, warns that averaged tables do not automatically preserve the multiplication of a common coordinate algebra.
- `05-knowledge/results/binary_lfsr_tournament_clocks_20261020/BINARY_LFSR_PRIMARY_CLOCKS.md` supplies the linear-feedback context; the new statement does not replace its primary-exponent sidecar by a cycle graph.
- The earlier session packet `SIMPLEX_SPHERE_CUBE_COMPLETION.md`, Theorems 1 and 3, supplies the consumer-rank rule and the exact full-character positivity test.

For external context, [Anantharaman and Sule, *Koopman operator approach for computing structure of solutions and observability of non-linear finite state systems*, arXiv:2010.10752v1, §§1.1–1.2 and 2](https://arxiv.org/pdf/2010.10752), defines finite-state pullbacks and studies reduced invariant function spaces and observability. Its function spaces are finite-field-valued; the probability algebra below is real-valued. We use the paper for context, not as a dependency for the elementary proofs below.

The contribution here is the explicit bridge between a transition-stable finite quotient, its full simplex of laws, and the minimal character products needed to recover that simplex from linear output statistics.

## 2. General finite-state theorem

Let Ω be a nonempty finite set, let T₁,…,T_k:Ω→Ω be specified total deterministic maps, and let O⊆R^Ω be a finite list of observed functions containing the constant 1. The empty transition word is allowed. For a word w, write T_w for its composed map.

Define the linear pullback closure and its generated pointwise algebra by

    W = span_R {o∘T_w : o∈O, w any finite word},
    A = alg_R {o∘T_w : o∈O, w any finite word}.             (1)

Here alg means a unital real subalgebra under pointwise multiplication. Define

    x∼y iff o(T_w x)=o(T_w y)
            for every o∈O and every finite word w.          (2)

**Theorem 1 (observable-algebra completion).** A is the smallest unital pointwise algebra containing O and stable under every pullback T_j^*:f↦f∘T_j. Moreover,

    A = {f∈R^Ω : f is constant on every ∼-class}.            (3)

Thus dim_R A equals the number q of observational equivalence classes. Every T_j descends to the quotient Q=Ω/∼. This quotient is the coarsest deterministic quotient on which all original observations and chosen transitions are well-defined.

**Proof.** Pullback preserves constants, sums, scalar multiplication and products. Applying a pullback to a word-observable gives another word-observable. Hence A is stable. Every stable algebra containing O contains all word-observables and their products, proving minimality.

Every generator in (1) is constant on the classes of (2), so A lies in the right-hand side of (3). For the reverse inclusion, fix a class C and representative x_C. For every other class D choose a word-observable g_CD with g_CD(x_C)≠g_CD(x_D), which exists by (2). Then

    e_C = product over D≠C of
          (g_CD−g_CD(x_D))/(g_CD(x_C)−g_CD(x_D))           (4)

belongs to A. All factors are class-constant. The product is one on C and zero on every other class, so it is the class indicator. If q=1 the empty product is 1. These indicators span all class-constant functions, proving (3).

If x∼y, applying an arbitrary future word after T_j gives the same observations at T_jx and T_jy, so T_jx∼T_jy. Thus quotient transitions are well-defined. Finally, any quotient retaining observations and transitions must identify only states whose observations agree along every word. Its fibers refine ∼, and it factors onto Q. ∎

This is a statement about the entire declared Ω. A separate restriction to states reachable from a distinguished initial state changes the problem and must be stated before minimization.

### A finite stopping rule

Begin with the partition by the original observed value tuples. Repeatedly split classes according to the tuple of successor classes under T₁,…,T_k. The partition after h refinements records observations along every word of length at most h. A round that makes no split is stable forever. Every strict round increases the number of classes, so at most |Ω|−1 strict rounds suffice. Formula (4), or rational linear algebra when the inputs are rational, then produces the exact quotient indicators.

The stopping rule is finite because Ω is finite. It supplies no finite bound for an infinite arithmetic state space merely by applying a similar notation.

## 3. What the algebra determines about probability

For a law p on Ω, knowing E_p f for every f∈A is equivalent to knowing the entire pushforward law λ on Q: its atom masses are E_p e_C. Information about how mass is split inside one class is unobservable and is not recovered.

The q quotient atoms form a (q−1)-simplex. Its exact square-root image is a spherical orthant in R^q. Quotient dynamics pushes this simplex forward by the deterministic transition matrices, or, dually, pulls back functions in A. This is a genuine cube/character → law/simplex → spherical-amplitude connection with a declared operation.

Knowing only expectations of a basis of W can be much weaker. W is sufficient to evolve its own means, because it is pullback-stable, but may omit products required to recover λ. The minimal-observer theorem from the earlier packet gives the exact number of additional arbitrary scalar expectations needed to identify every quotient law after all W-expectations are known:

    q−dim_R W.                                             (5)

Indeed W⊆A, 1∈W, and A has dimension q. Choosing any basis of A modulo W attains (5); a smaller family cannot distinguish all laws by linear expectations. Positive small perturbations of the uniform law provide the necessity witness.

Importantly, W and A induce the **same deterministic state partition**. A adds no new distinction between individual states already identifiable from their full output records. It instead linearizes all functions of that record, including joint statistics, so that expectations identify arbitrary distributions over the states.

The extra product expectations must actually be available from the same initial random state, or be supplied by a proved model. Multiplying separate means does not provide them. For multiple controllable transitions, products of observations under different words refer to their common initial state; separate experiments with unrelated initial samples do not automatically measure such joint quantities.

## 4. Binary linear and affine dynamics

Now let Ω=F₂^n, including n=0 if desired, and let the observed nonconstant functions be **individually available Walsh characters**

    χ_a(x)=(-1)^(a·x),       a∈S⊆F₂^n.                    (6)

Constants are also included. Consider affine transitions

    T_j(x)=M_jx+b_j,       M_j over F₂.                     (7)

No M_j is assumed invertible. Let L be the smallest F₂-linear subspace containing S and stable under every M_jᵀ. Write r=dim_F₂ L and

    L⊥={v∈F₂^n : a·v=0 for every a∈L}.                    (8)

**Theorem 2 (character closure and quotient).** In this setting,

    A = span_R {χ_a : a∈L},
    dim_R A = 2^r,
    x∼y iff x−y∈L⊥.                                      (9)

There are exactly 2^r distinguishable classes, each of size 2^(n−r). The affine shifts b_j do not change L or the equivalence relation. They do affect the induced quotient transitions and signs of moment evolution.

**Proof.** The elementary identities are

    χ_aχ_c=χ_(a+c),
    χ_a∘T_j=(-1)^(a·b_j)χ_(M_jᵀa).                       (10)

Products add labels, and pullbacks act by transpose matrices up to a nonzero scalar sign. The algebra therefore contains every character labelled by the smallest additive, transpose-stable set containing S, which is L. Conversely the span indexed by L is a unital, multiplication-closed, pullback-stable algebra, proving equality.

Distinct Walsh characters are linearly independent over R, so its dimension is |L|=2^r. Equality of all its character values at x,y is exactly a·(x−y)=0 for every a∈L. The annihilator has dimension n−r, giving the class counts. Shifts in (10) contribute only a common sign to the compared observations and do not change their separating power. ∎

The induced quotient is the affine finite system on F₂^n/L⊥. It is well-defined because M_j(L⊥)⊆L⊥: for v∈L⊥ and a∈L, a·M_jv=(M_jᵀa)·v=0. Translation by b_j then acts on quotient classes normally.

### Compute the labels without materializing the whole algebra

Starting with L₀=span_F₂ S, iterate

    L_(h+1)=L_h + Σ_j M_jᵀL_h.                            (11)

The first equality L_(h+1)=L_h proves stability. Each strict step raises rank, so no more than n−dim L₀ strict steps occur. Gaussian elimination over F₂ computes the r-dimensional label space using n-bit vectors. Materializing all its 2^r characters or every quotient law can still be exponential in r.

For one matrix M, L is the ordinary finite-field observability/Krylov label span of (Mᵀ)^tS. Cayley–Hamilton also shows that exponents 0,…,n−1 suffice when n≥1. This label rank is not the real dimension of the full probability algebra.

### Full probability moments and their evolution

Let μ_a=Eχ_a(X), a∈L, with μ₀=1. The quotient law is recovered exactly by

    λ_[x] = 2^(−r) Σ_(a∈L) μ_aχ_a(x).                    (12)

The right-hand side is constant on each coset. The matrix indexed by a,c∈L,

    H(μ)_(a,c)=μ_(a+c),                                   (13)

is PSD if and only if the λ_[x] are nonnegative. This is the earlier full-character test on the smaller group F₂^n/L⊥; its eigenvalues are 2^r times the quotient masses.

Under transition j the exact moment update is

    μ'_a=(-1)^(a·b_j) μ_(M_jᵀa).                          (14)

Thus (9), (12), (13), and (14) supply state quotient, probability reconstruction, realizability, and dynamics in one compatible carrier. They do not recover a law inside unobserved cosets.

## 5. Three dimensions that must not be conflated

Let K contain 0 and every label obtained from S by any finite word of transpose maps, **without** additive closure. Then

    W=span_R{χ_a:a∈K},       dim_R W=|K|,
    L=span_F₂ K,             dim_F₂ L=r,
    A=span_R{χ_a:a∈L},       dim_R A=2^r.                  (15)

These are respectively the real linear closure of output statistics, the binary state-observation rank, and the real full-law algebra dimension. The minimum missing scalar expectations after W is known is exactly 2^r−|K|. When added observations must themselves be individual Walsh moments, the missing labels L\K supply that minimal family.

### A three-bit LFSR hostile

Take the cyclic register

    T(x₁,x₂,x₃)=(x₂,x₃,x₁),       observed output χ_(1,0,0).

This is the binary linear-feedback register with characteristic polynomial t³−1=t³+1 over F₂. Its first three outputs identify the three initial bits exactly. Here

    K={0,e₁,e₂,e₃},   r=3,   dim_R W=4,   dim_R A=8.

Let p_+ and p_− be uniform on the four words with even and odd total bit parity. They give the same expected output at every time, and hence the same expectation for every element of W. But

    E_(p_+) χ_(1,1,1)=+1,
    E_(p_−) χ_(1,1,1)=−1.

The product of the first three sign outputs distinguishes these laws. A complete law needs the three pair characters and the triple character in addition to W. There is no contradiction with state observability: one full three-output record identifies its initial state, while three marginal means discard the correlations between those outputs.

### Affine shifts preserve the quotient but change its action

On F₂³, observe χ_e₁ and use

    T₁(x₁,x₂,x₃)=(x₂+1,x₁,x₃),
    T₂(x₁,x₂,x₃)=(x₁+x₂,x₂,x₃+1).

Then L=span{e₁,e₂}, so there are four observable classes, indexed by the first two bits; the third bit is invisible. The shifts produce the signs in (14), but neither exposes e₃. The first map flips the sign of χ_e₁ after pullback; ignoring that sign gives the wrong dynamics even though it gives the right support space.

### Guard: a Fourier support list is not the observed function

Theorem 2 assumes the characters in S are individually recorded. If only the sum f=χ_e₁+χ_e₂ is recorded under identity dynamics on F₂², f has three values 2,0,−2 and identifies three classes. Its generated algebra has dimension three. Treating the two support labels as separately observed would incorrectly predict four classes. The general Theorem 1 handles the sum correctly; the label shortcut does not apply.

Likewise, affine shifts of a general sum can change interference among its terms. The shift-independence claim in Theorem 2 is scoped to individually available characters, not arbitrary observables with that Fourier support.

## 6. Tiny exact controls and reproducibility

The standard-library script `verify_observable_algebra.py` compares two independently organized constructions: rational row-space closure under pullbacks/products, and finite partition refinement. For binary cases it also computes label closure over F₂, verifies annihilator classes, Fourier inversion, and the signed affine moment update.

The complete declared universe consists of eight specified systems:

1. A five-state nonlinear map with four observable classes;
2. Identity on three states observed by values 0,1,2: linear dimension two, algebra dimension three;
3. The summed-character identity hostile above;
4. Identity on two bits with the two individual coordinate characters: linear dimension three, algebra dimension four;
5. The three-bit cyclic register: dimensions four and eight;
6. The two-map affine quotient above: four observable classes;
7. The zero map on three bits with one initially observed character: two classes, because the empty word is retained;
8. No nonconstant observed character on a three-bit cyclic register: one class.

Normal and optimized runs agree byte-for-byte:

    python -B verify_observable_algebra.py
    python -B -O verify_observable_algebra.py

The JSON receipt records **701 exact checks** and each case's dimensions and refinement depth. These finite controls support the proofs; they are not a census of all registers or finite systems. The earlier geometric packet and its audit remain unchanged. This addendum is ready for its separate independent review; no theorem ID is reserved and no Git write was made.
