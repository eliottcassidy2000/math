# Native spectral boundary and an atomic Abel certificate

Session106, 2026-10-29. **PROVED elementary transfer/obstruction theorems; finite exact checks.** No claim of Collatz convergence or literature priority. The main new synthesis is a sharp refresh-rate threshold, a native ℓ¹ spectral obstruction, and a fully source-weighted resolvent interface.

**Signed-residual extension, Oct29:** [proof-target refinement, §4](../../proof_targets_refinement_20261029/README.md)
allows a summable signed trial flow with a quantified residual in place of
an exact positive supersolution. It retains negative residual credit and
proves monotonicity under paid native refinement, with the original source
atoms and integer clock unchanged. Independent arithmetic tail bounds,
and hence universal Collatz entry, remain open.

## 1. Fix the operator before using its spectrum

Let X be countable, H the genuine terminal set, and T deterministic on X\H. Set Y=X\H and define the column substochastic killed operator

    (Qv)(y) = sum_{x in Y:T(x)=y} v(x),  y in Y.

Transitions into H are discarded, not continued or refreshed. Choose a probability w on X with every w(x)>0 and restrict it to Y when multiplying Q. Let τ be first hitting time, including τ=0 on H. Then exactly

    A_k = ||Q^k w||₁ = sum_x w(x) 1{τ(x)>k}.

All sums are nonnegative; no change of source law is allowed. The inherited Session104 atomic theorem says A_k→0 iff every source terminates. Equivalently Q is strongly stable on ℓ¹(Y): first prove decay on each point mass, then use finite-support approximation and ||Q||≤1. It does NOT say ||Q^k||→0.

For accelerated Collatz take X the positive odd integers, H={1}, and T(n)=(3n+1)/2^{v₂(3n+1)}. A particularly transparent full-support law is

    w(2j−1)=1/[j(j+1)], j≥1;
    sum_{j>J} w(2j−1)=1/(J+1).

These are source weights, not stationary weights, not Haar cylinder masses, and not weights reinitialized after a move.

## 2. A uniform native ℓ¹ spectral gap is impossible

**Theorem.** For this killed accelerated Collatz operator, ||Q^k||₁→₁=1 for every k≥0, and its ℓ¹ spectral radius is 1. This holds without assuming or refuting convergence.

For k≥1, choose n₀=2^{k+1}−1. For 0≤i≤k the actual orbit is

    n_i=2^{k+1−i}3^i−1.

For i<k, 3n_i+1=2(2^{k−i}3^{i+1}−1), with the bracket odd, so the exact valuation is 1. All these states exceed 1. Hence Q^k δ_{n₀} has mass 1. Substochasticity gives the reverse norm inequality. The spectral-radius formula yields radius 1.

Thus a bound ||Q^k||≤Cρ^k with finite C and ρ<1 is impossible on ordinary ℓ¹, even if every orbit eventually stops. A weighted norm or source-specific estimate may evade this theorem, but its weights and atom comparison must be proved. The ordinary countdown map T(j)=j−1, H={0}, already shows strong stability coexisting with spectral radius 1.

## 2a. The entire closed unit disk is already forced by terminating sources

The [source-response construction](../observability/SOURCE_RESPONSE_RANK.md), §1, supplies a strictly increasing inverse chain x₀=1,x₁=5,x₂=13,... with T(x_j)=x_(j−1) and exact τ(x_j)=j. For completeness, given odd y not divisible by3, choose x=(2^a y−1)/3, using a=4 at y=1; otherwise choose a from {2,4} if y≡1 mod3 or {3,5} if y≡2 mod3, rejecting a candidate divisible by3. Both candidate numerators cannot vanish mod9, since multiplication by4 sends residue1 to4. The chosen x is integral, odd, greater than y, not divisible by3, and has exact valuation a. Thus the chain exists without a convergence assumption.

Complexify ℓ¹(Y). For every complex |λ|<1, the nonzero vector

    v_λ = sum_{j≥1} λ^(j−1) δ_(x_j)

belongs to ℓ¹, with norm 1/(1−|λ|). Boundedness of Q permits termwise application. The first basis vector is killed and all later ones shift down the chain, so

    Q v_λ = sum_{j≥2} λ^(j−1) δ_(x_(j−1)) = λ v_λ.

This includes λ=0, for which v₀=δ₅. Therefore the open unit disk is contained in the point spectrum. The spectrum is closed and is contained in the closed unit disk by ||Q||≤1 (outside it the Neumann resolvent converges). Hence **the spectrum of native killed accelerated Collatz on complex ℓ¹ is exactly the closed unit disk**. No assertion is needed about which boundary points are eigenvalues. This entire spectral obstruction already occurs on the closed invariant subspace supported on the displayed, individually proven-terminating chain.

For real 0<λ<1, normalize to the probability p_λ=(1−λ)v_λ. Then

    ||Q^k p_λ||₁ = λ^k,
    E_(p_λ) τ = 1/(1−λ),
    (1−z)||(I−zQ)^{−1}p_λ||₁ = (1−z)/(1−zλ).

Thus arbitrary rates approaching1 occur for source laws supported entirely on known finite trajectories. These laws are NOT full support on all positive odd integers, and no failure of universal convergence follows. For a finite truncation v_(λ,N), the exact residual is

    Q v_(λ,N) − λ v_(λ,N) = −λ^N δ_(x_N).

The checker verifies this using exact rational real/imaginary parts, independently reconstructing the inverse chain and native valuations. The infinite eigenvector theorem follows from absolute summability, not finite residual numerics. This strengthens the norm-one obstruction and complements the unbounded source-response rank; it does not rule out source-specific weighted certificates.

### Finite cutoff spectra can be maximally misleading

On span{δ_(x₁),...,δ_(x_N)}, Q is the nilpotent backward shift: Q^N=0 but Q^(N−1)δ_(x_N)=δ_(x₁)≠0. Every finite restriction therefore has spectrum {0}, whereas their invariant infinite-ray limit has spectrum the entire closed unit disk. Finite cutoff eigenvalues cannot certify a global native spectral gap. Even for a fixed cutoff, ||Q^k||=1 through k=N−1. At any |λ|=1, the truncated eigenvector above has norm N and residual norm1, so normalization gives approximate boundary eigenvectors of residual norm1/N. These observations do not invalidate genuine finite-source hitting certificates; they identify the missing uniform control in a spectral-limit argument.

### Arbitrary vanishing survival tails on the same terminating support

More generally let 1=a₀≥a₁≥a₂≥...→0. Put mass p_j=a_(j−1)−a_j at x_j. These nonnegative masses telescope to1, and since τ(x_j)=j,

    Pr_p(τ>k)=sum_{j>k}p_j=a_k.

Thus every such prescribed survival tail is realized by a source law on this explicit terminating chain. Zero p_j are allowed. In particular p_j=1/[j(j+1)] gives A_k=1/(k+1) and infinite mean lifetime, although every supporting source has its complete finite trajectory supplied. This weight law uses the chain index j, NOT the earlier law w(2j−1)=1/[j(j+1)] on all positive odd integers. Arbitrarily slow decay for arbitrary chain-supported priors does not settle decay for that fixed, full-support original-source prior. It rules out a uniform decay class across arbitrary priors, not a specific scalar generating function claim for the chosen Collatz source law.

## 3. Refresh can manufacture convergence; its cost is exact

Let R be any positive substochastic kernel on the SAME active carrier, and 0<α<1. Set L=(1−α)Q+αR. Repeated positivity gives

    L^k w ≥ (1−α)^k Q^k w entrywise,
    A_k ≤ (1−α)^{−k} ||L^k w||₁.                 (1)

This uses no commutativity. Therefore a bound ||L^k w||₁≤Cρ^k transfers to native termination if ρ<1−α. More generally it suffices that (1−α)^{−k}||L^k w||₁→0. For a given source x, any k with C[ρ/(1−α)]^k<w(x) proves τ(x)≤k.

**Sharp hostile.** One active state b is a native fixed point, Q=[1]. Refresh sends b directly to H, R=[0]. Then L=[1−α] has exponential killed decay and finite randomized expected lifetime 1/α, while b never stops natively. The threshold ρ=1−α cannot be relaxed to ≤. Adding more transient states changes nothing about the obstruction. Random refresh has paid precisely for the event that no refresh has occurred.

The Session105 gap 3α/n concerns mixing of a different, finite reversible tournament process toward its stationary law. It is neither killed survival decay nor a bound beating the no-refresh cost. A stationary spectral gap needs additional hitting-set hypotheses even to imply absorption after killing. Those hypotheses still concern the randomized process, not deterministic Collatz.

## 4. Exact Abel resolvent: a finite certificate at one discount

For 0<z<1 define on ℓ¹(Y)

    G_z = (I−zQ)^{−1} w = sum_{k≥0} z^k Q^k w,
    D_z = (1−z)||G_z||₁.

The Neumann series always exists since ||zQ||<1. Mere existence says nothing about termination. Tonelli and the geometric sum give the exact identity

    D_z = sum_x w(x)(1−z^{τ(x)}),                 (2)

where z^∞=0. Therefore lim_{z↑1}D_z is exactly the total weight of nonterminating sources. Universal termination is equivalent to D_z→0, but need NOT imply bounded ||G_z|| as z↑1: that stronger property is finite weighted mean lifetime, the Session104 budget condition.

**Atomic extraction at one z.** If a proved upper bound D_z≤d satisfies d<w(x), then τ(x)<∞. More sharply every integer m with

    w(x)(1−z^m)>d

satisfies τ(x)<m. This is an exact rational test when z,d,w(x) are rational; logarithms and limiting exchanges are unnecessary. The upper bound may be large relative to other atoms and still certify the chosen x.

**Noncircular positive supersolution.** Supply a nonnegative u_z∈ℓ¹(Y), verified independently of unknown hitting times, such that

    u_z ≥ w + z Q u_z.                            (3)

Iterating (3) proves u_z≥sum_{k=0}^N z^k Q^k w and hence u_z≥G_z. Consequently (1−z)||u_z||₁ is a valid d. A family with (1−z)||u_z||₁→0 proves universal termination. This is a genuine certificate format: one must give a formula/finite description for u_z, prove its incoming-flow inequality for every state, and certify its total mass. Defining u_z through unknown lifetimes or simply restating G_z with an unproved improved norm is not a solution.

An undiscounted summable u≥w+Qu is a stronger sufficient condition giving the finite-mean budget. With fixed polynomial source weights it may be obstructed by expanding words; none is asserted here.

### Finite computations and their exact missing mass

With the telescoping law above, simulate only sources j≤J, with horizon K. For each resolved source use 1−z^τ; for each unresolved source use 1; add the exact omitted tail 1/(J+1). The result is a rigorous rational upper bound on D_z. For survival mass the corresponding upper bound is

    A_K ≤ sum_{j≤J, τ(2j−1)>K} 1/[j(j+1)] + 1/(J+1).

A lower bound drops the tail. There is no hidden zero-density approximation. If the selected source is unresolved and included, its contribution 1 already costs its entire atom, so this crude Abel certificate cannot certify it. If it is omitted, the tail likewise contains its atom. Finite simulation alone does not create new convergence information; improvement requires a separately proved tail majorant such as (3).

## 5. Proposal clocks preserve paths, not gaps

A fixed-probability lazy native clock P=(1−p)I+pQ on active states, 0<p≤1, inserts waits but does not refresh the integer. Its random lifetime starting from a finite native lifetime τ has mean τ/p; an infinite native lifetime remains infinite. Exact occupation resolvents obey

    I−zP = [1−z(1−p)] [I−ζQ],
    ζ = zp/[1−z(1−p)],
    (1−z)||(I−zP)^{−1}w||₁ = (1−ζ)||(I−ζQ)^{−1}w||₁.  (4)

Thus compatible fixed clocks are a reparameterization of the same Abel obstruction, not a source of decay. With state-dependent acceptance p(x)>0 every fixed native edge is eventually accepted almost surely, but there need be no uniform physical-time bound; p=0 at a source breaks even almost-sure equivalence. The tournament all-permutations clock repairs reversibility because each guarded proposal is an involution. Native Collatz transitions do not acquire inverse edges from waiting.

## 6. Why Weyl harmonics and average projections do not supply this bridge

Session105 has degree-changing tournament maps D,E with [D,E]=I, and an exact resampling kernel ED/n. There is no supplied source-faithful Collatz realization of those maps. Equality of a formal commutator, sequence dimensions, pole locations, or mean drifts cannot identify the respective killed operators. Signed harmonic projections in particular need not preserve the positive cone or the source-atom lower bound in (2).

The finite lumpability obstruction is already visible without number theory. Let active states a,b be aggregated into one cell, with T(a)=H and T(b)=b. Uniformly resampling a representative on each visit gives the scalar killed kernel [1/2] and apparent decay 2^{−k}; actual source b survives forever. The original source law (half on each active point) has A₀=1, A_k=1/2 for k≥1. Endpoint averaging has silently switched dynamics. Exact native quotienting requires identical killed block transition totals for every representative; deterministic transitions force the terminal/nonterminal distinction to survive.

This is the same structural warning as Session105's authenticated strip-extension obstruction, but it is not claiming an isomorphism between the two examples.

## 7. Provenance and concrete next target

Read inputs: Session105 `updown/UPDOWN_HARMONIC_AUDIT.md`, `metagraph/WEIGHTED_NATIVE_RECURSION.md`, `harmonic/HARMONIC_GF.md`; Session104 `descent/POSITIVE_ATOMS_AND_INTEGER_DESCENT.md`, `positive/POSITIVE_INCIDENCE_AND_ATOMS.md`; Session100 `drift/DRIFT_TARGET.md`, `controllers/RETURN_CONTROLLER_TRAP.md`.

Inherited: fixed positive atoms, source-preserving stopped clocks, finite-mean rank extraction, Haar/integer distinction, expanding valuation-one words, and tournament refresh spectrum. Proved here: their sharp native spectral and mixed-kernel consequences, exact discounted certificate packaging, and explicit hostile controls. All infinite statements have proofs above; finite checks are not evidence for an unproved infinite quantifier.

**Next arithmetic obligation:** construct an explicitly described summable u_z satisfying (3) for actual stopped Collatz, with independently certified mass o((1−z)^{−1}), or supply a weaker source-specific mass bound sufficient to cross one new atom threshold. The ordinary ℓ¹ uniform-gap route is ruled out; the atomic strong-stability route remains open.

Reproduce with `python3 check_spectral.py` and `python3 -O check_spectral.py`. Standard library only; checks use explicit exceptions, not removable assertions.
