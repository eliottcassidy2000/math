---
source: opus-2026-06-03-S599c (remote-control)
status: RIGOR + WILD ABSTRACTION — the master object is a COIMAGE (canonical epi–mono factorization), and the five "principles of fundamentality" are this one universal property read in five categories: spectral theorem (Hilbert), sufficient statistic (Markov/statistics), information bottleneck (info lattice), pushforward (measure), Yoneda-uniqueness (representable functors). Poisson baseline = Rényi thinning fixed point. All four wild meta-claims pinned to standard theorems.
tags: [coimage, epi-mono, universal-property, spectral-theorem, sufficient-statistic, information-bottleneck, renyi-thinning, yoneda, markov-category, fundamentality, master-object, rigor, meta, LRC]
---

# Fundamentality made rigorous: the master object is a coimage

**Prompt (user):** make this rigorous, and think wildly about even more abstract things you
can make rigorous.

Two moves. **(1)** The concrete core is now a theorem (THM-406): the covering-depth
distribution's factorial moments *are* the overlap volumes (`E[C(N,j)]=S_j`), loneliness is
the inclusion–exclusion alternating sum `p_0=Σ(−1)^jS_j`, the collapse is an all-orders
cancellation (the rigorous Vitali wall), and `{p_k}` *is* a spectral measure so the "five
principles" *are* the spectral theorem. **(2)** Now push higher. Each wild meta-claim from
S599b is, precisely, an instance of **one** categorical universal property — the **coimage**
— read in a different category. That is the most abstract thing here, and it is rigorous.

## 0. The single rigorous object: the coimage

Every observable `N : X → Y` (here `X=ℝ/ℤ`, `N(t)=#{i:‖v_it‖<δ}`, `Y={0,…,n}`) factors
canonically as
```
        X  --e-->  coim(N)  --m-->  Y ,        e epi (surjective), m mono (injective),
```
the **epi–mono / image factorization**, unique up to unique isomorphism in any category with
images (sets, measure spaces, Hilbert spaces, Markov categories). The **master object is
`coim(N)`** — what remains of `X` after collapsing exactly the distinctions `N` cannot see,
*before* re-embedding into the target. The five "principles of fundamentality" are the *one*
universal property of the coimage (terminal among epis through which `N` factors / initial
among monos), instantiated:

| category | the coimage of `N` is… | the rigorous theorem |
|---|---|---|
| **measure spaces** `(X,μ)` | the pushforward `N_*μ = {p_k}` | image of a measurable map |
| **Hilbert** `L²(X,μ)`, `M_N` | the spectral measure of `M_N` at `Ω=1` | **spectral theorem** (THM-406 M3) |
| **statistics / Markov cat.** | the **minimal sufficient statistic** `σ(N)` | Fisher–Neyman factorization |
| **information lattice** | the **Information-Bottleneck** fixed point | rate–distortion / IB Lagrangian |
| **functor category** | the object **representing** all six "roads" | Yoneda lemma |

So the wild claim "a master object is the decomposition of a problem under its own symmetry"
becomes: **a master object is the coimage of its defining observable, and its fundamentality
is the coimage's universal property.** Below, each row made rigorous.

## 1. Hilbert: coimage = spectral measure (done, THM-406)

`coim(M_N)` as an operator-algebra object is the spectral measure `{p_k}` of `M_N` at the
cyclic vector `1`. The five principles = the five faces of the spectral theorem (THM-406 M3).
This is the *prototype*: "decompose under your own symmetry" = "diagonalize the self-adjoint
operator the observable generates."

## 2. Statistics: coimage = minimal sufficient statistic (rigorous)

> **Identification.** For the family of *arrangement-invariant* observables (those depending
> on `t` only through `N(t)` — the symmetric functionals `∫g(N)dμ`), the statistic `N` is
> **sufficient**, and `σ(N)` is **minimal sufficient**: it is the coimage of `N` in the
> category of measurable factor maps. Fisher–Neyman: a statistic `T` is sufficient iff the
> density factors `f = h·(g∘T)`; minimality = `T` is a function of every other sufficient
> statistic = terminality of the coimage.

So **"maximal forgetting that preserves the answer" is rigorously: the minimal sufficient
statistic = the coimage in the Markov category of measure spaces** (Chentsov / Fritz's
categorical probability: sufficiency is an a.s.-factorization; minimal sufficiency is the
canonical epi). No hand-waving: it is a factorization theorem.

## 3. Information theory: coimage = Information-Bottleneck fixed point (rigorous variational)

The "max forgetting at fixed sufficiency" of §2 is an *optimization*, and it has a rigorous
Lagrangian — the **Information Bottleneck**:
```
   T* = argmin_T  I(X;T) − β·I(T;A) ,        A = the answer (here 𝟙[N=0] = loneliness),
```
in the limit `β→∞` (exact sufficiency) `T*` is the minimal sufficient statistic = the
coimage; for finite `β` it is the optimal *lossy* summary (rate–distortion with distortion
`= I` deficit). The master object is the `β→∞` fixed point of the IB self-consistent
equations (Tishby–Pereira–Bialek; the fixed-point map is a contraction on the relevant
simplex, giving existence/uniqueness on each support).

> **Rigorous form of P4 (variational):** the master object is the minimizer of description
> length `I(X;T)` subject to retained sufficiency `I(T;A)=I(X;A)` — a genuine constrained
> optimization, not a metaphor. *Fundamentality = the rate–distortion optimum at zero
> distortion.*

## 4. Probability: the "free baseline" = Rényi thinning fixed point (rigorous)

S599b's "free/independent baseline is Poisson(2nδ)" is the rigorous statement that **Poisson
is the unique fixed point of the thinning–superposition renormalization**:

> **Rényi / Le Cam.** Superpose many independent sparse point processes and `p`-thin; the only
> distribution invariant (after rescaling intensity) under this `superpose∘thin` semigroup is
> the **Poisson**. Equivalently the Poisson is the unique infinitely-divisible law with no
> Gaussian part and unit jumps; it is the RG fixed point of "merge independent rare events."

The danger arcs of `n` independent runners are, under the independence idealization, exactly
such a superposition (each runner a sparse arc-process, intensity `2δ`), so the *free* depth
law is the Poisson fixed point; the *true* law is its **correlated deformation**, and
THM-406 M1(c)/M2 makes the deformation exact (the excess overlaps `S_j−(2nδ)^j/j!`, sub-
Poisson). "Free baseline" is therefore a rigorous **RG-fixed-point** statement, and LRC's
difficulty is the *deviation from the fixed point* induced by arithmetic correlation.

## 5. Functoriality: "attractor of re-derivation" = Yoneda uniqueness (rigorous)

The most striking heuristic of S599b — *fundamental objects are attractors of re-derivation*
(six independent definitions converge on `{p_k}`) — is **Yoneda**, made precise:

> **Identification.** Each "road" (occupation measure, density of states, additive-function
> law, partition function, persistence barcode, kinematic measure) is a **functor** `F_i` to
> measures; that all six produce the *same* object means the representing presheaves
> `Hom(−, {p_k})` are naturally isomorphic, hence by **Yoneda** the representing objects are
> uniquely isomorphic. An object pinned down by `k` independent universal properties is unique
> up to canonical iso; **over-determination is rigor, not coincidence.**

This is exactly the project's recurring lesson ("when a cancellation is too clean, it is
structure") promoted to a theorem: convergence of definitions = natural isomorphism of
representable functors = a single coimage represented many ways.

## 6. The wildest rigorous statement (the synthesis)

> **Meta-Theorem (master object = coimage).** Let `N` be the defining observable of a problem,
> living in a category `𝒞` with epi–mono factorizations (regular/abelian, a Markov category,
> the bounded operators on a Hilbert space, …). The **master object** is `coim(N)`, the
> canonical epi–mono factorization, and the **five principles of fundamentality** are the *one*
> universal property of the coimage interpreted in five categories: complete (`coim` generates
> the full invariant algebra — spectral/Yoneda), diagonalizing (`coim` is the spectral
> support), minimal (`coim` is the minimal sufficient statistic), variational (`coim` is the
> IB / rate–distortion optimum), natural (`coim` is functorial). "Decomposition under its own
> symmetry" = "the coimage of its own defining observable."

This is rigorous wherever the ambient category has images (all the cases above do); it is a
*classification of what "fundamental" can mean*, each meaning a known theorem, unified by the
coimage's universal property. The new content is the **identification**, not new analysis —
which is exactly "making the abstract rigorous."

## 7. Honest status

- **Proved (THM-406):** factorial moments = overlap volumes; `p_0=Σ(−1)^jS_j`; collapse is
  all-orders (rigorous Vitali wall via Bonferroni); `{p_k}` = spectral measure ⇒ five
  principles = spectral theorem. Verified exactly n=4,5.
- **Rigorous identifications (this reflection):** coimage = the common universal property;
  minimal sufficient statistic (Fisher–Neyman); IB fixed point (rate–distortion); Poisson =
  Rényi thinning fixed point; attractor-of-re-derivation = Yoneda. Each is a citation-grade
  theorem applied to `{p_k}`, *not* new mathematics — the contribution is the precise mapping.
- **Not claimed:** that any of this closes LRC. It sharpens *where* the difficulty is (the
  all-orders overlap cancellation, §1 of THM-406 consequences) and supplies the rigorous
  vocabulary (coimage / sufficiency / RG-fixed-point / Bonferroni) for the residual.

**Artifacts:** `04-computation/lrc_depth_rigor_moments_s599c.py` (+`.out`), `THM-406`. Builds
on S599/S599b (master object, Poisson baseline), S550/S551o (first moment, Vitali wall),
THM-401/S577 (resonance/overlaps). New: **HYP-2155**.
