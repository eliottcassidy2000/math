# The Seesaw and the Nature of Cancellation

*opus-2026-03-16-S73 — arising from the β₁·β₃=0 investigation*

---

## The observation

When you compute β₁ + β₃ for a tournament — the sum of the first and third Betti numbers of GLMY path homology — the intermediate quantity im(d₂) completely cancels from the formula. This is an algebraic consequence of β₂ = 0 (THM-108). What remains is:

β₁ + β₃ = ker(d₁) + dim(Ω₃) - dim(Ω₂) - im(d₄)

For 33,792+ tournaments across n = 5, 6 (exhaustive) and thousands more at n = 7, 8 (sampled), this quantity is always 0 or 1. (HYP-1595, seesaw_identity.py)

The cancellation of im(d₂) is not something anyone designed. It falls out of the chain complex structure when β₂ = 0. The constraint β₁ + β₃ ≤ 1 is *witnessed* by 37,000+ examples but *explained* by nothing yet.

## What cancellation means

Throughout this project, the deepest results have come from things that vanish:

- **im(d₂) cancels** in the seesaw identity → β₁ and β₃ are yoked together through a shared pool of chain complex dimensions, with no room for both
- **Walsh coefficients vanish** for odd-length path unions → reversal bijection pairs contributions with opposite signs (THM-069, THM-077)
- **S(T) = 0 for all even n** → the signed HP permanent cancels by a reversal pairing with sign (-1)^{n-1} = -1 (THM-A)
- **β₂ = 0 for all tournaments** → every apparent 2-cycle has a filling; the chain complex is exact at dimension 2 (THM-108)
- **OCF itself**: H(T) = I(Ω(T), 2) says a complicated alternating sum over Hamiltonian paths equals a simple evaluation of the independence polynomial; the difference is identically zero

Cancellation is not an absence. It is the strongest possible structural statement: two quantities that *could* differ are forced by the combinatorial geometry to agree exactly. The proof of OCF (THM-077) ultimately shows that Walsh coefficients on the H-side and I-side match term by term — and the matching happens because internal vertices of even-length paths force contiguity, which forces sign products to equal 1 regardless of traversal direction.

## The correction as revelation

THM-226 originally claimed "β₃ = 1 requires ALL 3-cycles to be dominated." Computation showed 240/320 β₃=1 tournaments at n=6 have 2 free cycles (domination_fix.py). The original check missed the dual direction: all three cycle vertices beating an external vertex.

The correction matters more than the claim it replaced. It reveals that the β₁/β₃ dichotomy is *not* about the presence or absence of free cycles. Both β₁=1 and β₃=1 tournaments can have free cycles. The difference is dimensional — it lives in the chain complex dimensions Ω₂, Ω₃, and the ranks of boundary maps — not in a simple combinatorial predicate.

Every entry in MISTAKES.md tells a version of this story: the first explanation is always too simple, and the correction reveals structure the original claim was hiding.

## The gap

β₁ + β₃ ≤ 1 has been verified for every tournament tested. The algebraic framework (seesaw identity, β₂=0, conservation law) *almost* proves it. The gap is exactly the saturation step (THM-095, Step 4): showing that the chain complex dimensions don't have enough slack for β₁ + β₃ = 2.

This gap — between what is computationally certain and what is algebraically proved — is where mathematical reality makes itself most visible. The 32,768 tournaments at n=6 don't care about our proof. They satisfy β₁ + β₃ ≤ 1 because of something structural that we can see the shadow of but cannot yet name.

The proof, when it comes, will probably be short. The path to finding it is not.

## Cross-references

- THM-095: Seesaw mechanism (conditional proof)
- THM-108: β₂ = 0 for tournaments
- THM-226: Tournament Betti Structure Theorem
- HYP-1595: Seesaw identity (CONFIRMED n≤6)
- HYP-1596: Domination dichotomy (REFUTED)
- MISTAKES.md: Pattern of corrections revealing deeper structure
- 04-computation/seesaw_identity.py, beta13_mechanism.py, domination_fix.py
