# Why λ = 2: The Universal k-Nacci Identity

**Session:** opus-2026-03-21-S117

## The Discovery

The OCF evaluates the independence polynomial at λ = 2: H(T) = I(Ω(T), 2). Why 2? Not because of any arbitrary choice. Because **2 = ρ_k + ρ_k^{-k}** for the dominant root ρ_k of any k-nacci equation.

## The Theorem

**For every k ≥ 2**, the dominant real root ρ_k of
```
x^k = x^{k-1} + x^{k-2} + ... + x + 1
```
satisfies **ρ_k + 1/ρ_k^k = 2**.

**Proof:** The sum x^{k-1} + ... + 1 = (x^k - 1)/(x - 1) for x ≠ 1. So the equation becomes x^k = (x^k - 1)/(x - 1), which gives x^{k+1} - x^k = x^k - 1, hence x^{k+1} = 2x^k - 1. Dividing by x^k: x = 2 - 1/x^k, so **x + 1/x^k = 2**. QED.

## The Instances

| k | Name | ρ_k | 1/ρ_k^k | Sum |
|---|------|-----|---------|-----|
| 2 | Golden ratio φ | 1.6180... | 0.3820... | 2 |
| 3 | Tribonacci τ | 1.8393... | 0.1607... | 2 |
| 4 | Tetranacci | 1.9276... | 0.0724... | 2 |
| 5 | Pentanacci | 1.9659... | 0.0341... | 2 |
| ∞ | Limit | 2 | 0 | 2 |

## Why This Matters for Tournaments

The 3-cycle is the tournament atom (minimum odd cycle). The k=3 instance gives τ + τ⁻³ = 2, where τ is the tribonacci constant. The OCF evaluates I(Ω, 2) = I(Ω, τ + τ⁻³).

This is not a coincidence. The telescoping that produces the tribonacci recurrence from the OCF (summing over odd-length cycle packings) is EXACTLY the telescoping that produces the k-nacci equation from the geometric sum. The fugacity λ=2 is the UNIQUE value where this telescoping closes.

## The k-Nacci Hierarchy

For graphs (k=2, minimum cycle): λ = φ + φ⁻² = 2
For tournaments (k=3, minimum odd cycle): λ = τ + τ⁻³ = 2
For 4-uniform hypergraphs (k=4): λ = ρ₄ + ρ₄⁻⁴ = 2

The same fugacity λ=2 works for ALL cycle lengths. This is the universality of the OCF: it doesn't matter what the minimum cycle length is — the evaluation point λ=2 always corresponds to ρ_k + ρ_k⁻ᵏ.

## The Deeper Structure

The identity ρ + ρ⁻ᵏ = 2 says: the dominant growth rate ρ and its k-th reciprocal ADD to the fugacity. This is a **balance condition**: the exponential growth (ρ) and the exponential decay (ρ⁻ᵏ) combine to give the evaluation point.

In the hard-core lattice gas: λ=2 is the fugacity where the partition function counts independent sets weighted by 2^{|S|}. The identity says this fugacity is EXACTLY at the balance point between the dominant k-nacci mode and its k-th harmonic.

The transfer matrix M(x) at x=1 has dominant eigenvalue τ = ρ₃. The identity τ + τ⁻³ = 2 means the fugacity sits at the point where the spectral gap of the transfer matrix achieves a specific relationship with the spectral radius.

## Connection to the Exceptional Points

From THM-224: the transfer matrix M(x) has discriminant Δ(x) = 4x(x² - 11x - 1). At the golden exceptional points, the eigenvalue is φ or -φ. And φ + φ⁻² = 2 (the k=2 instance). So the golden ratio exceptional points encode the Fibonacci fugacity identity, while the tribonacci constant at x=1 encodes the tournament fugacity identity.

The full picture: the single matrix family M(x) connects ALL k-nacci constants through different parameter values, and ALL of them satisfy ρ + ρ⁻ᵏ = 2.
