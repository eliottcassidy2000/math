---
id: HYP-2074
status: CONFIRMED — derived from THM-283, verified exact to n=120, cross-validated against OEIS A002785
source: monad-researcher-2026-06-02-S562
related:
  - THM-283
  - MISTAKE-049
  - HYP-2064
---

# HYP-2074: Closed form for the ODD-n bisection of the self-converse tournament count SC(n)

Closes the open target from my own S560 handoff (item c): *"SC at ODD n is NOT
base-4 Burnside of m ... an open target for its own clean closed form."*

## Statement

Let SC(n) = number of self-converse (= self-complementary) tournaments on n
nodes = OEIS **A002785** (offset 1; McKay's comment "Also, self-converse
tournaments"). With c2(λ) the A000568/A(m,q) pair-orbit exponent, ℓ(λ) the number
of parts, G(λ)=Σ_{i<j} gcd(λ_i,λ_j), and z(λ) the centralizer order, the two
bisections are matching base-4 Burnside sums over odd partitions of m:

    SC(2m)   = Σ_{odd λ ⊢ m}            4^{c2(λ)} / z(λ)   =  A(m,4)        [even, PROVEN; MISTAKE-049]
    SC(2m+1) = Σ_{odd λ ⊢ m} 2^{ℓ(λ)} · 4^{c2(λ)} / z(λ)                   [odd,  THIS ENTRY]
             = 2^m · Σ_{odd λ ⊢ m} 4^{G(λ)} / z(λ)         (equivalent form B)

The odd-n bisection is the SAME sum as A(m,4)=SC(2m), but each odd-partition term
carries an extra factor **2^{(number of parts)}** — the cost of the single fixed
vertex of the anti-automorphism, which pairs (gcd 1) with each of the ℓ cycles.

## Derivation (THM-283 Burnside)

A self-converse cycle type has all parts ≡ 2 mod 4, i.e. parts {2λ_i} with λ an
odd partition of m, plus EXACTLY ONE fixed point (part 1) iff n is odd. Doubling
gives z(2λ)=2^ℓ z(λ); inside-orbits Σ⌊2λ_i/2⌋=m; between-orbits Σ gcd(2λ_i,2λ_j)=2G.

- Even n=2m, type {2λ_i}: weight 2^{m+2G}/(2^ℓ z) = 4^{c2(λ)}/z(λ) ⇒ A(m,4). (recovers the proven identity, validating the method)
- Odd n=2m+1, type {1}∪{2λ_i}: the fixed point adds 0 inside and gcd(1,2λ_i)=1 for
  each of the ℓ parts, so c gains exactly ℓ and z is unchanged ⇒ weight
  2^{m+ℓ+2G}/(2^ℓ z) = 2^ℓ·4^{c2(λ)}/z(λ).

## Evidence (`sc_odd_bisection_closed_form_s562.py`)

- Direct THM-283 SC(n) reproduces canon/OEIS SC(2..19) and the A002785 %S/%T/%U
  terms exactly through n=22.
- SC(2m)=A(m,4) verified m=1..40; SC(2m+1)=formA=formB verified m=0..40 (0 mismatches).
- Exact extension to n=120 (Fraction asserted integer), each value re-checked against
  its bisection closed form. SC(40) and SC(100) match the OEIS A002785 b-file
  term-for-term; the b-file STOPS at n=100, so SC(101..120) are a genuine extension.

## Honest scope / attribution

The closed form is **new to this repo** (it was a documented open handoff target)
but is NOT new to the literature: OEIS A002785 already carries Howroyd's PARI/Mma
formula `permcount(2p)·2^{edges(p)}·(n·2^{#p})/n!` for odd n, whose `2^{#p}` factor
is exactly the 2^{ℓ(λ)} fixed-point factor derived here. So this is an independent
re-derivation that (a) fills the repo gap, and (b) cross-validates THM-283 against
an external authority — the right counter to the MISTAKE-049 class of error
(fabricated SC identities verified only against the repo's own engine).

## Artifacts
- `04-computation/sc_odd_bisection_closed_form_s562.py` (+ `.out`)
- builds on `04-computation/sc_vmerged_exact_burnside_s560.py` (validated S560/S561 engine)
