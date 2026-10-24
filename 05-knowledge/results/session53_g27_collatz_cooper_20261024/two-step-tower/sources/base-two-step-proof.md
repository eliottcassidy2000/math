# A same-precision two-step Cartier law for the weight-one s7 root

Status: PROVED CONDITIONAL ON THE EXPLICIT INHERITED INPUTS BELOW. The new reduction is a proof, not an extrapolation of the finite tests. Upgrading the inherited modular inputs requires their independent audit. No claim of external literature novelty is made.

Normalize G(0)=s(0)=1. Let G(t)=sum c(n)t^n, s²=(1+t)(1−27t), H=G/s, and C_p(sum a(n)t^n)=sum a(pn)t^n. Then

**Conditional theorem.** For every prime p and every n≥0,

c(p²n) ≡ c(n) (mod p²), equivalently C_p²G ≡ G (mod p²).

For split primes and p=7, this is an immediate iteration of inherited Theorem D. The content here is the inert complement p≡3,5,6 modulo 7. It is emphatically not the open strengthening c(np^r)≡c(np^(r−1)) modulo p^(2r), which only makes sense with that one-step target at split/ramified primes.

## Inherited inputs, explicitly separated

Source: cooper_s7_weight_one_root_frobenius_20261018.md (2026-10-18), sections 1, 5.7, 6; this artifact's manifest records its hash.

1. t(q), θ(q), and s satisfy G(t(q))=θ(q) and dt/t=θ²s dq/q. Put f=θ³s. Thus c(n)=CT_q(f t^(−n)), while [t^n]H=CT_q(θ³t^(−n)).
2. The exact identity f=1−7A1/8−49A2/8, where A1=sum_{n≥1}(sum_{d|n}χ(d)d²)q^n and A2=sum_{n≥1}(sum_{d|n}χ(n/d)d²)q^n, χ(d)=(d/7). These are weight-three Eisenstein series, with the constant term on the A1 component.
3. The excellent lift T=t(q^p)=t^p(1+p ε_p(t)), ε_p∈t Z_p[[t]]. For inert p≥5, inherited B and the proof of E provide H≡H_p G(t^p) mod p and ε_p H_p≡P_p mod p, where deg P_p≤d_p−1, P_p(0)=0, d_p=floor(p/3)+floor((p+1)/3)<p. For p=3 they provide ε_3≡0 mod3.
4. Inherited Theorem D handles split p and p=7.

## A new useful weight-three companion identity

Define f*=1−7A1/8+49A2/8. Then exactly

f*=θ³(1−3t/4).

Indeed, both sides are holomorphic weight-three forms on Γ0(7) with character χ: θ³t=(η(τ)η(7τ))³ is the weight-three cusp form; f* differs from the inherited holomorphic f by the Eisenstein form 49A2/4. The Sturm bound is floor(3[SL2(Z):Γ0(7)]/12)=2. Both sides have coefficients (1,21/4,105/4) through q². This proves equality. The independent integer check multiplies the identity by 8 and verifies through q^80.

For p≠7, the divisor-sum identities, with A_i(q^p) interpreted as V_p A_i, give

U_p A1=(1+χ(p)p²)A1−χ(p)p² V_p A1,
U_p A2=(χ(p)+p²)A2−χ(p)p² V_p A2.

Consequently for inert (necessarily odd) p,

U_p f≡f* (mod p²), and U_p²f≡f (mod p²).

All inert primes are odd, so the denominators 4 and 8 are p-adic units. These congruences include constant coefficients. The latter follows by applying U_p twice to the two Eisenstein components, whose respective eigenvalues modulo p² are 1 and −1. No assumption that U_p f≡f is made.

## The two-step constant-term calculation

Write T2=t(q^(p²)). Substitution of q^p in the definition of ε gives

T2=t(q^p)^p(1+p ε_p(t(q^p)))
   ≡t^(p²)(1+p ε_p(t(q^p))) (mod p²).

Here (1+p ε)^p≡1 mod p². This is valid for Laurent powers as well, after removing the common t-adic monomial. For n≥0,

c(p²n) ≡ CT_q(f T2^(−n))
          +pn CT_q(f ε_p(t(q^p)) T2^(−n)) (mod p²).

Use CT(A(q)B(q^p))=CT((U_p A)(q)B(q)); all Laurent products involved have finite negative tails. The first term is CT((U_p²f)t^(−n))≡c(n). In the correction term apply U_p once, and reduce its remaining factors modulo p:

CT(f ε_p(t(q^p))T2^(−n))
 =CT((U_p f) ε_p(t) t(q^p)^(−n))
 ≡CT(f* ε_p(t)t^(−pn))
 =[t^(pn)] H(t)(1−3t/4) ε_p(t) (mod p).

For p=3 this is zero because ε_3≡0. For inert p≥5, inherited input 3 changes the series in the final line to

(1−3t/4)P_p(t)G(t^p) (mod p).

The polynomial (1−3t/4)P_p has degree at most d_p<p and constant term zero. Therefore it has no terms of degree divisible by p, and neither does its product with G(t^p). The correction vanishes. This proves the inert claim and hence the all-prime theorem under the listed inputs.

## Boundaries and failed shortcuts

- C_pG≡H modulo p² is false: c(3)=336 while [t]H=15, difference321≡6 mod9. The proof never asserts it.
- U_p f≡f modulo p² is false at inert primes; the companion f* is essential.
- At p=3, c(9)−c(1)≡18 modulo27, so the all-prime two-step statement cannot uniformly be strengthened to modulus p³. The precise integer value is in normal.json.
- The new polynomial multiplier has degree one. This is exactly why the inherited bound deg P_p≤d_p−1 is sufficient; an unspecified formal-series error would not suffice.

## Finite verification and search scope

verify.py computes c(0)..c(10000) with exact integers and checks division at every recurrence step. For every prime p≤100 it checks every positive n with p²n≤10000. It independently checks the companion identity through q^80 and the exact divisor-sum Hecke identities in that range. It contains no Python assert statements. Normal and python -O outputs are byte-identical, both PASS.

A text search in the inherited clean-repo .md and .gp files for c(np^2), c(p^2n), two-step Cartier and U_p^2 f found only the existing root note's stronger open split-prime lift, not this same-precision inert assertion. This bounded repository search does not establish literature novelty.
