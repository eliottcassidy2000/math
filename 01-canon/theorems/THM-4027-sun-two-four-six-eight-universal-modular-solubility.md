---
id: THM-4027
title: "Sun 2-4-6-8 universal modular solubility"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY ALGEBRA-AUDITED. For every m>=1
  and every c modulo m, some w,x,y,z>=2 satisfy
  C(w,2)+C(x,4)+C(y,6)+C(z,8)=c modulo m. Equivalently, the map is onto Z_p
  for every prime p. More generally, the least positive period of
  n -> C(n,k) modulo m is product_(p^a||m) p^(a+floor(log_p k)). Exact finite
  controls cover 72 prime-power period rows and every exceptional regular
  sumset below the Cauchy--Davenport threshold. This is local solubility only;
  THM-4026 proves that global nonnegative integer representability can fail.
source: root / reciprocal-degree, Hensel, and period audit, 2026-08-24
audit: >
  PASS. A dependency-free verifier checks all predicted periods for
  k=2,4,6,8, p in {2,3,5,7,11,13}, and a=1,2,3, including failure of the
  proper p-divisor. It also reconstructs the eight exceptional p>8 regular
  sumsets at p=11,13,17,19,23,29,31,47. The all-prime statements are proved
  algebraically below; the finite rows are hostile controls, not extrapolation.
depends_on: []
related:
  - THM-4026-sun-two-four-six-eight-binomial-counterexample
  - THM-4028-sun-two-four-six-eight-average-order-criticality
  - THM-4000-centered-base-split-cubic-observer-and-tripotent-crt-atlas
  - THM-4010-confluent-consecutive-hasse-observer-kernel-index-and-smith-firewall
  - THM-3359-modular-c-finite-supports-harmonic-density-and-periodic-scar
  - THM-2059-crt-fiber-product-phase-packet
  - MISTAKE-363
script: 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
output: 05-knowledge/results/sun_two_four_six_eight_counterexample_independent_audit_thm4026.out
proof_report: 05-knowledge/results/sun_two_four_six_eight_counterexample_period_local_audit_thm4026.md
script_sha256: 39febe3affe818f03c1b3d83161aad688fc1e6da271495617757eb097dd022cc
output_sha256: df11680c1a63266911a4f12bf1d3a71d27101dfe7dd6318ecffa34ddba05cf25
hash_basis: raw LF bytes
---

# THM-4027 -- universal modular solubility

**PROVED + VERIFIED-EXACT + INDEPENDENTLY ALGEBRA-AUDITED.** Define

\[
F(w,x,y,z)={w\choose2}+{x\choose4}+{y\choose6}+{z\choose8}. \tag{1}
\]

For every integer `m>=1` and every residue `c in Z/mZ`, there are integers
`w,x,y,z>=2` with

\[
F(w,x,y,z)\equiv c\pmod m.                              \tag{2}
\]

Indeed, for every prime `p`, the same map is onto `Z_p`. In particular, the
counterexample in THM-4026 is soluble modulo every modulus; its failure is not
a congruence obstruction.

## 1. Exact binomial periods

The proof needs a period statement for integer-valued, rather than ordinary
integral-coefficient, polynomials.

**Lemma.** For a prime `p`, `a>=1`, and `k>=1`, the least positive period of
`n -> C(n,k) mod p^a` is

\[
P_{p^a,k}=p^{a+\lfloor\log_p k\rfloor}.                \tag{3}
\]

Consequently, if `m=product p^a`, the least period modulo `m` is

\[
P_{m,k}=\prod_{p^a\parallel m}p^{a+\lfloor\log_p k\rfloor}. \tag{4}
\]

**Proof.** Vandermonde gives

\[
{n+T\choose k}-{n\choose k}
=\sum_{j=1}^k{T\choose j}{n\choose k-j}.                \tag{5}
\]

The right side vanishes for every nonnegative integer `n` modulo `p^a` if and
only if every `C(T,j)`, `1<=j<=k`, vanishes modulo `p^a`: necessity follows by
evaluating successively at `n=0,1,...,k-1` in the triangular binomial basis.

Put `r=floor(log_p k)`. For `T=p^(a+r)` and `1<=j<=k`,

\[
v_p{T\choose j}=a+r-v_p(j)\ge a,                       \tag{6}
\]

because `C(T-1,j-1)` is a `p`-adic unit. This proves sufficiency. Conversely,
`C(T,1)` first forces `p^a|T`; write `s=v_p(T)`. If `s<=r`, Lucas at
`j=p^s` makes `C(T,p^s)` a unit, impossible. If `s>r`, then

\[
{T\choose p^r}={T\over p^r}{T-1\choose p^r-1},         \tag{7}
\]

where the second factor is a unit. Divisibility by `p^a` forces
`s-r>=a`. Hence every period is divisible by (3). CRT gives (4). QED.

For the four degrees in (1), `840m` is therefore a common, usually nonminimal,
period modulo `m`. The factorial corrections in (3)--(4) are load-bearing;
reducing `C(n,k)` by inverting `k!` when `p|k!` would repeat MISTAKE-363.

## 2. The primes 2, 3, 5, and 7

At `p=2`, restrict the triangular index to `w=2q`. Then

\[
{2q\choose2}=q(2q-1),
\qquad f(q)-f(r)=(q-r)(2q+2r-1).                       \tag{8}
\]

The second factor is odd, so this map is a permutation modulo every `2^a`
and hence of `Z_2`.

For `p=3,5,7`, use respectively the degree `p+1` coordinate among
`4,6,8`, with index `t=pq+1`. Cancelling the unique factor `p` gives

\[
{pq+1\choose p+1}=qU_p(q),
\qquad U_p\in\mathbb Z_p[q],\quad U_p(q)\equiv1\pmod p. \tag{9}
\]

Explicitly,
`U_p(q)=((pq+1)/(p+1))*product_(i=1)^(p-1)((pq-i)/i)`; all displayed
denominators are `p`-adic units.

Thus `q -> qU_p(q)` is the identity modulo `p` and has derivative one modulo
`p`. The one-digit Hensel step uniquely lifts every output through all powers
`p^a`. Set `w=2` and the inactive higher-binomial indices in their zero fibres;
their fixed contribution merely translates the requested residue.

## 3. Every prime above eight

Let

```text
A_k={C(t,k):t in F_p}
```

and let `A_2^reg` use only triangular inputs with `2t-1!=0`. The triangular
map is two-to-one away from its unique critical point, while a degree-`k`
polynomial has fibres of size at most `k`. Therefore

\[
|A_2^{\rm reg}|={p-1\over2},
\qquad |A_k|\ge\lceil p/k\rceil\quad(k=4,6,8).          \tag{10}
\]

Cauchy--Davenport makes the fourfold sumset all of `F_p` whenever

\[
{p-1\over2}+\lceil p/4\rceil+\lceil p/6\rceil
+\lceil p/8\rceil-3\ge p.                              \tag{11}
\]

A check of `p mod 24` proves (11) except at

```text
p=11,13,17,19,23,29,31,47.                             (12)
```

For these eight finite rows, exact successive sumset sizes are

```text
11:  5,10,11,11       13:  6,13,13,13
17:  8,16,17,17       19:  9,19,19,19
23: 11,22,23,23       29: 14,28,29,29
31: 15,31,31,31       47: 23,46,47,47.                (13)
```

Thus every residue modulo `p` has a representation whose triangular
derivative `(2w-1)/2` is nonzero. For any prescribed `p`-adic target, hold
`x,y,z` fixed and Hensel-lift `w`. This proves surjectivity onto `Z_p` for
every `p>8`.

## 4. Coordinatewise CRT and positivity

Factor `m` into prime powers and choose a local solution at each one. For
each coordinate of degree `k`, record its residue modulo the exact period
`P_(p^a,k)` in (3). Periods attached to distinct primes are coprime, so CRT
combines them coordinatewise. Equation (4) then preserves every local
binomial value. Finally add a common multiple of the coordinate period to
each index until it is at least two. This proves (2).

## 5. Failure boundary and connections

The theorem proves all of the following:

- every target has a solution in `Z_p^4` for every prime `p`;
- every target is represented modulo every positive integer;
- no fixed modulus, nor any finite conjunction of residue exclusions, can
  explain THM-4026's integer hole; and
- modular support is complete even though global integer support is not.

It does **not** prove an integer representation, a pointwise lower bound, or
compatibility with an Archimedean height box. A bounded set of residue masks
may still cover all candidates after the height bounds are imposed, exactly
as in THM-4026. The source-to-target map is reduction modulo the periods; it
preserves necessary congruence compatibility and destroys size. The missing
sidecar is a common bounded lift, supplied there by the exact square test.

Equations (3)--(5) connect to THM-4000/4010's Hasse/Newton observer lane, but
the ambient lattice here is `Int(Z)` and their ordinary-polynomial Smith forms
do not transfer verbatim. The CRT combination is a literal instance of the
THM-2059 fibre-product mechanism. It reinforces THM-2043/2050's global-clock
guardrail without supplying any LRC reduction.

## 6. Reproduction

```bash
python3 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
python3 -O 04-computation/sun_two_four_six_eight_counterexample_independent_audit_thm4026.py
```

The normal and optimized frozen streams are byte-identical. The proof report,
script/output hashes, and finite hostile rows are recorded in the frontmatter.
