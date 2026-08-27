---
id: THM-4255
title: "Specialization kernels, restricted sources, and transverse Hasse/Kronecker repair"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Evaluation on one formal
  arc has principal kernel and need not reflect coefficientwise f-order. The
  exact restricted-source obstruction is a short-jet kernel; finite normal
  Hasse jets, enough sections, a universal slope, and finite-box Kronecker
  substitution give sharp repairs. EXTERNAL APPLICATION CORRECTED: the
  p-adic-zeta density draft does not display the torsor specialization needed
  by the proposed u-f counterexample. Its Propositions 6.2/12.3 remain
  PREPRINT CLAIMS / UNDER SPECIALIST AUDIT, but are not refuted on that ground.
source: codex-padic-zeta-density-cartier-20260826
audit: >
  PASS. Two primary exact companions and a dependency-free finite-field
  referee verify graph kernels, restricted short jets, Hasse recovery,
  Cartier hostiles, interpolation, universal-slope boxes, and sharp finite-box
  Kronecker boundaries. Normal and optimized runs match their frozen outputs.
  MISTAKE-527 records the corrected source application.
depends_on: []
related:
  - THM-3476-rule30-depth-four-transverse-jet-barrier-and-slack-pascal-atlas
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
  - THM-4091-integral-coordinate-change-lcm-depth-boundary
  - THM-4120-jc23-extremal-target-degree-twenty-one-response
script: 04-computation/specialization_kernel_hasse_jet_repair_thm4255.py
output: 05-knowledge/results/specialization_kernel_hasse_jet_repair_thm4255.out
script_sha256: 8f60050e6f52afd469622d3ee53bff2c65332871951b5c0c0a4d8b96881b5979
output_sha256: fa0d507204427ee929e79f3144670350936cf7056d4f46450b8bc80b7306f32a
companion_script: 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
companion_output: 05-knowledge/results/padic_specialization_kernel_cartier_firewall_thm4255.out
companion_script_sha256: f2f6b5edf6cdc61d73b0c535aa1995b366ddc0734a0587a24dbf9805f54e2f7b
companion_output_sha256: 99c664cebf53cb4e10a068a11ef35f91f5dd06ead69f48890a16354cb43a59bd
independent_script: 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
independent_output: 05-knowledge/results/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.out
independent_script_sha256: e4fa806af39fead477eb51e88864c9af667d3afc1c16b3411e6071f597c811c5
independent_output_sha256: cade91738365f0b79901584981677821eeb253763299ccaf3a28c8cf01601169
hash_basis: raw LF bytes
external_source: https://github.com/octonion/p-adic-zeta-density/commit/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57
---

# THM-4255 -- specialization kernels, restricted sources, and transverse repair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. External density claims
remain PREPRINT CLAIMS / UNDER SPECIALIST AUDIT.**

The algebraic warning

```text
u-f != 0 in Q[u][[f]], but (u-f)|_(u=f)=0                 (1)
```

is correct. It diagnoses a chosen-section calculation only after that
section map has actually been exhibited. A universal identity on a torsor
chart keeps `u` free. These are different typed statements. This theorem
gives the exact kernel calculus and repairs, then audits which map occurs in
the external density draft.

## 1. One formal arc has a principal kernel

Let `A` be a nonzero unital commutative ring and

```text
R=A[u][[f]]={sum_(n>=0) f^n p_n(u): p_n in A[u]}.       (2)
```

For any `g in A[[f]]`, substitution defines a continuous map

```text
ev_g:R -> A[[f]],              P(f,u) |-> P(f,g(f)).    (3)
```

### Theorem 1.1 (kernel and every order tail)

For every `N>=0`,

```text
ker(ev_g)=(u-g)R,                                        (4)
ev_g^(-1)(f^N A[[f]])=f^N R+(u-g)R.                     (5)
```

Hence `R/(u-g)R ~= A[[f]]`. The quotient filtration is strict, but one
section is not order-reflecting on `R`.

### Proof

Write `P=sum_n f^n p_n(u)`. Polynomial division gives

```text
p_n(u)-p_n(g)=(u-g)q_n(u,g).                            (6)
```

After powers of `f` are regrouped, `sum_n f^n q_n(u,g)` lies in `R`: only
the finitely many `n<=M` can affect the coefficient of `f^M`. Thus

```text
P(f,u)-P(f,g) in (u-g)R.                                (7)
```

Equation `(4)` follows by setting `ev_g(P)=0`; `(5)` follows by subtracting
the scalar tail `f^N h`. Both reverse inclusions are immediate. QED.

Multiplication by `u-g` is injective even if `A` has zero divisors. At the
least nonzero `f`-coefficient of a putative annihilator, one multiplies a
nonzero polynomial by the monic polynomial `u-g(0)`.

### Corollary 1.2 (unbounded order inflation)

For every `N>=1`,

```text
P_N=u-g+f^N                                             (8)
```

has coefficientwise `f`-order zero and `ev_g(P_N)=f^N`. Thus no fixed bound
controls the order gained under one specialization, even in `u`-degree one.
If

```text
P=f^m p_m(u)+O(f^(m+1)),        p_m!=0,                 (9)
```

then the leading coefficient after evaluation is `p_m(g(0))`. A one-section
argument must retain this gate or control the kernel summand in `(5)`.

## 2. The exact restricted-source criterion

Let `V` be any `A`-submodule of `R`. Define

```text
S_(V,N):V/(V intersect f^N R) -> A[[f]]/(f^N),
        [v] |-> ev_g(v) mod f^N.                       (10)
```

Then

```text
ev_g(v) in f^N A[[f]]  implies  v in f^N R for every v in V
```

if and only if `S_(V,N)` is injective. Its exact obstruction is

```text
ker S_(V,N)
 = [V intersect ((u-g)R+f^N R)]/[V intersect f^N R].   (11)
```

This is immediate from `(5)`. It is the sharp repair when geometry restricts
the admissible coefficient space. Dimension or degree bounds alone do not
prove injectivity; one must calculate `(11)` for the actual `V`.

## 3. Transverse Hasse jets recover the normal direction

Let `D_u^[j]` be the `j`th Hasse derivative,

```text
D_u^[j](u^i)=binom(i,j)u^(i-j),                         (12)
```

and define

```text
J_(g,r)(P)=(D_u^[j]P(f,g(f)))_(0<=j<=r).               (13)
```

### Theorem 3.1 (exact jet kernel)

For every `r>=0`, in every characteristic,

```text
ker J_(g,r)=(u-g)^(r+1)R.                              (14)
```

If `R_(<=d)` denotes the elements of uniformly bounded `u`-degree at most
`d`, then

```text
J_(g,d):R_(<=d) -> A[[f]]^(d+1)                        (15)
```

is an isomorphism, with inverse

```text
(b_0,...,b_d) |-> sum_(j=0)^d b_j(u-g)^j.              (16)
```

Consequently, for every `N`,

```text
P in f^N R_(<=d)
 iff D_u^[j]P(f,g) in f^N A[[f]] for all 0<=j<=d.      (17)
```

For `d>=1`, the depth is sharp: `J_(g,d-1)` misses `(u-g)^d`.

### Proof

Hasse--Taylor gives the locally finite identity

```text
P(f,u)=sum_(j>=0) D_u^[j]P(f,g)(u-g)^j.                (18)
```

Vanishing of its first `r+1` normal coefficients is exactly divisibility by
`(u-g)^(r+1)`. On `R_(<=d)` the sum ends at `d`, proving `(14)--(17)`. QED.

Ordinary derivatives are not a characteristic-safe replacement: in
characteristic `p`, the ordinary derivatives through order `p` of `u^p`
vanish at zero, whereas `D_u^[p](u^p)=1`.

## 4. Three alternative lossless observers

### Theorem 4.1 (unit-separated interpolation)

If `g_0,...,g_d in A[[f]]` have unit pairwise differences, then

```text
P |-> (P(f,g_i(f)))_(0<=i<=d)                           (19)
```

is an isomorphism on `R_(<=d)`. For `d>=1`, only `d` sections miss the
hostile `product_(i=0)^(d-1)(u-g_i)`.

When `A=K` is a field, `K[[f]]` is an `f`-adic DVR. If a Vandermonde
determinant has order `v`, values modulo `f^(N+v)` suffice by the adjugate
formula to recover coefficients modulo `f^N`. This is sufficient, not always
sharp; the largest Smith exponent is the exact worst coordinate loss.

### Theorem 4.2 (universal-slope pencil)

The map

```text
Phi:R -> A[lambda][[f]],             u |-> lambda f    (20)
```

is injective without a `u`-degree bound, because

```text
[f^N]Phi(P)=sum_(n+i=N)c_(n,i)lambda^i                 (21)
```

recovers every diagonal coefficient polynomial. Over a field, one
transcendental slope valued in a field extension also works. No finite list of algebraic slopes is
injective on unrestricted `R`. The power of `lambda` is essential sidecar
data: scalarizing it discards the bidegree again.

### Theorem 4.3 (finite-box Kronecker repair)

Let `d>=1`, `J>=0`, `M>=1`, and suppose

```text
P=sum c_(n,i)f^n u^i,       0<=i<=d, 0<=n<=J.          (22)
```

Then `u |-> f^M` is injective on this box if and only if `M>J`. Sufficiency
follows because the exponents `n+Mi` have unique base-`M` decomposition.
If `M<=J`, the nonzero polynomial `u-f^M` is an admissible hostile.

Normal jets, several sections, the universal pencil, and Kronecker separation
preserve different sidecars. The cheapest valid observer depends on the map
and support actually supplied by the geometry.

## 5. Cartier and a second arithmetic filtration

No target-side operators can resurrect the kernel of one section. For any
family `C_alpha:A[[f]]->T_alpha`,

```text
ker(ev_g) subset ker((C_alpha o ev_g)_alpha).           (23)
```

For the special graph `g=f`, source Cartier channels make the loss explicit.
If `C_r^f` selects source monomials by their explicit `f`-exponent modulo
`ell>=2`, treating `u` as a coefficient, then

```text
C_0^f(u-f)=u,             C_1^f(u-f)=-1,
C_r(ev_f(u-f))=0 for every r.                           (24)
```

Thus coefficientwise source channels do not descend through the graph
quotient. Cartier after evaluation is well-defined, but it is a different
grading: `u^a f^i` moves from source channel `i` to scalar channel `i+a`.

Graph-normal strictness also does not imply strictness for another
filtration. Over `Z_(ell)`, under `u |-> 1+ell f^N`, the primitive element
`u-1` maps to `ell f^N`. An arithmetic application therefore needs both the
normal criterion `(11)` and strictness for its coefficient lattice.

For `I=(u-g)`, the canonical carrier of the discarded direction is

```text
gr_I R ~= A[[f]][V].                                    (25)
```

## 6. Corrected audit of the p-adic-zeta density draft

The external source is commit
[`4c87bcdf`](https://github.com/octonion/p-adic-zeta-density/commit/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57).
Equations `(4)--(24)` would expose a flaw if the manuscript first evaluated a
torsor coordinate on an `f`- or `q`-dependent section. It does not display
such a map in Propositions 6.2 or 12.3. The earlier contrary repo conclusion
is **SUPERSEDED** by MISTAKE-527.

### Proposition 6.2

The source works
[`on one fixed flat torsor chart`](https://github.com/octonion/p-adic-zeta-density/blob/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57/multiweight_arithmetic_holonomy_all_primes.tex#L1040-L1049)
and writes

```text
F_lambda(f)=H_lambda(f,u)+sum Q_(s,e,lambda)(f,u)R_(s,e)(f), (26)
```

with `u` free. The predecessor manuscript says explicitly that the intrinsic
scalar `F_lambda` is
[`independent of the torsor point`](https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/hybrid_arithmetic_holonomy.tex#L1528-L1546)
after pullback and that fibre variables remain coefficients. A zero scalar
coefficient therefore pulls back to literal zero in the coefficient algebra;
it is not obtained from `u=f`.

The actual substitution is the separate block variable `z=f^ell`. The stated
bound `deg_f Gamma<ell` makes `j+ell*r` a unique block encoding. Hostiles such
as `z-f^ell` violate this anti-aliasing bound.

**Verdict:** the `u-f` objection does not refute Proposition 6.2. The proof
still owes an explicit completed-chart presentation such as `A_ell[[f]]` and
remains an external specialist-audit obligation.

### Proposition 12.3

Section 12 says digitization is
[`coefficientwise in the finite free global coefficient lattice`](https://github.com/octonion/p-adic-zeta-density/blob/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57/multiweight_arithmetic_holonomy_all_primes.tex#L1997-L2016).
Again, no chosen torsor section `u=u(q)` occurs. The displayed block
specialization `z=q^ell` is protected conditionally by Lemma 12.2: if

```text
Gamma(q,z)=sum_(0<=r<ell) w_r(q)z^r,
ord_q(w_r)<=D+C<ell for every nonzero coefficient,       (27)
```

then the first nonzero block ends before the next begins. The hostiles
`z-q^ell` and `z(z-q^ell)` fail `(27)`. The intended identity

```text
Pi_b(q)=Gamma_b(q,q^ell)                                (28)
```

is asserted module-valued; `(27)` is its anti-aliasing sidecar.

**Verdict:** Proposition 12.3 is a **CONDITIONAL / AUDIT-OPEN** pass on this
objection. The draft must name the common finite-free module, its vector
`q`-order, the scalar pole-symbol inclusion, and prove that all block digits
inherit `(27)`. A hidden contraction such as

```text
sigma(a,b)=a+qb                                         (29)
```

would kill the nonzero vector `(q,-1)` and revive the objection, but no such
map is stated. Downstream density conclusions are not retracted unless that
universal-module assertion fails.

The graph-firewall companions' pivots `38/46` remain valid **abstract hostile
controls for a chosen-section template**. They are not countermodels to the
two source proofs, because the necessary torsor specialization is absent.

## 7. Cross-frontier consequences

- **Planar Jacobian:** [THM-4120](THM-4120-jc23-extremal-target-degree-twenty-one-response.md)
  is a positive control. On its `DE` edge, restriction has kernel `(w-T)` and
  the proof retains re-entry orders `ord_z(w-T)=1` or `2`. A future edge
  restriction with `w`-degree cap `d` can be made lossless with the `d+1`
  Hasse jets in `(15)` or `d+1` lawful fibres in `(19)`. This gives an audit
  protocol, not a Keller-map exclusion.
- **LRC(14):** a complete local/Frobenius packet can still forget the global
  owner-height coordinate. The theorem reinforces the need for that sidecar
  but supplies no loneliness bound.
- **Rule 30:** THM-3488's parity/Hasse channel is a positive control for
  `(14)--(18)`; it supplies no prize resolution.
- **Cerf percolation preprint:** arXiv:2608.23661 is a source mismatch, not a
  p-adic-density paper. Its local/nonlocal obstruction split is only a method
  analogy and has no theorem map into `(2)`.

## 8. Exact finite controls

The general companion checks over `F_p[f]/(f^N)` for
`p in {2,3,5,7}`, `N=2,...,8`, and `u`-degrees `1,...,6`: `504` principal
kernels, `504` full Hasse-jet isomorphisms, `2,016` order hostiles, `91`
section families, `168` universal-slope boxes, `168` sharp Kronecker
boundaries, and four characteristic-`p` derivative hostiles.

The graph-firewall companion checks triangular truncations through degree
nine, `54` graph-plus-tail cells, five Cartier controls, `14` Hasse controls,
and two abstract hostile windows. Its clean-room referee exhausts all `512`
elements of a `3x3` monomial box over `F_2` and independently checks graph
division and restricted short-jet ranks. These computations verify the
algebraic template only; source applicability is settled by the typed audit
in Section 6.

Reproduce with

```bash
python3 -B 04-computation/specialization_kernel_hasse_jet_repair_thm4255.py
python3 -B -O 04-computation/specialization_kernel_hasse_jet_repair_thm4255.py
python3 -B 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
python3 -B -O 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
python3 -B 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
python3 -B -O 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
```
