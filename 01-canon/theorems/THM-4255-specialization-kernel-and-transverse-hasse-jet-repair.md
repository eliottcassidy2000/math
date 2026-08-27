---
id: THM-4255
title: "Specialization kernel, restricted short jets, and transverse Hasse repair"
status: >
  PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED. Graph specialization
  u -> f has kernel (u-f), and the inverse image of the order-n tail is
  (u-f)+(f^n). A restricted source admits the coefficientwise order inference
  exactly when its induced short-jet map modulo its own f^n tail is injective.
  Explicit-f Cartier channels do not descend coefficientwise through the graph
  quotient; the full transverse Hasse tower restores the lost normal data.
  This diagnoses, but does not repair or refute the conclusions of, the
  external p-adic-zeta density manuscript audited in the linked intake.
source: codex-padic-zeta-density-cartier-20260826
audit: >
  PASS. The primary exact companion verifies triangular truncations through
  degree nine, 54 graph-plus-tail cells, five Cartier controls, 14 transverse
  Hasse controls, and two hostile Cartier windows. A dependency-free referee
  exhausts all 512 elements of a 3-by-3 monomial box over F_2, independently
  checks graph division, short-jet ranks, both windows, and restricted
  ell-adic grade loss. Normal, optimized, and distinct-hash-seed outputs agree.
depends_on: []
related:
  - THM-4091-integral-coordinate-change-lcm-depth-boundary
  - THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification
script: 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
output: 05-knowledge/results/padic_specialization_kernel_cartier_firewall_thm4255.out
independent_script: 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
independent_output: 05-knowledge/results/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.out
script_sha256: f2f6b5edf6cdc61d73b0c535aa1995b366ddc0734a0587a24dbf9805f54e2f7b
output_sha256: 99c664cebf53cb4e10a068a11ef35f91f5dd06ead69f48890a16354cb43a59bd
independent_script_sha256: e4fa806af39fead477eb51e88864c9af667d3afc1c16b3411e6071f597c811c5
independent_output_sha256: cade91738365f0b79901584981677821eeb253763299ccaf3a28c8cf01601169
hash_basis: raw working-tree bytes (LF)
---

# THM-4255 -- specialization kernel, restricted short jets, and transverse Hasse repair

**PROVED + VERIFIED-EXACT + INDEPENDENTLY AUDITED.** Specialization along a
graph forgets a normal direction. Vanishing after specialization therefore
means graph-ideal membership, not coefficientwise vanishing in the source.
The exact replacement is a restricted short-jet injectivity condition; the
lost normal coordinate is recovered by transverse Hasse jets.

## 1. The graph kernel and every order tail

Let `R` be a commutative ring and let

```text
S=R[u][[f]],                  ev:S -> R[[f]],
ev(u)=f,                      ev(f)=f.                (1)
```

Here `R[u][[f]]` means an `f`-adic series whose coefficient at each fixed
`f`-degree is a polynomial in `u`. For every `F in S` and every `n>=0`,

```text
ker(ev)=(u-f)S,                                      (2)
ev^(-1)(f^n R[[f]])=(u-f)S+f^nS.                    (3)
```

To prove this, write `F=sum_i a_i(u)f^i`. Polynomial division gives

```text
a_i(u)=a_i(f)+(u-f)b_i(u,f).                         (4)
```

The series `Q=sum_i f^i b_i(u,f)` belongs to `S`: any fixed `f`-coefficient
receives contributions from only the finitely many indices below it. Hence

```text
F=(u-f)Q+ev(F).                                      (5)
```

Equation (2) follows immediately. If `ev(F)=f^n h(f)`, equation (5) gives
`F in (u-f)S+f^nS`; the reverse inclusion is immediate, proving (3).
Multiplication by `u-f` is injective even when `R` has zero divisors: at the
least nonzero `f`-coefficient of a putative annihilator, multiplication by the
indeterminate `u` cannot vanish. Thus (5) is the unique graph-normal
decomposition.

In particular,

```text
ev(F) in f^n R[[f]]  does not imply  F in f^nS.       (6)
```

The minimal witness is `F=u-f`, whose image is zero to every order.

## 2. Exact restricted-source repair

Let `V` be any `R`-submodule of `S`. The induced short-jet map is

```text
J_(V,n): V/(V intersect f^nS) -> R[[f]]/(f^n),
         [v] |-> ev(v) mod f^n.                      (7)
```

Then the desired coefficientwise inference on `V`,

```text
for every v in V,
ev(v) in f^nR[[f]]  implies  v in f^nS,              (8)
```

holds if and only if `J_(V,n)` is injective. Its exact obstruction is

```text
ker J_(V,n)
 = [V intersect ((u-f)S+f^nS)]/[V intersect f^nS].   (9)
```

This quotient is essential. Plain injectivity of `V -> R[[f]]/(f^n)` is
equivalent only in the special case `V intersect f^nS=0`. Consequently a
dimension, degree, or global-section bound on an unspecialized coefficient
space does not by itself prove (8); one must calculate (9) for the actual
restricted source.

## 3. Why coefficientwise Cartier does not descend

For an integer `ell>=2`, let `C_r^f` extract the terms whose **explicit source
`f`-exponent** is congruent to `r mod ell`, treating `u` as a coefficient.
Then

```text
C_0^f(u-f)=u,                  C_1^f(u-f)=-1,        (10)
ev(C_0^f(u-f))=f,              ev(C_1^f(u-f))=-1,   (11)
C_r(ev(u-f))=0 for every r.                           (12)
```

Thus these pre-specialization channels are not well-defined on the quotient
by `(u-f)`. Cartier extraction performed **after** the quotient is perfectly
well-defined; it is not the same channel decomposition. Indeed the monomial
`u^a f^i` lies in source channel `i mod ell`, while its specialization lies
in scalar channel `i+a mod ell`. The discarded `u`-degree is precisely the
missing channel sidecar.

## 4. Transverse Hasse jets recover the lost direction

Set `v=u-f`. Substitution gives an isomorphism between `S` and the row-finite
series in `f` with polynomial `v`-coefficients. Every `F` has a unique
expansion

```text
F=sum_(j>=0) v^j F_j(f),
F_j(f)=ev(D_u^[j] F),                               (13)
```

where `D_u^[j]` is the `j`-th Hasse derivative in `u`. Hence the full tower

```text
(ev(F), ev(D_u^[1]F), ev(D_u^[2]F), ...)            (14)
```

is injective and reconstructs `F`. More sharply, if

```text
F=(u-f)^m G,                 ev(G) != 0,             (15)
```

then the jets of orders below `m` vanish and

```text
ev(D_u^[m]F)=ev(G) != 0.                             (16)
```

For a source of bounded `u`-degree, only finitely many transverse jets are
needed. For an unrestricted source, no fixed finite tower is automatically
complete.

## 5. A separate arithmetic strictness gate

Graph-normal injectivity does not imply strictness for another filtration.
Over `Z_(ell)`, under the moving section

```text
u |-> 1+ell f^N,                                    (17)
```

the primitive source element `u-1` maps to `ell f^N`. Thus a restricted
coefficient subspace can lose its `ell`-adic associated grade. The full graph
quotient is still a split strict surjection; the failure is the inference on
the restricted pre-specialization subspace. Any arithmetic application needs
both the short-jet condition (9) and strictness for its coefficient lattice.

## 6. Consequence for the external density manuscript

The supplied manuscript's Propositions 6.2 and 12.3 infer unspecialized
coefficientwise order from scalar vanishing after a torsor specialization and
then invoke coefficientwise Cartier/no-backflow. Equations (2)--(12) show that
this inference is invalid without an additional restricted-source theorem.

Two abstract hostile controls fit the syntactic window template:

```text
Prop. 6.2 template:  D=10, C=0, ell=23, N=15,
                     ell^(-1)(u-1)z,
                     u=1+f^N, z=f^ell,
                     scalar pivot n=38;
                     n=j+ell*r has no r>=1, 0<=j<=10.

Prop. 12.3 template: D=12, C=1, ell=29, N=17,
                     ell^(-1)(u-1)z,
                     u=1+f^N, z=f^ell,
                     scalar pivot n=46;
                     n=j+ell*r has no r>=1, -1<=j<=13.  (18)
```

These are countermodels to the stated algebraic inference, not counterexamples
to the manuscript's modular geometry or density conclusions. The missing
repair could be supplied by either:

1. proving (9) vanishes for every actual finite digit space, together with
   `ell`-adic associated-grade strictness; or
2. working after actual specialization and proving that each nonzero scalar
   digit has order at most `D+C<ell` before applying scalar no-backflow.

Until one of those theorems is proved, the narrow Cartier windows and results
depending on them are **OPEN / unsupported by that manuscript**, not refuted.

## 7. Scoped transfers

- **Older p-adic-zeta irrationality draft:** its 22 claims remain
  `AUTHOR-CLAIMED / UNREFEREED`. The new decisive audit is to compute (9) for
  its actual `(p,s)=(5,5)` digit module and test Cartier preservation there.
  THM-4089's formula optimization and THM-4091's LCM transport remain safe.
- **Planar Jacobian conjecture:** specialization to a wall such as `W=0`
  should retain the `W`-adic strict transform or conormal Hasse jet. For each
  residual relation `R(W,t)`, compute `nu_W(R)` and
  `W^(-nu_W(R))R|_(W=0)`, then compare `gr_W(ker)` with `ker(gr_W)`. This is a
  lawful entry diagnostic, not an exclusion of any remaining incidence.
- **LRC(14):** local Hasse/Frobenius packets may be complete yet owner-height
  blind. The graph-kernel result reinforces the existing requirement to retain
  the height/owner sidecar; it supplies no loneliness bound.
- **Rule 30:** THM-3488's retained parity/Hasse channel is a positive control
  for (13)--(16). It supplies no prize resolution.

The percolation preprint at arXiv:2608.23661 is a source mismatch and has no
typed map into this theorem: the shared word “density” preserves no object,
predicate, or operation.

Reproduce the exact companions with

```bash
python -B 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
python -B -O 04-computation/padic_specialization_kernel_cartier_firewall_thm4255.py
python -B 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
python -B -O 04-computation/padic_specialization_kernel_cartier_firewall_independent_audit_thm4255.py
```
