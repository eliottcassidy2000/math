# `p-adic-zeta-density` specialization-kernel audit (2026-08-26)

**Status: PREPRINT CLAIM / UNDER SPECIALIST AUDIT. The proposed `u-f`
objection is valid algebra but does not hit the maps written in Propositions
6.2 or 12.3. Genus-zero Cartier receives a conditional type pass; the
general-curve coefficient-module diagram remains audit-open. No density or
irrationality theorem is promoted to proved canon here.**

## Frozen source

- Christopher D. Long,
  [*Multiweight arithmetic holonomy and density-one irrationality of p-adic
  zeta values*](https://github.com/octonion/p-adic-zeta-density), frozen at
  [commit `4c87bcdf4d7d62d0f1981f16e228901f02cd9f57`](https://github.com/octonion/p-adic-zeta-density/commit/4c87bcdf4d7d62d0f1981f16e228901f02cd9f57),
  commit time 2026-08-26 19:22:27 -0400.
- The source first appears in commit `234434b8`; the later commits through
  `4c87bcdf` change README equation typography, not the mathematical TeX.
- Frozen SHA-256 at `4c87bcdf`:
  - TeX: `abc8da8c5209d20d19c1ac6efc66e9d3f38d4d18da8f2bad9a370c96d651e48e`;
  - PDF: `eed1b53981e12967ee8db0781c7bc784fec2a443b27a57543fe549eae204d62a`;
  - supplied certificate: `f9877e4c64433d30736f4bf6a2e8fcae9689643cf41ac45ecbf5908289ea3c56`;
  - supplied output: `065eab6919eca798365560c3b7ae54802c4bf029221af8282d5597754ebd6bf6`.

The paper claims quantitative density-one irrationality at every fixed prime,
genus-zero and positive-genus finite-packet counts, and the large-prime
Cartier windows used in their slope inequalities. These remain author claims.

## The objection and its exact algebraic content

For

```text
ev:Q[u][[f]] -> Q[[f]],        u |-> f,
```

one has

```text
ker(ev)=(u-f),
ev(u-f+f^N)=f^N although u-f+f^N is not in f^N Q[u][[f]].
```

Thus vanishing or high order after a chosen-section specialization does not
imply coefficientwise vanishing before specialization. The exact theorem,
normal Hasse-jet repairs, universal-slope repair, and finite controls are
[THM-4255](../../01-canon/theorems/THM-4255-specialization-kernel-and-transverse-hasse-jet-repair.md).

This objection is decisive only if the manuscript uses the map `u->f` or an
equivalent `q`-dependent torsor contraction. A line-by-line search found no
such map. The distinct substitutions `z=f^ell` and `z=q^ell` are present and
must be audited through the manuscripts' block-separation hypotheses.

## Proposition 6.2: universal chart, not chosen section

The current source says at lines 1042--1049 that the equality is **on one
fixed flat torsor chart** and writes

```text
F_lambda(f)=H_lambda(f,u)+sum Q_(s,e,lambda)(f,u)R_(s,e)(f),
```

with `u` left as free fibre coordinates. The products and pole symbols at
lines 1094--1121 still retain `u`. The predecessor draft makes the intended
typing explicit at
[lines 1528--1546](https://github.com/octonion/p-adic-zeta-irrationality/blob/b46a1770901551961710e155d775aae7c5ea39e7/hybrid_arithmetic_holonomy.tex#L1528-L1546):
`F_lambda` is pulled back to a flat integral torsor chart, is independent of
the torsor point, and the fibre variables are coefficients receiving no
Frobenius action.

Consequently, if a lower coefficient of the intrinsic scalar contraction is
zero on the base, its pullback is the zero element of the torsor coefficient
algebra. No converse to specialization is used. The nonzero reduced leading
scalar remains nonzero in any nonzero unital `F_ell`-algebra. Neither step
requires faithful flatness; that would matter only for descending an equality
in the reverse direction.

The current density draft is terser than its predecessor. It should name the
integral chart algebra `A_ell` and explicitly state a separated completed
presentation `A_ell[[f]]` with fibre coordinates in `A_ell`. Subject to that
standard torsor-trivialization statement, Proposition 6.2 passes this
specific audit. The `u-f` hostile is not a counterexample to it.

## Proposition 12.3: intended universal module, diagram still owed

Lines 1911--1941 put multiplier coefficients into a fixed global coefficient
bundle and prove short-jet injectivity for its reduced sections. Lemma 12.2 is
formulated for `q`-germs in that coefficient space. Lines 2005--2014 then
digitize **coefficientwise in the finite free global coefficient lattice**
and assert, before scalarization,

```text
Pi_b(q)=Gamma_b(q,q^ell).
```

Again, there is no displayed torsor specialization `u=u(q)`. There is an
explicit block specialization `z=q^ell`, but Lemma 12.2 supplies its intended
no-cancellation gate: every nonzero coefficient `w_r(q)` has order at most
`D+C<ell`, so the first nonzero block begins strictly before the next one.
The formal hostile `z-q^ell` violates that order bound, as does
`z(z-q^ell)` in its first nonzero `z` coefficient. To make the argument
publication-grade, the paper should name the module containing both sides,
define vector-valued `q`-order, prove that each digit coefficient inherits
the stated bound, and embed the invariant scalar pole symbol before any
`q`-dependent contraction.

The binary hostile is exact. If the actual map were instead

```text
sigma:k[[q]]^2 -> k[[q]],       sigma(a,b)=a+q b,
```

then the nonzero vector `(q,-1)` would vanish after contraction despite both
components having small order. Under that hidden-section branch Proposition
12.3 would fail. The written source asserts the universal branch, so the
proper status is **audit-open conditional pass**, not refuted.

## Conditional dependency fork

If a future source audit proves that either identity exists only after a
chosen-section specialization, then the fallout is large:

```text
Prop 6.2 -> Props 6.3/6.6 -> Thms 1.1/1.2
         -> genus-zero finite counts and density;

Prop 12.3 -> Prop 12.4 -> Prop 13.5 / Thm 13.6
          -> exact/robust slopes -> all-prime density and
             positive-genus finite counts.
```

The pole-depth bound `1<=a<=s_*`, divided-Frobenius identities, zero
estimates, standalone Fitting-length cap, capacity calculations, and
conditional numerical inequalities would survive this particular objection.
This fork is recorded for audit planning only; the antecedent has not been
established and no downstream theorem is retracted here.

## Reproducibility findings independent of the objection

The repository tracks only the TeX/PDF, README/LICENSE, and one script/output
pair. Seven files cited by the manuscript are absent:

```text
hybrid_rational_interval_certificate.py
hybrid_rational_interval_certificate_output.txt
hybrid_rational_interval_certificate_sha256.txt
multiweight_all_primes_interval_certificate.py
multiweight_all_primes_interval_certificate_output.txt
multiweight_all_primes_interval_certificate_sha256.txt
allprime_capacity_finite_packet_certificate_sha256.txt
```

The supplied positive-genus script runs under normal and optimized Python and
both outputs byte-match the tracked output. Its actual hashes are the frozen
values above. The manuscript instead prints `52e404ad...` and `3e2a4104...`,
so its embedded manifest is stale. This is a real reproducibility defect but
does not change the conditional arithmetic computed by the script.

There are no tests, CI workflows, tags, or issue/PR audit trail in the frozen
source. As in the earlier
[`p-adic-zeta-irrationality` audit](p-adic-zeta-irrationality-source-audit-20260825.md),
a numerical `PASS` checks final formulas, not the new geometric and adelic
bridges.

## Cheapest decisive source repairs

1. Add the completed universal torsor-chart diagram for Proposition 6.2 and
   state explicitly that no section `u=u(f)` is chosen.
2. Add the finite-free coefficient-module diagram for Proposition 12.3,
   including the scalar inclusion and vector-valued order.
3. Supply one worked torsor/pole-grade example whose universal coefficient
   identity can be checked before the substitutions `z=f^ell` or `z=q^ell`.
4. Restore the seven cited certificate artifacts and regenerate the printed
   manifest.

Until those tasks and the other gates in the prior source audit receive
specialist review, the density claims remain **PREPRINT CLAIM / UNDER
SPECIALIST AUDIT**.
