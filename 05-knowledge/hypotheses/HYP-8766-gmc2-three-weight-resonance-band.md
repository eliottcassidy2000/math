# HYP-8766 — close the finite three-weight resonance band

**Status:** OPEN, reduced by THM-2017 and THM-2018; exact low-degree and
proportional-central models are closed by THM-2014/2018.

**Owner:** codex-2026-07-21-gmc2-degree-gap.

**Target:** NC2, hence GMC(2), for

```text
P=Z^p a(s)+b(s)+Zbar^q c(s),      s=Z Zbar,
```

with nonzero `a,b,c in C[s]` and arbitrary positive charges `p,q`.

Put

```text
g=gcd(p,q), p0=p/g, q0=q/g, r=p0+q0,
h=s^(pq/g) a^q0 c^p0, d=deg b, e=deg h, lambda=e-rd.
```

THM-2017 proves the target whenever `|lambda|>=r+1`, and at
`|lambda|=r` unless one explicit hyper-Bessel value vanishes. Therefore the
entire unresolved degree geometry consists of

```text
lambda in {-r,-r+1,...,r-1,r},                                     (1)
```

with the two endpoints already reduced to discrete exceptional
leading-coefficient loci. THM-2017 removes even those exceptions when
`p0=q0=1` and `b,h` are monomials: the universal `1/m` correction is a
nonzero Bessel derivative at every zero of the leading term.

## Conjecture

For every fixed `(p,q)` and every offset in (1), the exact channel identity

```text
E[P^m]=sum_k m!/((q0k)!(p0k)!(m-rk)!) L(h^k b^(m-rk))               (2)
```

cannot vanish for all positive `m` unless one of `a,b,c` vanishes; after the
degenerate two-weight/neutral cases are resolved, the moment-null members are
exactly the strict one-sided ones.

## Why this is a finite asymptotic problem, not a domination problem

At the boundary `rd-e=r`, each fixed `k` has a finite nonzero scale relative
to the neutral endpoint. If `alpha,beta` are the leading coefficients of
`h,b`, then formally (and rigorously under the same mixed-factorial bounds as
THM-2017)

```text
E[P^m]/L(b^m) -> Phi_(p0,q0)(xi),
Phi_(p0,q0)(xi)=sum_{k>=0} xi^k/((q0k)!(p0k)!),
xi=alpha/(beta^r d^r).                                                (3)
```

THM-2017 proves (3), not merely formally. Thus boundary nonvanishing is not
literal endpoint dominance: the endpoint is resummed into a hyper-Bessel
function, which can have complex zeros. At `e-rd=r`, on `m=rn`, removal of
`j` primitive returns analogously and rigorously
produces

```text
sum_{j>=0} eta^j/(rj)!,
eta=(q0^q0 p0^p0/e^r)(beta^r/alpha).                                 (4)
```

For `0<|lambda|<r`, the active boundary layer has
`k or j ~ m^((r-|lambda|)/r)`, so a large-argument generalized Bessel saddle
replaces (3)-(4). At `lambda=0`, channels proportional to `m` share the same
factorial degree and a full entropy saddle is unavoidable.

This explains the sharpness of THM-2017: the first charged channel has
multinomial gain `m^r`, while a degree deficit `Delta` supplies only
`m^Delta` of factorial loss. Asking endpoint domination to cross
`Delta<=r` contradicts the scale of (2).

## Proposed closure program

1. **Boundary offsets `lambda=+-r`.** Extend the proved limits (3)-(4) to a
   full `1/m` expansion. If the leading hyper-Bessel value is zero, show the next
   coefficients form differential operators in `theta=xi d/dxi` with no
   common nonzero zero. THM-2017 carries this out for monomial `b,h` when
   `(p0,q0)=(1,1)`; THM-2014's linear exact EGF is the nonmonomial resonant
   model of the same principle.
2. **Intermediate offsets.** Apply saddle analysis to the hyper-Bessel
   boundary layer. Separate the finitely many exponential phases by a
   subsequence in `m`; all-moment vanishing would require every saddle
   amplitude to vanish.
3. **Central offset `lambda=0`.** Use the proportional-channel entropy
   functional, not a top-degree order. Prove that some arithmetic subsequence
   has a unique maximal real part after a harmless rotation/scaling of `P`,
   or convert equality of saddles into an exact holonomic identity and rule
   that identity out from its first coefficients. THM-2018 carries out the
   second option on the arbitrary-charge hypersurface `h=kappa*b^r`: the
   channel sum factors as `H_m(kappa)L(b^m)`, and root-of-unity symmetry proves
   that the hyper-Bessel coefficient sequence `H_m` is not eventually zero.
4. **Coefficient descent.** Once leading coefficients are eliminated, repeat
   on the next Newton layer. This is the analytic analogue of HYP-8765's
   radial-channel resultant tower.

## Assumption challenge

- First-return atoms are not independent equations (MISTAKE-211).
- Equal factorial degree does not imply equal asymptotic contribution;
  multinomial entropy selects a boundary layer or interior saddle.
- A hyper-Bessel leading constant can vanish at isolated complex parameters,
  so proving only its nontriviality as a function is insufficient. The proof
  needs a derivative/transseries tower with empty common zero locus.
- The useful quotient vertices are radial channels. It preserves factorial
  slope but destroys phase, so Tournament Analysis locates regimes and cannot
  certify noncancellation.

## Evidence

- THM-2014 closes the full constant-endpoint `(p,q)=(1,1)` model, including
  the resonant degrees zero and one by exact EGFs.
- THM-2017 closes every offset outside (1) and verifies (2) independently
  against direct Wick expansion.
- THM-2018 closes the complete arbitrary-complex slice `p=q=1`,
  `deg b<=1`, `deg(a*c)<=1`. Its exact Catalan resummation is
  `F(t)=exp(A(t))/sqrt((1-b1*t)^2-4*alpha*t^2)` with `A` algebraic; nullity
  would make the exponential of a nonconstant algebraic germ algebraic. It
  also closes `h=kappa*b^r` for arbitrary charges and arbitrary radial degree.
- HYP-8765 supplies exact finite-support evidence for a multilevel
  radial-channel tower and two counterexamples to naive low cutoffs.

The boundary `1/m` calculation is now complete in the symmetric monomial
model (THM-2017), while THM-2018 supplies exact non-asymptotic closures inside
the central band. The remaining decisive work is a phase-sensitive saddle or
holonomic separation theorem away from `h=kappa*b^r`, followed by coefficient
descent; on the general sharp boundary it must also remove the discrete
hyper-Bessel zero loci for nonmonomial `b,h`.

## Exact resonance closures from THM-2018

THM-2018 turns two of the proposed mechanisms into theorems.

1. For `p=q=1`, `b=b0+b1*s`, and `a*c=delta+alpha*s`, it evaluates the
   ordinary moment EGF exactly. The Catalan identity
   `[x^j]C(x)^n/sqrt(1-4x)=binom(n+2j,j)` resums every return channel into an
   algebraic-exponential germ. If the EGF were one, an algebraic germ and its
   exponential would both be algebraic, which is impossible unless the germ
   is constant. This forces `b=0` and `a*c=0` with no sign or reality
   hypothesis.
2. For all `p,q`, if `h=kappa*b^r`, every return channel has the identical
   radial shadow `b^m`. The remaining scalar sequence has EGF
   `exp(t)Phi_(p0,q0)(kappa*t^r)`. If it were eventually zero, this EGF would
   be a polynomial, but its transformation under a nontrivial `r`-th root of
   unity would make a nonconstant exponential rational. EMP then supplies a
   large index at which both factors are nonzero.

The second argument is the exact central-band observability principle that
the proposed entropy program was seeking: root-of-unity symmetry restores the
phase coordinate destroyed by the tied radial-channel quotient. What remains
is to deform this separation away from exact proportionality.

## Incoming Sheffer/positivity connection: HYP-8769

HYP-8769 identifies the `p=q=1`, constant-negative-endpoint specialization
of (2) with the earlier S64 Sheffer channel sequence. Under that dictionary,
the differential/transseries tower proposed above is the Sheffer analogue of
the S62 Hermite no-common-root mechanism. This is a genuine unification of
the algebraic and asymptotic attacks, but not yet a proof: the required
no-common-nonzero-zero statement is precisely what remains open beyond the
base Hermite case.

It also isolates the phase obstruction cleanly. If every coefficient is real
and nonnegative, each Wick/channel contribution is nonnegative and at least
one is positive, so moment nullity is impossible. Any counterexample in the
resonance band must therefore exploit mixed signs or genuinely complex phase.
A successful Hankel/SOS/Bargmann reformulation would have to encode those
signed channel sums themselves; positivity of the raw norm `E[|P|^2]` is
orthogonal to nullity and cannot close the band by itself. This makes
phase-sensitive Sheffer resultants, rather than another top-term comparison,
the concrete algebraic companion to the saddle program.

## Relation to HYP-8770

HYP-8770's symmetric-top Watson problem is the many-shell analogue of the
same entropy-versus-factorial obstruction. THM-2017 settles the separated
three-channel slopes and HYP-8766 isolates their exact resonance band;
HYP-8770 asks for uniform noncancellation when many near-top shell levels
populate the Watson boundary layer. Neither target subsumes the other:
HYP-8766 supplies explicit hyper-Bessel/transseries models for a controlled
slice, while HYP-8770 supplies the global cross-shell geometry that a full
NC2 proof must eventually absorb.
