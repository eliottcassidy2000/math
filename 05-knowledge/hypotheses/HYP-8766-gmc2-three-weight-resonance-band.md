# HYP-8766 — close the finite three-weight resonance band

**Status:** OPEN, reduced by THM-2017; one boundary model closed by THM-2014.

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
   that identity out from its first coefficients.
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
- HYP-8765 supplies exact finite-support evidence for a multilevel
  radial-channel tower and two counterexamples to naive low cutoffs.

The next decisive calculation is the explicit first two `1/m` coefficients
at `lambda=+-r`, expressed as differential operators on (3)-(4), followed by
a common-zero resultant for small `(p0,q0)`.

## Relation to HYP-8770

HYP-8770's symmetric-top Watson problem is the many-shell analogue of the
same entropy-versus-factorial obstruction. THM-2017 settles the separated
three-channel slopes and HYP-8766 isolates their exact resonance band;
HYP-8770 asks for uniform noncancellation when many near-top shell levels
populate the Watson boundary layer. Neither target subsumes the other:
HYP-8766 supplies explicit hyper-Bessel/transseries models for a controlled
slice, while HYP-8770 supplies the global cross-shell geometry that a full
NC2 proof must eventually absorb.
