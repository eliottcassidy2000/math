# Cubic Russell seeds: the residual is gauge, the density debt is not

**Status:** FINITE-EXACT SCOUT, 2026-08-23.  This note records an exact
negative classification and its quotient audit.  It is not a theorem file,
constructs no Keller map, and has no consequence for `JC(2)`, which remains
open.

## 1. Inheritance and concept board

The closest proved mechanism is THM-3868: a monic hidden-control relation in
the normal Russell ring restores the controls, after which the seed different
would have to be a unit.  The canonical hostile is the alternative quadratic
seed

```text
p=-4e^2,                     q=-4e^2u+e,                    (1)
```

whose symplectic pole chart makes `A,C` regular while `u,x` both have simple
poles.  Its boundary seed Jacobian is `16e^2v`, so it is not Keller.  The
corrected near miss is the earlier phrase “first genuine cubic residual”:
that residual is genuine only relative to a chosen quadratic section.  The
least-used sidecar is the determinant cocycle of the source change.

The active board was:

| Object | Predicate | Lost coordinate | Cheapest test |
|---|---|---|---|
| filtered normal seed jets | quotient-invariant residual | source determinant | compute every graded action matrix |
| affine cubic seed | seed Jacobian in `k*` | none in raw space | solve all six Jacobian buckets |
| exceptional quadratic pole chart | regular cubic correction | valuation/effectivity | inspect the first two Laurent coefficients |
| fixed-seed density ODE | polynomial algebraization | squareclass | resum the unique truncation |

The exact companion is
[`jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py`](../04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py).

## 2. The quotient boundary

In normalized coordinates put

```text
A_0=u^2-x,
C_0=u^3-u-(3/2)ux,
J_(x,u)(A_0,C_0)=1+(3/2)x.                                  (2)
```

The ideal

```text
I_0=x k[u][[x]]^2                                           (3)
```

is the full source-addition ideal vanishing on the arm `x=0`.  Fixing the
first normal row restricts it to

```text
I_1=x^2 k[u][[x]]^2.                                        (4)
```

At every associated-graded normal degree, the new source coefficients
`(xi_n,eta_n)` act on the target coefficients by

```text
M(u)=[[-1,       2u],
      [-(3/2)u, 3u^2-1]],                    det M=1.         (5)
```

Thus `M` is invertible over `k[u]`, not merely over `k(u)`.  At quadratic
order the unique section for target coefficients `(p,q)` is

```text
eta_2=(3/2)pu-q,                 xi_2=2u eta_2-p.             (6)
```

Relative to this section, the displayed cubic residual is

```text
rho=(r, s+(3/2)eta_2).                                      (7)
```

But a cubic source addition realizes every such residual because

```text
(xi_3,eta_3)=M^(-1)rho
 =((3u^2-1)rho_A-2u rho_C, (3/2)u rho_A-rho_C).              (8)
```

Consequently **no cubic residual descends to the full cosmetic quotient**.
After each lower-order choice, the fresh coefficients at the next normal
grade again enter through `(5)`, so the unrestricted formal quotient has no
finite-grade residual of this kind.  This is a formal statement; it does not
give a polynomial automorphism, a rational source change, or effectivity at
global divisors.

The quotient loss ledger is:

```text
raw seed jets  ->  quotient by I_1 source additions
preserved:         arm values and first normal packet
destroyed:         cubic coefficient pair and source determinant
needed sidecar:    J(source change), with its complete divisor
decisive hostile:  (x,u)->(x+x^3,u), which stays in one formal orbit but
                    multiplies density by (1+3x^2).
```

Hence constant Jacobian is not a predicate on the bare quotient.  It becomes
typed only after retaining the determinant cocycle or restricting the source
action to a proved determinant-one/effective subgroup.  The phrase “first
genuine residual” must therefore mean “first residual in the chosen
quadratic section,” not “first quotient invariant.”

## 3. Raw affine cubic classification

Now do not quotient.  Search the raw polynomial cell

```text
A=u^2-x+p(u)x^2+r(u)x^3,
C=u^3-u-(3/2)ux+q(u)x^2+s(u)x^3,                            (9)
```

with `p,q,r,s` affine in `u`.  Its Jacobian has six buckets `1,E_1,...,E_5`.
The first row is

```text
E_1=3/2+(6u^2-2)p-4uq.                                     (10)
```

Setting it to zero forces uniquely

```text
p=3/4,                         q=(9/8)u.                    (11)
```

The second row then forces uniquely

```text
r=-9/8,                        s=-(27/16)u.                 (12)
```

These coefficients do not give a Keller seed.  They give the exact
one-variable precomposition

```text
chi_3=x-(3/4)x^2+(9/8)x^3,
A=u^2-chi_3,
C=u^3-u-(3/2)u chi_3,                                      (13)
```

whose determinant is

```text
1+(135/16)x^3-(405/64)x^4+(729/128)x^5.                    (14)
```

Thus the raw affine constant-J cell is empty, with first failure exactly
`E_3=135/16`.

This remains empty when the chosen quadratic section is restored to its full
polynomial fibre.  The complete `E_1=0` fibre is

```text
p=3/4+2u eta,
q=(9/8)u+(3u^2-1)eta.                                      (15)
```

Writing affine relative cubic residuals `alpha,beta` gives

```text
E_2=eta'-27/8+3(3u^2-1)alpha-6u beta,                       (16)

E_3=-beta'+27/16+9u eta+(12u^2+4)eta^2
    -(9/2)alpha+(3/2)u alpha'.                              (17)
```

Equation `(16)` makes polynomial `eta` have degree at most four.  In `(17)`,
the highest coefficients of degrees `10,8,6,4,2` successively equal

```text
12 eta_4^2, 12 eta_3^2, 12 eta_2^2, 12 eta_1^2, 12 eta_0^2. (18)
```

Hence `eta=0`.  Equation `(16)` then returns `(12)`, while `(17)` becomes
`135/16`, a contradiction.  This is the section-independent affine-relative
classification.  It is separate from the quotient statement: the former is
a statement in raw seed space with the determinant retained; the latter says
the cubic coefficient pair itself is not an invariant.

## 4. The debt is the first omitted Catalan coefficient

The unique continuation `(13)` is the cubic truncation of the exact density
ODE

```text
(1+(3/2)chi)chi'=1,                    chi(0)=0.             (19)
```

Integration and the chosen origin give

```text
chi+(3/4)chi^2=x,
chi=(2/3)(sqrt(1+3x)-1)                                    (20)
   =x-(3/4)x^2+(9/8)x^3-(135/64)x^4+... .
```

The derivative of the first omitted term is `-(135/16)x^3`, exactly
cancelling the debt in `(14)`.  This explains the coefficient structurally,
not as an isolated failed solve.

The exact solution is outside the classified polynomial cubic section.  Its
radicand `1+3x` has a simple prime valuation, so it is not a square in the
rational function field.  This is precisely the normalized form of
THM-3846's canonical Catalan/Hensel algebraization wall.  Thus the cubic
search rediscovers the same boundary from the opposite direction:

```text
finite polynomial truncation  -> nonzero next density bucket;
exact formal resummation       -> constant density but nonsquare algebraic sheet.
```

## 5. Recovery and effectivity separate the two near misses

The raw near-candidate `(13)` retains the original monic hidden-`u` equation

```text
u^3+(2-3A)u+2C=0.                                           (21)
```

After recovering `u`, one gets `chi_3=u^2-A`, and `x` obeys the monic cubic

```text
x^3-(2/3)x^2+(8/9)x-(8/9)chi_3=0.                          (22)
```

Hence this candidate lies completely on the integral-recovery side.  It
cannot supply hidden-control nonintegrality or a pole chart in a normal host.

Conversely, the known nonintegral seed `(1)` has

```text
E_1=3/2+8e^2-4eu-8e^2u^2,                                 (23)
```

and cubic terms cannot change `E_1`.  On its exact pole chart

```text
u=t^(-1),
x=(2et)^(-1)-(8e^2)^(-1)-t/(4e)+vt^2,                      (24)
```

an affine cubic correction `(a_1u+a_0)x^3` has leading coefficients

```text
(a_1/(8e^3))t^(-4),             (a_0/(8e^3))t^(-3),         (25)
```

and the same statement holds in the second coordinate.  Since the old
coordinates are already regular, effectivity forces every affine cubic
correction to vanish.  The boundary density therefore remains `16e^2v`.

The simultaneous target is empty for two independently visible reasons:

```text
constant-density affine cubic side -> monic recovery and E_3 debt;
nonintegral effective pole side     -> immutable nonconstant E_1 debt.
```

This does not classify new pole charts unrelated to `(24)` or rational
coefficient corrections which decay at `u=infinity`.

## 6. Cross-thread connections and next tasks

The strongest connection is to the bounded normal-strip program.  THM-3861
already proves every polynomial Keller pair of transverse degree at most
three is an automorphism; the self-identifying arm in `(2)` therefore cannot
occur.  THM-3867 and THM-3871 extend that exclusion through degrees four and
five.  The present calculation supplies the local coefficient anatomy, but
it also shows that a **polynomial constant-density alternative seed should
jump directly to transverse degree six**.  Searching another cubic section
would duplicate a proved global closure.

The full-addition-ideal warning was also suggested by the three-cusp
polarization work, but no theorem from that lane is used here.  Incoming
THM-3881 is currently only a **RESERVED** cusp-ideal residual-transport signal;
it is neither a dependency nor evidence for any assertion in this note.

The generated tasks are therefore:

1. **Polynomial anchor:** start at transverse depth six and retain every
   determinant bucket; the arm collision supplies the hostile control.
2. **Rational niche:** on `(24)`, test coefficients with at least cubic decay
   at `u=infinity`, then audit every other denominator divisor.  Polynomial
   coefficients are now exactly excluded on this chart.
3. **Nonmonogenic lane:** replace one hidden primitive by two integral
   sidecars whose finite order is nonmonogenic; test the index/different unit
   rather than another resultant leading coefficient.
4. **Method wildcard:** if formal source additions are quotiented, carry the
   determinant cocycle and its divisor as part of the state.  A residual
   coefficient without that sidecar is not a lawful scheduler for Keller
   tasks.

## 7. Reproduction

```bash
python3 -B 04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py
python3 -B -O 04-computation/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.py
```

Both executions byte-match
`05-knowledge/results/jc2_russell_cubic_seed_quotient_constant_density_scout_20260823.out`.
The companion performs `56` optimization-safe exact checks.

```text
script SHA-256   3d67cc6674c0f432c4018d9cd03a0f646a9f9599b4cb4c00335c160049854a49
output SHA-256   8883668f8c8666bc30adcc1a89fc7e95bd80c295041b05f0e8c13096d0641d8f
semantic SHA-256 ad0debdf06adf279ee0f7e84c9e89b8257d2b14575f394c00b519a2d2f9121c9
```

The honest frontier remains `JC(2)` open.  This scout closes only the raw
affine cubic constant-density cell, the affine-relative cell over the full
quadratic section fibre, and affine cubic effectivity on the named
nonintegral pole chart.
