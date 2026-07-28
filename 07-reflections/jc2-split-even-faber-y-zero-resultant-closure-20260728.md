# The split `y=0` boundary is a degree-23 constant-field trap

**Status:** PROOF-COMPLETE + VERIFIED-EXACT CANDIDATE, awaiting independent
hostile audit as
`THM-2724-split-even-faber-y-zero-resultant-closure.md`.
This note concerns only the chosen-sheet split polynomial exact-square-prefix
degree-twenty-two **even-Faber** bank with `lambda!=0` and `y` identically
zero.

## 1. Inheritance pass

The closest proved algebra is THM-2411's pair of exact degree-twenty-two flux
numerators.  The closest corrected near miss is the temptation to use the
prime-23 scaled curve at `y=0`: the coordinate `t=rho/y` does not represent
that boundary at all.  Treating it as a point at `t=infinity` would silently
assume the very extension needing proof.

The least-used sidecar is again the third Faber observable.  It is not needed
to compute the response fibre; it becomes decisive only after the first two
fluxes force every response coordinate into the constant field.

The live board is

```text
y=0; q^2=T; Z=q^4; the linear u-equation G1;
the cubic u-equation N2; the degree-23 q-resultant;
the constant field C; the third flux R_Q'=kappa/U.
```

## 2. Why elimination is uniform here

At `y=0`, the chosen-sheet first flux is linear in `u`, while the second is
cubic.  Their exact elimination has the form

```text
Res_u(G1,N2)=nonzero_integer * P23(q).
```

Two coefficients carry the whole uniform conclusion:

```text
[q^23]P23=96059601,
[q^0]P23=-992137445376 lambda^3.
```

The first is independent of every coefficient parameter, so no exceptional
choice of `B,C,D,E,W` can make the eliminant identically zero.  The second
shows the chosen sheet never reaches `q=0` when `lambda!=0`.  The remaining
32 monomials describe the exceptional finite root positions but cannot
restore a nonconstant trajectory.

This is stronger than a generic resultant statement.  There is no leading-
coefficient divisor and no saturation chart to recover.

## 3. The two-stage constant descent

The eliminant first gives

```text
q algebraic over C => q in C*.
```

With `q` constant, the second flux retains its fixed nonzero cubic pivot in
`u`, so

```text
u algebraic over C => u in C.
```

Thus `T=q^2`, `d_ctr=u/T`, and `s_ctr=y/11` are all constants.  Only now is
the third flux spent:

```text
R_Q constant,                 but R_Q'=kappa/U!=0.
```

The mechanism is a response-to-continuation descent, not a genus argument.
It remains valid on every coefficient specialization in its exact chart.

## 4. Why 23 returns

The prime `23` reappears without introducing the scaled prime-23 curve.  The
chosen split sheet contributes one factor of `q` to the weight-five first
flux, `Z=q^4` contributes the fourth-power ladder, and eliminating the cubic
`u`-coordinate leaves the literal degree-23 polynomial `P23(q)`.  This is the
same weight invoice seen by THM-2704 from a different boundary chart, not an
identification of the two response curves.

That distinction matters:

```text
y!=0:  projectivize by t=rho/y and compare five q-pole fibres;
y=0:   eliminate u directly and force q into the constant field.
```

The two charts have complementary proof carriers.

## 5. Loss ledger

| operation | preserved | lost | repair |
|---|---|---|---|
| set `y=0`, `Z=q^4` | exact chosen-sheet flux equations | the signed scale `t` | eliminate in `(q,u)` directly |
| eliminate `u` | every necessary `q` value | the actual `u` root | return to the cubic `N2` after `q` is constant |
| conclude `q,u` constant | full response coordinates | continuation information | spend `R_Q'=kappa/U` |
| inspect only generic roots of `P23` | degree and typical separability | exceptional coefficient collisions | use the fixed leading pivot instead |

## 6. Exact boundary and next connection

The result closes

```text
y identically zero, lambda!=0, split exact-square prefix,
degree 22, even-Faber bank.
```

It leaves `y!=0`, `lambda=0`, the odd bank, and the broader split branch
open.  The natural synthesis is with the source pole-capacity theorem and
the arbitrary-`B` scale extension on `y!=0`.  Those ingredients should remain
separate until each has its own exact package and hostile audit; only then is
the reserved unified closure justified.
