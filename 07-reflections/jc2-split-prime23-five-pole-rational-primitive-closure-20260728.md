# The prime-23 response curve is closed by a source-side pole budget

**Status:** PROOF-COMPLETE + INDEPENDENTLY HOSTILE-AUDITED CANDIDATE under
the slug
`THM-2718-split-prime23-five-pole-rational-primitive-closure.md`.
The bare theorem ID is temporarily quarantined by a later duplicate
reservation.  This note concerns only the `B_0!=0`, `y` not identically zero,
`lambda!=0`, split polynomial exact-square-prefix degree-twenty-two
**even-Faber** scale chart.  It does not include the odd bank or the full
split branch.

## 1. Inheritance pass

The closest proved response-space mechanism is THM-2713.  Every prime-23
curve on the chart is integral, and its normalized chosen-sheet function has

```text
div(q)=5N-3O,
```

with five distinct old poles of order three.  The canonical hostile member
is still an abstract curve with rational normalization and exactly `89`
extra delta units; THM-2713 deliberately does not exclude it.

The closest proved source-space mechanism is THM-2214's rational-primitive
lemma.  If a polynomial `U` satisfies

```text
U S'=nonzero constant,                 S in C(x),
```

then `U` is constant or is supported at one finite point.  This was used in
older nonsplit closures only after explicit normalization of a singular
spectral curve.

The underused sidecar was the third Faber flux.  On a split sheet it becomes
more, not less, elementary:

```text
q=A_src/U in C(x),             U R_Q'=kappa.
```

Nonsplit parity is irrelevant to this differential equation.  It had been
used only to replace the constant first flux by zero.

## 2. The missing comparison was across spaces

The response-curve and source-line statements initially look incomparable:

```text
X=normalization(C):       q has five pole points;
P1_x=source line:         the denominator U has at most one finite root.
```

A physical trajectory supplies the map between them.  Once the map is shown
nonconstant, it is finite and surjective.  The five old points therefore
pull back to five disjoint nonempty fibres.  At most one fibre can contain
the source infinity, so at least four finite source points are poles of
`q=A_src/U`; all four must be roots of `U`.  The primitive lemma permits at
most one.

This bypasses the large perfect-power system

```text
t=alpha beta^4/H,       q=alpha^5/beta^3,
zeta=alpha^23/H^3,      F1=beta^23/H^5.
```

Those identities remain the correct description of an abstract rational
member, but a physical Keller lift is killed before they need to be solved.

## 3. Why the constant-map case matters

One cannot silently call the source map dominant.  If its projective image
were constant, the pole pullback would say nothing.  The signed scale repairs
this point:

```text
t=rho/y,        v=u/y^2,        zeta=Z/y^3,        rho!=0.
```

Constancy of `[1:t:v,zeta]` makes `y,u,Z` constant.  The split first flux
then makes `q,T`, and the centered quartic coordinate `d_ctr=u/T` constant.
The third observable `R_Q` is consequently constant, contradicting
`R_Q'=kappa/U`.  Dominance is therefore a theorem, not an implicit physical
assumption.

## 4. Loss ledger

| operation | preserved | lost | repair |
|---|---|---|---|
| project to the prime-23 curve | two constant flux equations and the `q` divisor | the Keller differential | retain `R_Q'=kappa/U` |
| normalize the response curve | the five distinct old pole points | the source polynomial denominator | pull back along the physical map |
| write `q=A_src/U` | every finite source pole must lie over a root of `U` | cancellation multiplicities | use valuations `ord(U)>ord(A_src)>=0` |
| solve abstract rationality by `alpha,beta,H,V` | all perfect-power divisors | whether the curve lifts physically | compare pole capacity first |

The cancellation line is worth retaining.  It is not necessary to claim
that `A_src/U` is reduced.  At an actual pole, negative valuation directly
forces strictly more vanishing in `U` than in `A_src`.

## 5. Structural lesson

The successful move is a **cross-space divisor budget**:

1. extract an intrinsic divisor with many distinct support points on a
   response curve;
2. retain the differential equation that constrains the support of a source
   denominator;
3. prove the physical trajectory is dominant;
4. pull distinct target support back to disjoint source fibres; and
5. compare support cardinality before solving either space explicitly.

This pattern is stronger than genus when the physical sidecar is rigid.  A
rational response state may exist while no continuation state can traverse
it.  It should be tested in the degree-eighteen discriminant strata and in
the remaining degree-twenty-two odd-bank charts whenever a flux quotient has
three or more forced poles.

## 6. Exact boundary and next work

The proof closes only

```text
B_0!=0, y!=0 as a function, lambda!=0,
split exact-square prefix, degree 22, even-Faber bank.
```

The next high-value lanes are orthogonal:

1. restore the eleven odd Faber coefficients and ask whether a divisor with
   at least three forced `q`-poles survives;
2. close `B_0=0` without the `rho` scale;
3. close the invariant boundary `y=0`; and
4. analyze `lambda=0` on the split deck without importing nonsplit parity.

The abstract rationality question for special prime-23 curves remains
mathematically interesting, but it is no longer a physical obstruction on
the chart just closed.
