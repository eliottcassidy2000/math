# The Pell-predicted packet survives the next finite-sheet gate

**Current status: PROMOTED AS THM-3527 AFTER A SPLIT-OUTER AND FROZEN-`H`
AUDIT.**  At the canonical finite old-`L` inverse
point

```text
q=(2,5/6,-7/8),
```

the complete six-layer recursive norm evaluator gives

```text
R_7(q)!=0.
```

All three primary reductions and the disjoint split reduction are nonzero.
Together with the proved `R_7` packet `A(66907,24255)`, THM-3527 gives

```text
v_L(N(R_7))=-66907,
R_8=L^66907 N(R_7) polynomial and L-coprime,
R_8 has A(419839,152211).
```

The final pair is exactly the next row predicted by the Pell-57 recurrence,
but the recurrence did not prove this finite-sheet gate.  No image,
irreducibility, all-level, arbitrary-map, or general Jacobian-conjecture claim
is made here.

## Inheritance pass

- **Closest proved mechanism:** THM-3523 splits the old-`L` valuation of
  `N(R_6)` into two complete-face divergent sheets and one finite sheet,
  proving `R_7` polynomial and `L`-coprime.
- **Canonical hostile:** the finite branch could vanish even when the two
  divergent valuations are completely known.  In that case nominal
  denominator clearing would overstate the exact valuation.
- **Corrected near miss:** THM-3522 propagates packets only after polynomiality
  has been proved.  The Pell-57 sidecar predicts the arithmetic row but pays
  no denominator gate.
- **Least-used sidecar:** one more transitive norm factor `N^6(L)`, represented
  in a complete `729`-dimensional inverse algebra.

## Exact norm-orbit identity

The fixed definitions are

```text
H   =2^6 L N(L),
J   =2^35 L^7 N(H),
G   =L^43 N(J),
R_5 =L^271 N(G),
R_6 =L^1699 N(R_5),
R_7 =L^10663 N(R_6).
```

Norm multiplicativity gives, in the localized function field,

```text
R_7
 =2^4293 L^10663 N(L)^1699 N^2(L)^271 N^3(L)^43
          N^4(L)^7 N^5(L) N^6(L).                   (1)
```

The scalar exponent is not guessed from the packet recurrence:

```text
4293=35*3^4+6*3^5.                                   (2)
```

Thus the new work relative to THM-3523 is exactly one additional inverse
cubic and one additional orbit factor.

## Complete 729-sheet calculation

The implementation hash-pins the recursive-adjugate engine independently
audited in THM-3526.  For each of `p=101,103,107`, it starts at the reduction
of `q`, constructs six nested cubic algebras of dimensions

```text
3,9,27,81,243,729,
```

and checks all six inverse graph identities by direct substitution.  At every
level it also certifies four units: the leading `L`, cubic derivative,
`y`-chart denominator, and `x^3` denominator.  All `72` displayed unit values
are nonzero.

The seven-factor norm orbits are

```text
p=101: (16,12,72, 9,49,97,71),
p=103: (12,53,22,85,76,94,17),
p=107: (38,45,28, 3,17,17,30).                       (3)
```

Substitution in (1) gives

```text
R_7(q) mod (101,103,107) = (72,44,53).                (4)
```

Any one entry of (4) proves nonvanishing if the representation is accepted;
all three are retained as hostile redundancy.  Omitting the factor `64` from
the definition of `H` removes it on all `243` bottom `H` leaves and changes
the residues to

```text
(66,78,19).                                           (5)
```

The exact exponent change is `64^-243`, equivalently
`2^4293 -> 2^2835`.  This checks both the nonmonic normalization and the leaf
multiplicity.

## Why the Pell agreement matters, and why it is not a proof

THM-3522 gives the packet transition

```text
(e,m)->(7e-2m,3e-2m).
```

The Pell-57 sidecar identifies this arithmetic recurrence with multiplication
by `(5+sqrt(57))/2`.  Starting from the proved `R_7` row gives

```text
(66907,24255)->(419839,152211).                       (6)
```

Equation (6) was known before the 729-sheet run.  What (3)--(5) add is the
missing geometric gate: the finite old-`L` branch does not vanish.  The two
pieces have different logical roles:

```text
Pell/face recurrence: predicts the next packet if polynomiality is paid;
finite-sheet norm:    pays the exact old-L denominator gate at this rung.
```

Their agreement is therefore a successful prediction/realization test, not
an induction.  The next rung would again require a new finite-sheet unit.

## Proved consequence and exact boundary

On the two divergent sheets above the generic old-`L` divisor, the complete
`R_7` face gives valuation `-66907/2` on each.  A nonzero finite specialization
shows that the regular third sheet has generic valuation zero.  The split
audit realizes all `3*243` sheets, so the same UFD argument as THM-3523 proves

```text
v_L(N(R_7))=-66907,
R_8=L^66907 N(R_7) in Q[a,b,c],
gcd(R_8,L)=1.
```

THM-3522 then turns (6) into the complete `R_8` packet.

The theorem does not prove that `R_5`, `R_6`, `R_7`, or `R_8` is an
image prime; does not add a nonproperness component; does not prove
irreducibility or squarefreeness of a cleared norm; and does not prove an
unconditional all-level tower or any general Jacobian-conjecture statement.

The disjoint audit uses the split-outer representation from the previous rung.
Over `F_71`, the outer cubic above `q` splits at `w=10,23,38`; the complete
243-sheet branch values are `56,65,10`, whose product gives
`R_7(q)=56 mod 71`.  Omitting the last branch gives `34`.  A frozen-`H`
bottom determinant independently checks each branch value.

## Reproduction

```text
python -B 04-computation/keller_R7_finite_sheet_recursive_norm_probe_20260816.py
python -B -O 04-computation/keller_R7_finite_sheet_recursive_norm_probe_20260816.py
python -B 04-computation/keller_R7_finite_sheet_split_outer_independent_audit_20260816.py
python -B -O 04-computation/keller_R7_finite_sheet_split_outer_independent_audit_20260816.py
```

The pinned semantic digest is
`82efb24e0c4a6e0df9671f0f5a5009dd0e77d1b0aa8ef2341780dfe23ea28c38`.
