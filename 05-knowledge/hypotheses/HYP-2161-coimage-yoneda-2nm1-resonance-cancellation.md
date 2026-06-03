---
id: HYP-2161
status: synthesis / proof-program
source: user-2026-06-03; codex-2026-06-03-S605
related:
  - HYP-2083
  - HYP-2084
  - HYP-2135
  - HYP-2138
  - HYP-2142
  - HYP-2153
  - HYP-2154
  - HYP-2155
  - HYP-2156
  - HYP-2157
  - THM-401
  - THM-406
---

# HYP-2161: coimage + Yoneda; `2n-1` resonances are the cancellation

## Claim

The `p_0=0` collapse should be read through two equivalent coimage
presentations:

```text
p_0 = sum_j (-1)^j S_j
    = sum_{c in L(V)} prod_i kappa(c_i).
```

THM-406 gives the overlap/moment coimage: `S_j` are the total `j`-fold overlap
volumes, and `p_0` is their full alternating sum.  The S604 resonance-lattice
work gives the harmonic coimage: `L(V)={c in Z^m : c.V=0}` indexes the Fourier
correction terms to the free value.

The new synthesis is:

```text
coimage = the minimal quotient preserving the observable;
Yoneda = enough probes determine that quotient;
C=2n-1 shell probes = the natural probes for the LRC floor stratum.
```

Thus "`2n-1` resonances are the cancellation" means: at the LRC floor, the
resonance corrections that cancel the free/Poisson ground cell are organized by
the canonical pair-sum shell modulus `C=2n-1`.  The shell ledger is not a
decorative arithmetic side condition; it is a candidate finite presentation of
the all-orders coimage cancellation.

## Why `2n-1`

THM-401 proves that `C=2n-1` is simultaneously:

```text
1. the Farey successor denominator of the floor: 1/n -> 2/(2n-1);
2. the odd summand-shell modulus P_a={a,C-a};
3. the odd-square companion sqrt(8*binom(n,2)+1).
```

So the pair-sum witness clock, the additive shell clock, and the triangular
pair-count clock all name the same object.  This is exactly the kind of
over-determination that the coimage/Yoneda reflection predicts: many probes
keep re-deriving one quotient.

## Cancellation Equation

In the harmonic presentation, the free term is positive:

```text
free = (1-2 delta)^m.
```

Tightness requires the nonzero lattice vectors to cancel it:

```text
sum_{c != 0 in L(V)} prod_i kappa(c_i) = -free.
```

Relation-poor `V` cannot do this.  AP rows, additive chains, and composite-shell
sporadics can.  HYP-2161 says the first finite arithmetic quotient worth testing
is not a low moment, but the `C=2n-1` shell probe system and its refinements:
unit shells, nonunit shells, lift denominators, endpoint owners, and bounded CRT
states.

## Yoneda Form

Let `Q(V)` be the depth coimage of a speed set at the floor, and let
`Shell_C(V)` be the labelled `C=2n-1` shell ledger.  There is a forgetful map

```text
Q(V) -> Shell_C(V)
```

that records the responses of the coimage to the `C`-shell probes.  The desired
Yoneda-style theorem is that this probe family is conservative on the worry
stratum:

```text
same Shell_C data + same lift/CRT refinements
  => same p_0-collapse decision.
```

Equivalently, the shell probes represent enough of `Hom(-,Q)` to determine the
ground-cell question.  If true, the `2n-1` ledger is a finite presentation of an
otherwise all-orders cancellation.

## Consequence for `n=14`

For `n=14`, `C=27=3^3`.  Prime `C` rows mainly test unit shell transversals, but
composite `C` creates nonunit strata and lift denominators.  The proof target is
therefore:

```text
prove the C=27 unit/nonunit/lift/CRT ledger is conservative for p_0 collapse.
```

This converts "control all overlap orders" into a finite certificate program:

1. Missing a visible unit shell gives a `2/(2n-1)` witness.
2. Nonunit holes must lift to denominator or CRT probes.
3. Additive-chain/shell rows are the only known ways to keep all probes silent
   while cancelling the free ground cell.
4. The desired `n=14` certificate is an upper-cap statement on the `C=27`
   ledger: every candidate either exposes a witness probe or lies in a
   classified cancellation family with the unit boundary floor intact.

## Assumption Challenge

Do not make runners or arcs the default vertices.  The possible vertex/probe
sets are:

```text
depth cells p_k,
overlap orders S_j,
relation-lattice vectors c in L(V),
C=2n-1 summand shells,
unit/nonunit shell gaps,
lift denominators,
endpoint-owner obligations,
bounded CRT states,
two-block determinant components,
proof obligations.
```

The quotient must preserve the `p_0=0` decision and, for LRC, the existence of a
floor witness.  It deliberately destroys raw phase order and individual runner
identity.  The challenged assumption is that the cancellation must be explained
by low-order density or by runner-level geometry; HYP-2161 instead asks whether
the `2n-1` probe ledger is the conservative Yoneda presentation of the coimage.

## See

`07-reflections/lrc-coimage-yoneda-2nm1-cancellation-s605.md`,
`07-reflections/coimage-yoneda-2n-minus-1-resonance-s605.md`,
`07-reflections/lrc-coimage-fundamentality-made-rigorous-s599.md`,
`07-reflections/lrc-resonance-lattice-pvsnp-s604.md`,
`05-knowledge/hypotheses/HYP-2157-overlap-order-partition-function-calculus.md`,
`01-canon/theorems/THM-401-pair-sum-sieve-modulus-is-2n-minus-1.md`,
`01-canon/theorems/THM-406-covering-depth-master-object-factorial-moments-and-spectral-identity.md`.
