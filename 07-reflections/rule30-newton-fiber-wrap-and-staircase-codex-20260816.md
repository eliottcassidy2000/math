# Rule 30: protected faces, failed gluing, dyadic wraps, and staircase entropy

**SYNTHESIS / IDEA PROVENANCE, 2026-08-16.**  Current truth is in the linked
theorems, not in this reflection.  None of the three Rule 30 prizes is claimed.

This session produced five complementary exact objects rather than one scalar
attack on the center sequence.

## 1. The inward current has a protected Newton face

[THM-3488](../01-canon/theorems/THM-3488-rule30-inward-slack-monicity-and-parity-cartier-ramification.md)
retains Green-transport slack `v`.  After the first three boundary strips are
removed, the complete inward tail has a unique lowest ballistic face

```text
t-v=4,       Face_4=z^4/(1+zq).
```

Consequently every target slack polynomial is monic.  The protection is
genuinely graded: setting `q=1` can cancel it, and small physical targets
already do.  The useful analogy is a lowest tropical/Newton face, with the
same warning learned from factorial boundary currents: an unlabelled boundary
sum may vanish while its transverse jet is nonzero.

## 2. The marked center is a different transverse face

[THM-3489](../01-canon/theorems/THM-3489-rule30-packed-restart-and-pointed-pascal-face.md)
shows that a proper terminal arc reads the center from a sparse set of high
phase-Hasse coordinates.  Those coordinates lie strictly above the entire
low repeated-root sector lost by arc integration.  At wrapped depths the
packed orbit restarts and the center closes from the lift bit; at hard depths
the calibrated high face remains.

Slack Hasse order and phase Hasse order are therefore not two names for one
coordinate.  One belongs to Green transport, the other to cyclic phase.

## 3. The two exact marginals do not canonically glue

[THM-3492](../01-canon/theorems/THM-3492-rule30-fiber.md) gives the minimal
common object.  The terminal profile `F` and center-completed slack polynomial
`D` form a fiber product over their common marked bit.  Bivariate lifts fit
the exact sequence

```text
0 -> (q+1)P_(N-1) tensor ker(mark)
  -> P_N tensor V_p
  -> V_p x_(F_2) P_N -> 0.
```

Even after every coefficient profile is required to lie in the physical arc
image `I=Y^ell V_p`, the mixed ambiguity has dimension

```text
N(p-ell-1).
```

Every independently load-bearing pointed Pascal coordinate can carry the same
monic slack face.  Thus the missing datum is not another marginal statistic:
it is a raw eventwise slack-by-phase coupling, or an equivalent chain-homotopy
sidecar.  Depth five is the first two-carrier hostile; depth six is the first
one with a nontrivial arc-image constraint.

## 4. Wraps have an exact dyadic geometry

[THM-3493](../01-canon/theorems/THM-3493-rule30-dyadic-wrap-atlas.md) starts
from the persistent physical left boundary pair `11`.  It sharpens the period
floor to

```text
P_k >= 2^floor(log_2 k),
```

so the nominal endpoint `k=2P_k` never occurs.  In each dyadic block
`[2^m,2^(m+1)-1]`, wrapped depths are either absent or form one initial prefix
whose complete center word is `0...01`.  Hence wrapped ones contribute at
most one per dyadic scale; balance, if true, must come from the hard depths.

The finite-exact scout certifies every depth from `5` through `2^28-1` as
hard.  This large interval is a test surface, not asymptotic evidence.
Moreover, the wrapped-one reciprocal mass is at most `2`, so any divergent
harmonic signal in the center support transfers entirely to the hard core.

## 5. A smaller universal boundary improves the compiler

[THM-3491](../01-canon/theorems/THM-3491-rule30-seven-four-staircase-compiler.md)
pushes the exact staircase factor language to height 21.  A positive integer
certificate `4Av<=7v` proves staircase entropy at most `log_2(7/4)`, improving
the charged one-word macroblock query tariff from `7/2` to `13/4`:

```text
Q(n)=(13/4+o(1)) n^2/log_2(n)^2.
```

This is an all-`n` upper compiler in a declared word-RAM model.  It is not a
fixed-seed lower bound.  The finite `CABC/BCBC` block code is the live signal
for the opposite entropy direction, but no all-height encoder is known.

## 6. Updated next targets

1. **Physical coupling:** construct the raw space-time refinement mapping to
   both the slack Newton grading and calibrated phase face; determine which
   mixed-kernel class Rule 30 selects.
2. **Innovation valuations:** control whether
   `v_m=nu_2(R_(2^m)-1)` is eventually below `2^m`, or at least its frequency.
   The wrap atlas makes this the exact scalar gate for hard-depth density.
3. **Hard-face correlation:** evaluate the Mersenne pointed face on actual
   Rule 30 currents, rather than on the ambient arc module.
4. **Staircase lower encoder:** prove or refute an all-height synchronizing
   code extending the finite `CABC/BCBC` family.
5. **Compiler composition:** combine the staircase macrostep with dyadic wrap
   detection without charging the same simulated row twice.

The common lesson is now precise: preserve a protected grade, preserve the
physical mark, and invoice their coupling separately.  Any method retaining
only two of those three coordinates is stopped by an exact hostile above.
