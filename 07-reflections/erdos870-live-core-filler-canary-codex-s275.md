# Erdos-870 Live-Core/Filler/Canary Deletability Audit

Anchors: HYP-3148, HYP-3147, HYP-3146, HYP-3145, HYP-3144, HYP-3143, HYP-3142,
HYP-3141, HYP-3140, HYP-3134, HYP-3133, HYP-3124, HYP-3054, HYP-2534,
T1213, LTI-274, LTT-172, OPEN-Q-108.

The useful abstraction is not that the n=4 tournament table has four colored
classes.  It is that two equivalent-looking tables expose two different
proof economies.

S276/HYP-3143 already extracts the exact-order lesson: fixed filler plus a
free obstruction basis must exclude lower-order representations.  Incoming
S274/HYP-3144 adds the pair-function lesson: unordered quotients keep only the
symmetric payload unless ordered sidecars are retained.  Incoming HYP-3145
adds the filler-core interface: fixed-path tilings are atlases/alarms, while
terminal rows should expose the small retained core.  Incoming S274/HYP-3146
adds the shift-package/canary policy: keep redundant fiber mass for deletion
stability or scaffold it away for quotient congruence.  Incoming S277/HYP-3147
adds the local edge-flip/Worpitzky/function sidecar.  This note adds the
adjacent HYP-3148 audit: the same quotient must also say which coordinates are
live, filler, canary, and deletable.

The fixed-path tiling cube keeps three live skips `a,b,c`.  It sees the four
classes, but its full cube is skewed: `T:+:-:S = 1:1:1:5`, and `c` is
class-cover-deletable because `{a,b}` already reaches every class.  The
two-bit anchor freezes the long diagonal as filler/canary data.  Then the
live `x,y` square has uniform class counts `1:1:1:1`, and both live bits are
load-bearing.

That is the Erdos-870 warning in miniature.  Many representations, many
witnesses, or many class encodings do not imply minimal proof support.  A
quotient is proof-facing only after it says which coordinates are live, which
are deterministic filler, which are canary/shift controls, and which are
deletable.

For the LRC14 frontier this suggests a practical ledger rule.  Every
tournament, A000568, edge-witness, fiber-PGF, or k=8 sidecar row should add:

```text
live_core_bits
filler_bits
canary_bits
deletable_coordinates
class_distribution
minimal_cover_subbasis
edge_bounded_core_floor_exit
terminal_exit_or_named_debt
```

The creative reframe is to stop asking whether a small table is the right
model and ask which proof economy it certifies.  Scheme A is exploratory: it
shows where witnesses proliferate and which coordinate is cheap.  Scheme B is
the certificate form: it names a minimal live core and makes filler visible.
The next LRC14 test is to run the same audit on HYP-3141 edge packets and
HYP-3140/HYP-3142 coefficient/moment packets: if a coordinate is deletable
there, the proof should not pay for it; if it is not deletable, it must be
kept until a named sidecar closes the obligation.
