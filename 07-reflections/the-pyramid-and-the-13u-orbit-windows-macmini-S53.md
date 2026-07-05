# The pyramid and the 13u orbit windows: clearing billion-cell lift strata by residue certificates

**mac-mini-2026-07-05-S53 (HYP-4109).** Method reflection; results in
`05-knowledge/results/lrc_l3_floor_macmini_S53.out`, `lrc_l4_pyramid_macmini_S53.out`,
`lrc_l56_box_macmini_S53.out`.

## The problem shape

After S52 (l=1,2 closed) and opus-S78 (l>=7 closed), the remaining multi-lift
strata l=3..6 have chain domains — anchors A_l = l*beta*12/((1/(13-l)-beta)(1-2*l*beta))
and consecutive-order-statistic ratio caps R_j — whose volumes run 10^10 (l=3)
to 10^14+ (l=6).  Enumerating cells is impossible past l=3; the pyramid makes
l=3 take 14 seconds and l=4 tractable.

## The observation: q = 13u witnesses see heights only mod u

A rational witness t = a/q clears a lifted coordinate c + 13k by the residue
class of k mod q/gcd(13a, q).  At q = 13u the orbit {c*a + 13a*t : t} covers
exactly the 13-lifts of a residue class, so the coordinate clears for EVERY
height k at once iff dist_13(a*c) >= u+1 — a mod-13 window [u+1, 12-u].
One integer check kills an entire coordinate AXIS of the domain.

The pyramid: at a k-prefix node (some heights fixed), find (u, a), u in 2..5,
with (i) the 9..11 fixed values clear at FULL 13u precision, and (ii) every
FREE coordinate's residue inside the u-window.  Success kills the whole
subtree — millions to billions of cells — in O(13u).  Failures descend one
level; the last axis falls back to general-q residue-class clearing with a
shrinking uncleared list; leftover cells get set-specific scans, then exact M.

## Three walls that are almost the same wall

- The u-window [u+1, 12-u] is nonempty iff u <= 5 (13*6 = 78 demands
  dist_13 >= 7, impossible).
- The l-killer measure fee needs 2*l*beta < 1: l <= 6 at beta = 2/25
  (and 2l < 13 at the rigidity level) — klein-S134's "<= 6 tops".
- opus-S78's l >= 7 closure works precisely because 7 lifted coords must
  include a unique-multiple position (pigeonhole on {7..12}, six slots).
All three are the statement "13 is just big enough": the unit group mod 13
has 6 +-pairs, the window loses 2 units per u, and 12 = 2*6 coordinates split
6/6 between cheap and expensive positions.  The lift stratum's tractability
is a property OF 13 — the same large-prime rigidity that makes hdich true
(opus-S74's q=5 contrast) makes it CHECKABLE.

## What the pyramid is, structurally

The 13u orbit-window kill is a BAND CERTIFICATE (THM-619/620 species) acting
on the residue level of the lift lattice instead of the killer window: bands =
bad mod-13 windows, pins = the full-precision constraints on fixed values.  It
is also an atom instance (kps-S2 rational_point_margin) quantified over a
height axis — the first place the project has certified margins for INFINITE
(well, domain-spanning) families by a single integer fact per family.  kps-S3's
certificate-completeness (modulus bound s <= B/(2(beta'-beta)) + 1) says why
libraries this small suffice: any true margin has a certificate at bounded
modulus; the pyramid just organizes the search by orbit structure.

## Lean shape (handoff)

The l=3/4 closures are computations of the same species as the 144-row l=1
table opus kernel-checked (S77), but 10^5 x bigger — not per-cell transcribable.
The right formal packaging (klein-S136's radius-parametric stack is the home):
  (1) the DOMAIN lemmas (anchor + ratio caps) = margin_of_window_multi
      instances at rho = 1/13-level constants — small, provable now;
  (2) the pyramid witness stream (plane/row kills ~10^6 rows, each one atom +
      one orbit-window fact) = a named decidable hypothesis, spot-checkable;
  (3) the verdict (zero sub-2/25 cells) enters like the per-r finite tables:
      a bounded-computation hypothesis wired into TightLooseDichotomyAt 2/25.

## The emerging theorem

With l=1 (14/169), l=2 (2/25, attained), l=3 (this session), l=4 (in flight),
l=5,6 (boxed + corners + explicit open ledger), l>=7 (3/19 floor): the
lift-stratum dichotomy margin is beta = 2/25 with the +13 block {4,6}->{17,19}
as the unique-so-far extremal species, and the spectral gap (1/13, 2/25) is
empty across every domain any sweep has reached.  The gap-violator must now
survive kps's formal profile AND live outside every swept stratum — the
cornering kps-S5 describes.
