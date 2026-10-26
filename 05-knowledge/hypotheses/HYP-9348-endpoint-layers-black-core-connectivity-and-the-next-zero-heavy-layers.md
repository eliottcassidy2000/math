---
id: HYP-9348
title: "Filtered endpoint-four detours and black-core connectivity; zero-heavy layers through the first weighted endpoint"
status: OPEN corrected endpoint/black-core questions; PROVED closed endpoint-four n>=8 with finite heads; REFUTED unfiltered (a); SUPERSEDED old deficit-(k+1) frontier
source: mac-mini-2026-10-20-debtcones
depends_on:
  - 01-canon/theorems/THM-4734-endpoint-three-layer-closed-at-order-ten-and-the-endpoint-four-layer-at-orders-nine-and-ten.md
  - 01-canon/theorems/THM-4736-the-delta-equals-k-zero-heavy-layer-exact-counts-missing-gaps-and-the-adjacent-matched-pair-collapse.md
  - 01-canon/theorems/THM-4733-divisor-plane-of-p2qr-versus-the-marked-path-tiling-cube-3-4-3-and-4-2-4-reflections.md
related:
  - 05-knowledge/results/endpoint_three_repairs_20261020h/README.md
  - 05-knowledge/results/charge_palettes_boundary_20261020g/zero_heavy.md
  - 05-knowledge/hypotheses/HYP-9345-merged-metagraph-open-laws-after-thm-4729-4730.md
---

# HYP-9348 -- the next layers

**Updated 2026-10-25:** the unfiltered statement in (a) is refuted below;
(b) retains its non-pure-blue domain; (d) now routes the already proved
corridor and the new weighted-endpoint branch. Original statements remain
in Git history; the failure mechanisms are recorded in `01-canon/MISTAKES.md`.

## (a) The structural (n-5,4) theorem

**REFUTED as originally unfiltered.** At `n=7`, source
`(6,4,8,19,103,79,15)` has endpoint scores `(2,4)` and a closed blue standard
word. Its literal target is `(122,124,121,112,32,64,16)=C3 -> 1 -> C3`,
whose nine Hamilton paths are all blue. Thus no black detour reaches it.
The [exact witness](../results/shape_ports_20261025/tournaments.md) also
repairs THM-4734's analogous order-six omission.

**PROVED closed part with sharp unfiltered threshold.** Every closed blue endpoint-(n-5,4)
marked path at `n>=8` now has three directed black moves to its exact
labelled K target, with no additional source filter. The
[head assembly](../results/lca_ports_reverse_20261025/HEAD_ASSEMBLY.md)
combines the structural tail, complete affine12..14 repairs, incoming
order11 certificate, and freshly replayed complete finite heads8..10.
Arbitrary shared/reflected ports and clone configurations are retained
literally. The order7 hostile above prevents lowering the unfiltered
uniform threshold. The proof's finite certificates remain explicit.

**REFUTED old length-three extension at7:** the inherited F7 distance-five
exception in [THM-4732, section3](../../01-canon/theorems/THM-4732-black-detours-of-closed-blue-lines-order-9-census-order-10-detours-pure-blue-layers-and-the-klein-dictionary.md)
already prevents that bound even after excluding pure-blue targets.

**OPEN repaired question.** Does every such line whose literal target is
not pure blue admit a black detour with an appropriate length at the
remaining small orders? A universal three-step statement at7 must also
respect the inherited F7 distance-five exception. Higher endpoint layers
remain open. Passing an open blue word to its closed K partner changes
the layer by `k -> n-2-k`; open-four becomes closed layer `n-6`, so this
fixed closed-four theorem does not settle arbitrary open-four words.

## (b) Black-core connectivity

**Claim** (the lane's conjecture; THM-4731 verifies it for `n <= 8`): all non-pure-blue
classes lie in one black component of the line metagraph. The endpoint layers remove one
family of potential bridges within their admitted domain: `k=2` has its
inherited theorem; `k=3` has all `n>=7` and eleven admissible order-six
presentations; `k=4` has all closed lines at every `n>=8` by (a).
Pure-blue targets cannot be removed by a black detour. In the trijection coordinates (graphs on ordered `[n-1]`; THM-4733 Proposition A)
a K move along `Q` is "change the marked path, then complement"; a three-black detour is a
closed walk of three coordinate changes whose complement bits cancel to `E(P)`. Question:
is there a `k`-independent movable-cut word, written in these coordinates, that repairs every
closed line whose port set is not an interval? The packet's section-9 words are the case
`k = 3`; the half-covering template of THM-4734 section 2 is their `k = 4` shadow.

## (c) Fixed-word obstructions and the next port layer

The old second fixed-word obstruction (`n=17`, reflection-product orbits
`3,3,5,3,3`) is bypassed by the movable cut at `k=3`. Closed `k=4` lines
at17 are covered by (a). The former six failures at12/14 are now repaired
on their full affine cells: a three-arc selector handles12 and a one-arc
selector handles14. The failed-word automorphisms force those guards.
Their original hostile sources remain in the new controls.

The next unrestricted layer has four incoming interior ports, endpoint
scores `(n-6,5)`. The old argument forcing at least four ports no longer
contradicts that budget. A new repair or reduction between layers is
needed. Separately, THM-4732 Lemma2 still requires connectivity of the
line metagraph after pure-blue vertices are deleted; local redundancy
in one degree layer does not establish that global condition.

## (d) The corrected zero-heavy boundary

**SUPERSEDED:** the old “next deficit k+1” question is closed by the
[robustness packet](../results/zero_heavy_robustness_20261021a/zero_heavy.md)
and [two-unit/first-owner extension](../results/zero_heavy_extension_20261021a/README.md).
They reach deficit `2k` in their stated parameter ranges, with exact
small-parameter controls. The old paragraph incorrectly allowed two
missing gaps at `k+1`: the exact identity

    delta=(k+1)d+kh+U

charges two missing gaps `2k>k+1` for every `k>=2`. Shared endpoints first
appear at `k+1`, while two missing gaps first appear at `2k`.

The [new weighted-endpoint branch](../results/shape_ports_20261025/tournaments.md)
closes `r<k` at deficit `2k+1`: for `k>=3` it has
`2(k-1)(2k-3)!!` classes generally and `2(k-1)I_(k-1)` strongly. At `k=2`
two marked records merge into one class. The literal pointed-triangle
collapse proves `H_k^s=2(k-1)I_(k-1)` while retaining its external port.
The other `r=k` families at this deficit, the unrestricted owner grammar,
and a smaller unrestricted structural A000568 recurrence remain **OPEN**.
The ordinary exact reverse transfer and forced-skip contraction are now
[implemented and cross-checked](../results/lca_ports_reverse_20261025/counting/INCOMING_COUNTING.md);
they preserve active defect distinctions and do not supply that owner grammar.

## (e) Arithmetic Klein sets and the fixed-tile-free cube

THM-4733 shows `D(p^2 q r)` is a fold and a sub-object of the four mirror-moved tiles of the
order-five cube and cannot embed in the full cube. **Question.** For `N = p^k q r` (the
packet's `(2a - k(b+c), b-c)` chart, `4k + 2` points) is the corresponding object the fold of
the `2k` mirror-moved tiles of a suitable order along `k` orbits, and which orders `n` have
`floor((n-1)/2) = f` fixed tiles compatible with the sign split `2k, 2, 2k` (even `k`) versus
`2k+1, 0, 2k+1` (odd `k`)? The odd-`k` case has no `CQ`-fixed divisor and therefore has no
fixed-tile obstruction: it may embed in genuine presentations. First test: `p^3 q r` against
the order-six cube.
