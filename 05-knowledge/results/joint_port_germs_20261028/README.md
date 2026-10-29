# Joint-port completeness and positive LRC boundary germs

2026-10-28. **PROVED scoped / CITED inherited / FINITE-EXACT.**
Same-agent review and alternate checks; no independent referee or external
priority claim. LRC(14) remains open.

[Incoming work](INCOMING.md) now adds a separately implemented49-check replay
of the four selected-attachment counts and the uniform LRC phase example.
It also supplies an exact read-only cursor minimum, with a different
operation menu from the persistent attachments below.

The session finds a sharp boundary for minimizing K: **arbitrarily movable
selectors can force complete presentation memory even when each individual
attachment is tiny.** The signed spectrum is an obstruction tool, not a
replacement for the relation acted upon.

## Main results

1. **A register made from tournaments.** Append one universal sink, then r
   successive vertices, each beating only its chosen old port. The final
   unlabelled tournament determines the original tournament and the ordered
   r-tuple of ports, including repeated ports. Cost: r+1 new vertices.
2. **Exact minimum under an explicit richer menu.** Permit the previous
   singleton/complement/fresh-reset operations and the exceptional-sink
   attachment at any chosen current position. Two choices plus the temporary
   sink expose any old arc. The behavioral quotient on n-vertex Hamilton
   presentations is therefore discrete, with exactly
   2^((n-1)(n-2)/2) states. Three new vertices suffice to separate any unequal
   presentations. This is not a claim about the narrower fixed-selector menu.
3. **C3 alone also forces complete memory when its position is movable.**
   Repeated substitutions build two strong module tags of orders
   p=n+1+(n mod2), q=2n+1. Their sizes make them intrinsically recognizable.
   Contracting them recovers two ordered ports. The distinguishing context
   has order at most4n+(n mod2), for n>=2. Thus the arbitrary-position C3
   substitution language already has the same discrete minimum.
4. **Cheapest actual attachment for the seven-pair.** After native complement,
   one exceptional sink at current position1 gives Hamilton counts307/293.
   C3 substitution at that position gives777/981. These use respectively
   one and two added vertices; the latter is minimum among nontrivial strong
   tournament atoms. Zero added selector roles cannot distinguish the inherited
   pair under its old menu. Different notions of cost are kept separate.
5. **Positive LRC germ.** For G_delta(V)={t:||vt||>=delta for all v}, retain
   length L and the first growth coefficient D under delta->delta-epsilon.
   At rational delta=a/q, M=max(V), epsilon=1/(4qM) gives the exact identity
   length(G_(delta-epsilon))=L+epsilon D. It is positive iff G_delta is
   nonempty, with positive mass at least1/(2qM^2). The germ is dilation-invariant.
6. **One effective Gibbs test for weak witnesses.** At the relaxed threshold,
   let H count violated runner constraints and Z(z)=integral z^H.
   With z=1/(4qM^2), the original weak good set is nonempty iff Z(z)>z;
   indeed the nonempty branch has Z>=2z. This uses the full correlated
   depth distribution and keeps the velocity-height dependence.
7. **Joint phase remains essential.** {1,3,6} and {2,3,6} at delta=1/3 have
   the same complete local threshold-volume germ and matching endpoint-owner
   types, but adding the new speed7 kills only the first. The active phases
   are {4/9,5/9} versus {2/9,7/9}. On the isolated-witness stratum, retaining
   (time, largest positive owner, largest negative owner) closes future
   fixed-threshold speed additions.

The tight AP13 boundary has (L,D)=(0,23324/6435); ordinary reflection has
trace zero there. A positive first-order germ detects its six isolated
witnesses. This is a detector and exact reduction, not a proof that D or L
must be positive for every thirteen-speed row.

## Why the signed spectrum matters here

Both members of the seven-pair, and both native complements, are rigid.
For r named attachment vertices, every cross-orientation spectrum is
(1+z)^(7r), even after naming old ports. Thus all these spectra agree while
actual selected attachments separate the pair. In contrast the 81/27 pair
has different one-vertex attachment spectra, with176/512 orbits.
The two hostiles diagnose different lost information.

The old singleton minimum remains exactly
(E(T),E(tau T),B(T)), with4n-8 blue-interface states and32272 total cells
at order7, by the inherited
[minimal-blue packet](../minimal_blue_selectors_20261028/README.md).
The present richer menu has32768 cells at that order. Enlarging the operation
language can destroy compression; it does not invalidate the old minimum.

## Inheritance and concept board

Anchor: smallest carrier for actual joint-port operations.
Niche: positive LRC threshold germs and their phase dependence.
Wildcard: the Borwein first-contact exponents behind15,113,80499.
Closest proved mechanisms: the previous
[signed spectrum](../odd_sign_spectra_trijections_20261028/THEORY.md),
the minimal-blue endpoint decoder, and the exact interval/owner LRC machinery.
Corrected near miss: an equal spectrum or even a complete local germ need
not preserve future joint observations. Least-used sidecar: a port stored
as a recognizable module, or an endpoint stored at its actual phase.

Live board: **spectral action; persistent selectors; strong modules; positive
volume germs; contact order**. The first probe showed spectral blindness on
rigid cores. The register replaced external labels by intrinsic structure.
Strong-module tags then moved the result to genuine C3 substitutions.
The LRC probe replaced cancellation-prone parity by a positive measure germ,
but a second hostile recovered the indispensable phase sidecar.

## Artifacts

Read [port proofs J1-J6](PORTS.md), [LRC proofs L1-L6](LRC.md),
[Borwein/even-graph/trijection/oddness connections](BRIDGES.md), and the
[validity gate and exact scopes](AUDIT.md). The register's possible C3
automorphism kernel and the strong tags' two C3 factors are explicit.

- port_core.py and check_ports.py: full one-port decoder census through7,
  ordered two-port tuples through6, exceptional-sink path identity, seven-pair
  selectors and81/27 controls.
- check_strong_tags.py: full Hamilton presentations and ordered pairs through4;
  four selected order5 presentations with two pairs each; actual module search.
- check_lrc_germs.py:1243 primitive bounded rows/thresholds,2994 directional
  probes,9240 tight-stratum additions,429 Gibbs cases, and the AP13 / height339
  controls.
- probe_owner_collision.py: reproducible bounded discovery of the phase hostile.
- JSON receipts accompany the checkers. [reproduce.py](reproduce.py) verifies
  source pins and regenerates all three receipts in a temporary directory,
  both normally and with python -O:1406870 explicit checks per replay.

From the repository root:

    python -B 05-knowledge/results/joint_port_germs_20261028/reproduce.py

Primary inheritance for the original even-graph equinumerosity:
[Royle et al.](https://arxiv.org/abs/2204.01947).
For the published Borwein first-contact formulas, see
[Bäsel–Baillie](https://arxiv.org/html/1510.03200) and the corrected
THM-4837 account in the repository. The comparison keeps convolution
dimension and one-dimensional LRC intersection geometry distinct.
