# LRC sections as proof objects: Hall packets, wall debt, and switch tournaments

**Source:** codex-2026-06-17-S3.  Prompt: think about regions of the loop, such
as quadrants for LRC4, instead of runners; investigate whether a proof can check
off the `n` subdivided sections and compress exotic cases to relations or
switches with a tournament analogy.

## External compass

The literature compass still supports this kind of quotient.  Tao's Bohr-set
view frames LRC as a covering/avoidance problem on the circle, so a section
ledger is a natural geometric shadow.  Malikiosis, Sgall, and Somer describe a
finite-checking/zonotope style reduction, which again points toward finite
cellular decompositions rather than runner labels alone.  Recent arXiv progress
has pushed the verified frontier upward toward LRC14, so the live problem is to
find the right retained finite state rather than a new scalar invariant.

References searched this session:

- Tao, "A view-obstruction problem of Bourgain, and the lonely runner conjecture":
  https://terrytao.wordpress.com/2017/05/08/a-view-obstruction-problem-of-bourgain-and-the-lonely-runner-conjecture/
- Malikiosis-Sgall-Somer, "Linearly-exponential checking is enough for the
  Lonely Runner Conjecture": https://arxiv.org/abs/2411.06903
- Rosenfeld, "The Lonely Runner Conjecture holds for Eight Runners":
  https://arxiv.org/abs/2509.14111
- Sungkawichai-Trakulthongchai, "Eleven, twelve, thirteen lonely runners":
  https://arxiv.org/abs/2604.23906

## The section graph

Fix the slowest-runner gauge and include the stationary runner:

`V=(0,v_1,...,v_{n-1})`.

Divide the loop into `n` sections.  For each runner `r`, add an edge from
runner `r` to section `s` if there is a time at which runner `r` is lonely and
lies in section `s`.  Boundary witnesses are counted in the compactified graph
for both adjacent sections; strict-open witnesses use only open cell interiors.

The user's dream is now a precise Hall problem:

> Can every runner be matched to its own lonely section?

This is a useful reduction because "four runners each use a different quadrant"
is exactly a perfect matching in the LRC4 runner-section graph.  Overlapping or
reused regions are not a separate pathology; they become Hall packets.

## What the scout found

The exact rational scout `lrc_section_checkoff_switches_codex.py` scans open
cells and wall events.

Small primitive scans:

| n | max speed | primitive rows | compact matching | open matching |
|---:|---:|---:|---:|---:|
| 4 | 14 | 325 | 325/325 | 0/325 |
| 5 | 10 | 205 | 205/205 | 0/205 |
| 6 | 8 | 56 | 56/56 | 0/56 |

This is the central lesson: compactified sections behave exactly like the
check-off dream in the tested range, but strict-open sections expose wall debt.
So the section proof does not fail at overlap; it fails at equality walls, which
are already the natural habitat of the dihedral endpoint machinery from
HYP-2569.

For LRC4 AP `(0,1,2,3)`, the compact matching exists, but only two sections have
positive open measure.  The open Hall packet is

`deficit=2, runners=(0,3), sections=()`.

This says the first and last runners are checked off only by wall witnesses in
endpoint sections.  For LRC4 uneven `(0,1,2,4)`, the open Hall packet is

`deficit=1, runners=(1,2,3), sections=(1,2)`,

with wall debt in section `3`.  That is exactly the "multiple runners want the
same region" exotic case, now compressed into a finite relation.

For the LRC14 AP row and the Goddyn-Wong row, compact matching also succeeds and
strict-open matching fails with open cover `12/14`; the debt sits at endpoint
sections.  This is reassuring: the section picture sees the same endpoint wall
phenomenon that the dihedral endpoint atlas sees.

## Tournament analogy

The tournament vertices are sections, not runners.

Pairwise observable: compare sections by

`(witnessed runner count, total open witness measure, private runners)`.

Switch: orient `a -> b` when section `a` has larger support vector; ties follow
cyclic section order, giving the Hamiltonian tie path.

The raw section tournaments in the scan were transitive, with one Hamiltonian
path and no directed cycles.  That is still information: section support is a
ledger, not yet the whole proof object.  The likely nontrivial tournament lives
one level down, on:

- Hall-deficit packets,
- wall-debt sections,
- boundary-crossing events,
- endpoint-mouth exchanges,
- proof obligations.

This mirrors the last LRC14 lesson: runner-level tournaments were too coarse,
endpoint-mouth orbit tournaments were better, and now section tournaments show
where the Hall debt is located.

## Proposed proof route

1. Define the compactified runner-section graph for an arbitrary speed set in a
   fixed observer gauge.
2. Prove the compactified section-checkoff lemma, or reduce its failure to the
   existing observer-source theorem.
3. Given a strict-open Hall packet `(R,S)`, attach its wall-debt sections and the
   local boundary events that support the compact matching.
4. Prove the wall-switch lemma: a deficient packet can cross a boundary to gain
   a section, unless the boundary packet is already an endpoint-mouth orbit.
5. In LRC14, identify those endpoint-mouth orbits with the dihedral mouth
   determinant from HYP-2569, so wall debt creates positive open interval mass.

The point is not to force every speed set into the beautiful "one runner per
quadrant" form.  The point is to replace the ugly cases by a small set of
checkable relations:

- Hall deficit,
- wall-debt support,
- boundary switch orientation,
- section-support tournament fingerprint,
- endpoint-mouth descent flag.

## Honest status

This is not a proof.  It is a reduction candidate with unusually clean finite
certificates.  The compactified matching result is only tested in bounded small
scans and named LRC14 rows.  The strict-open failure is systematic, and that is
the useful discovery: any real section proof must be a compactification plus
wall-switch proof, not a raw open-section matching proof.

Cross-links: HYP-2570, HYP-2024, S539, HYP-2486, HYP-2568, HYP-2569,
OPEN-Q-108, THM-523, T838.
