# Oddness, signed spectra, and operation-preserving correspondences

**Status: PROVED scoped deductions / CITED inherited mechanisms / FINITE-EXACT controls.
Same-agent proof audit and alternate computational paths; no independent referee.**
Session: 2026-10-28. Proposition labels below are local to this packet, not reserved
canon IDs. LRC(14), a general unit-distance extremum, and a universally smallest
attachment carrier remain OPEN.

The useful common object is a **representation together with its symmetry action
and the operations it must preserve**. Equal counts are its weakest shadow.
For a tournament completion cube, a signed support polynomial simultaneously records
class counts, complement-fixed classes, and the spectrum of single-edge flips.
It extends the repo's existing Walsh interpretation of sign-even graphs to arbitrary
marked oriented carriers. Parallel-path carriers make the distinction decisive.

## Results

| Local result | Exact consequence | Boundary |
|---|---|---|
| S1, signed support polynomial | F(1) counts completion classes; F(-1) counts free-complement fixed classes; coefficients give weighted flip multiplicities | Relative automorphism action; weighted cube quotient, not an unweighted Hamilton metagraph |
| S2, bigon obstruction | Two equal directed paths plus an isolate have weak trijections, but no smaller tournament action has the same flip spectrum | New operation certificate for an already known non-monomial family |
| P1-P3, paired oddness | A free Klein four action explains 4-divisibility of odd self-converse classes of orders at least 5; joint negative signs can cancel | Uses the inherited Eplett class correspondence; weaker scalar result than existing 2-adic formulas |
| L1-L2, LRC reflection | The trace on good-time components is 1 exactly when all velocities are odd | Mixed parity can have empty, isolated, or interval witnesses with the same zero trace |
| G1-G3, geometric carriers | Faithful planar unit-distance realization and marked trijectivity imply neither one another | Explicit witnesses; no extremal unit-distance claim |
| G4, geometric encoding capacity | A faithful planar graph-only encoding needs m log2(m) >= (1/4-o(1)) n^2 vertices-times-log-vertices for n-vertex tournament classes | Exact real coordinates, labels, and unbounded sidecars are excluded from this output contract |

The smallest bigon has 32 marked completion classes and signed spectrum
**(1,4,8,10,7,2,0)**. Its free complement fixes no class; the natural edge-set
tournament target has eight fixed classes. For path length k the latter number is
2^(k^2-1), while the source remains zero.

Complete aligned-pair calculations:

| Core n | Odd self-converse order | Classes | Joint sign modes lost by demanding separate positivity | Orbits of all three flips, sizes 4 / 8 |
|---:|---:|---:|---:|---:|
| 2 | 5 | 8 | 2 | 2 / 0 |
| 3 | 7 | 88 | 16 | 8 / 7 |
| 4 | 9 | 2752 | 300 | 64 / 312 |

Thus the occurrence of **15** at core three is precisely 15 orbits of a specified
three-involution action. No Fano, Hamming, or vector-space structure on those 15
orbits has been established.

## Inheritance and portfolio

- **Anchor:** operation-preserving tournament/sign-even-graph/trijection comparison.
  Closest proved mechanisms are THM-4721
  [signed Walsh duality](../../../01-canon/theorems/THM-4721-twisted-burnside-on-odd-fixed-point-sets-sylow-refinement-of-even-graphs-and-tournaments-components-walsh-basis.md)
  and THM-4835
  [fixed-subgraph Euler criterion](../../../01-canon/theorems/THM-4835-euler-characteristic-criterion-for-marked-excess-minus-one-carriers-in-tournaments.md).
- **Niche:** odd self-converse classes via the
  [Session89 aligned-pair lift](../session89_consecutive_tournament_pairs_20261028/encoding/ODD_SELF_CONVERSE_LIFT.md).
  The least-used sidecar here is the **joint stabilizer sign**: discarding it
  loses modes even when the two marginal sign tests look natural.
- **Wildcard / LRC anchor check:** reflection on exact good-time components, then
  unit-distance realizations of the same parallel-path carriers.
- **Canonical hostile:** the 12-vertex equal-cyclic-triples pair with automorphism
  orders 81/27; the inherited guarded reconstruction remains necessary.
  Bigons and odd-prime theta carriers provide smaller operation-specific hostiles.
- **Corrected near miss:** an odd-order fixed-subgraph defect is not sufficient
  for failure of weak trijectivity (THM-4835 audit / MISTAKE-689).
  The current two-endpoint singleton carrier also stops reconstructing at seven;
  no conclusion here upgrades its weaker observation language.

The live concept board was: **signed stabilizers; E=V-1 fixed-subgraph marks;
paired converse actions; good-time reflection; geometric realizability**.
The first pull unified counts and flips but exposed the relative/global sign distinction.
The second found sign cancellation on the paired side and free four-element orbits.
The third showed that LRC's natural reflection has no parity leverage on mixed speeds.
The final pull made geometric realization independent of trijectivity and connected
both sparse and geometric encodings through a quantified information bound.

## Reading and reproduction

- [THEORY.md](THEORY.md): S1-S2 and P1-P3, with all-size proofs.
- [LRC_UNIT_DISTANCE.md](LRC_UNIT_DISTANCE.md): L1-L2 and G1-G4.
- [HISTORY.md](HISTORY.md): inherited versus new claims, repaired historical language,
  typed connections to arithmetic codes and attachment carriers.
- [AUDIT.md](AUDIT.md): universes, controls, proof review, limitations.
- [reproduce.py](reproduce.py): runs normal and optimized Python, compares committed
  receipts, and rechecks the inherited 81/27 hostile.

From the repository root:

    python3 -B 05-knowledge/results/odd_sign_spectra_trijections_20261028/reproduce.py

Only the geometry checker needs SymPy (tested with 1.14.0). Other checks use the
standard library. The receipts are deterministic mathematical outputs; floating
distance gaps are diagnostics and never decide an equality.

## Primary literature

The original equinumerosity is
[Royle, Praeger, Glasby, Freedman and Devillers, *Tournaments and Even Graphs are Equinumerous*](https://arxiv.org/abs/2204.01947).
Here “even” concerns automorphism orientation signs, not even vertex degrees.
The odd self-converse lift is inherited from
[Eplett, *Self-converse tournaments*](https://doi.org/10.4153/CMB-1979-004-6)
and the repo's explicit Session89 codec.
For the distinction between faithful and nonfaithful unit-distance graphs and
their counting bounds, see
[Alon and Kupavskii, *Two notions of unit distance graphs*](https://www.tau.ac.il/~nogaa/PDFS/ramsey10.pdf).
No priority claim is made for standard representation theory, reflection pairing,
two-circle geometry, or these cited correspondences.

