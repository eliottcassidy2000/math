# Session143: magic-square and cyclic-structure bridges

**Historical checkpoint: 1 — five-vertex tournament correction and a genuine C5 operator bridge.**
Session began **2026-11-02 at 23:55 UTC** (user-local date November2);
source checks and this release continue on **2026-11-03 UTC**. The filename
suffix is the UTC release date, not a change to the session's starting date.

**The unrestricted 3-by-3 magic square of nine distinct integer squares remains OPEN.**
This checkpoint establishes no magic-square, lambda-cube, octonion, Hecke,
or sporadic-group representation equivalence. Subsequent lanes require their
own explicit maps and audits; a shared number is not such a map.

## What was established

- **PROVED / FINITE-EXACT:** the proposed order-five census is right:
  **12** tournament isomorphism classes, **8** self-converse, **10** after
  identifying converse. But there are **nine adjacency spectra**, and only
  **one regular class**. Two imaginary magnitudes belong to the same regular
  spectrum, not to two different regular tournaments.
- **PROVED:** for the regular tournament, `A=P+P^2`, where `P^5=I`.
  Its polynomial is `x^5-5x^2-5x-2=(x-2)q(x)`, with
  `q(x)=x^4+2x^3+4x^2+3x+1`. Its individual automorphism group is `C5`;
  `5:4` acts on a four-orientation family, not as one graph's automorphism
  group. A spectrum loses SCC position; an ordered-SCC sidecar repairs the
  complete order-five classification.
- **PROVED lemma + CITED character input:** in any complex representation,
  `B_g=rho(g)+rho(g)^2` for an order-five element has these same possible
  eigenvalues, with multiplicities obtained by Fourier inversion. For the
  Baby Monster's cited degree4371 character, the `5A` and `5B` operators have
  polynomials `(x-2)^891 q(x)^870` and `(x-2)^871 q(x)^875` respectively.
  These are not the representing matrix `rho(g)`, character values, or
  tournament adjacency matrices; their traces are42 and-8.
- **PROVED boundaries + CITED inputs:** the nontrivial eigenvalues generate
  `Q(zeta5)`, which does not embed over `Q` into `Q(zeta17)`.
  `Q(X0(13))` is a function field, not that number field. Conway's `M13`
  is a 13-point legal-move pseudogroup, not a five-point Mathieu group.

## Read and reproduce

- [Synthesis and explicit nonclaims](SYNTHESIS.md)
- [Concept board and connection contracts](BOARD.md)
- [Tournament census, proof, exact table and losses](tournaments/RESEARCH.md)
- [C5 operator lemma and cited Baby Monster specialization](sporadic/report.md)
- [Verification scope](VERIFICATION.md)

From any working directory, run `python3 /path/to/packet/replay_checkpoint1.py`.
The portable packet needs only the Python standard library. It replays from a
fresh relocation, uses an unrelated working directory, compares exact frozen
outputs, and tests positive and hostile cases under ordinary Python and `-O`.
The GAP script is supplied but **NOT EXECUTED**. Imported group-theory inputs
remain cited; portable arithmetic replay is not a reconstruction of their
character table. No priority claim or unrelated conjecture closure is made.
