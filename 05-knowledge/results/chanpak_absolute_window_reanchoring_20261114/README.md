# Moving absolute windows: free heights, two axes, and the carry that remains

**Partial constructive progress. Universal placement and Chan–Pak remain open.**
The session drew from [random archived notes](RANDOM_TRAILS.md) on knot
curvature, Rule 30, the planar Jacobian problem, Berggren calibration,
and other unrelated topics. The useful transfers were diagnostic:
take a mixed difference, keep a first carry, and distinguish a restricted
projection from the full physical data.

## What now reaches exact counts

- **[Four exact odds at every integer height](LINEAR_WINDOWS.md)**:
  908u+9g(u)+{752,754,756,758}, for every u>=1. No dyadic restriction
  on u. Literal changing cores and an exact inverse:
  [linear_window_replay.py](linear_window_replay.py).
- **[A two-axis changed-core return](TWO_AXES.md)**:
  both 336uv+141u+228v+101 and that count plus 2, with u,v independent
  powers of two at least 2. The actual cores have 9 and 33 paths,
  different pair-start/end counts, and different orders.
  [Exact selector and DP replay](two_axis_replay.py).
- **[Three interleaved switches](INTERLEAVED_CONTACT_RETURN.md)**:
  a wholly transitive puncture produces an exact two-count return
  on two infinite arithmetic lines of wing parameters.
  [Literal replay](interleaved_controls.py). The
  [complete classification](SWITCH_CLASSIFICATION.md) proves three
  is the maximum for independent disjoint two-site switches and
  corrects a preliminary statement: separated multi-bit patches
  interact, but one adjacent single-bit switch and a full carry can
  cancel. [Finite classification controls](orthog_probe.py).

## What prevents arbitrary placement

- **[Every complementary collar gap](COMPLEMENT_GAPS.md)**: the
  only three-point window in one slice is the original upper one.
  Even the union of all widths hits at most 2(r−1) of the 2^r−1
  odds just above the square of 2^r in the specified band.
  The new free-height bank hits under 4.41% asymptotically in these
  bands; the fixed dyadic banks cannot fill the rest.
  [Replay](complement_gaps.py).
- The [affine rigidity theorem](LINEAR_WINDOWS.md) shows that bounded
  alignment for all high addresses comes only from appending fixed
  binary digits. A finite adaptive family of non-dyadic affine
  reindexings also fails, with errors of order the height along a
  shared subsequence.
- The [exact ordinal carry test](ABSOLUTE_ODOMETER.md) separates
  CRT phase from the actual height. A binary neutral bank cannot
  make three consecutive odds. Full canonical residues alone still
  do not make a full absolute interval. But residue-wise redundant
  digits in all but one block **would** supply the missing carry
  and close a modulus-sized exact window, conditional on their
  physical construction. [Replay and literal controls](replay_odometer.py).

The clearest next proof target is to construct physical local digit sets
with the hidden quotient alternatives in the last result, or a changing
reader that supplies those alternatives without requiring complete local
banks. Neither a denser address family nor another congruence check by
itself resolves that obligation.
