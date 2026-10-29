# Independent replay of inherited attachment and phase controls

**Status: FINITE-EXACT / VERIFIED inherited controls.** These are not new
counterexamples, not the user's four large Hamilton-count witnesses, and not a
new Lonely Runner theorem. No claim here settles LRC(14) or a general memory
minimality problem.

## Provenance and implementation boundary

Source pin: `b07ec51b204e5c0e09e52eed206cb2bd3cc29ebc` in
`eliottcassidy2000/math`.
Source directory: `05-knowledge/results/joint_port_germs_20261028/`.
The inherited claims are in `README.md`, items 4 and 6; the literal tournament
inputs and operations are in `port_core.py`; the phase control is in
`check_lrc_germs.py` (the `phase_hostile` block).

`check_source_phase_controls.py` imports only Python standard-library modules,
not any producer module. It reconstructs adjacency matrices, uses a direct
subset Hamilton DP, cross-checks every reported Hamilton count by enumerating
all vertex permutations, and intersects closed intervals with exact rationals.
It is a bounded control replay, not a replay of the source's entire census.

## 1. Named attachment retains an observable distinction

Literal outgoing-neighbor bit rows, with vertices numbered 0 through 6:

- A: `(66,12,41,17,39,75,30)`;
- B: `(66,20,9,51,37,71,30)`.

Their native complement keeps consecutive path arcs `i -> i+1` and reverses
all other arcs. The checker verifies tournament validity, the distinguished
path, and complement involutivity on these inputs.

At **zero-based old port 1**, perform either:

1. append a vertex that beats only port 1 and loses to every other old vertex
   (the source calls this an exceptional sink); or
2. replace port 1 by a directed triangle, giving every triangle vertex the
   old port's external adjacencies.

| Operation after native complement | A | B | Resulting order |
|---|---:|---:|---:|
| Exceptional sink | 307 | 293 | 8 |
| Directed-triangle substitution | 777 | 981 | 9 |

All four direct DP counts agree with permutation enumeration. This verifies
these particular selected attachments; it does **not** reconstruct the
2,483,993 / 2,491,877 or 2,899,675 / 2,906,091 run-list witnesses.

## 2. The same local length does not determine continuation by a new speed

Put `G_V(d)={t in [0,1]: ||vt|| >= d for every v in V}`.
The implementation constructs this closed set by intersecting
`[(k+d)/v,(k+1-d)/v]`, for each `0 <= k < v`, retaining singleton components.
For `d=1/3`:

- `G_{1,3,6}(d)={4/9,5/9}`;
- `G_{2,3,6}(d)={2/9,7/9}`.

Write `e=1/3-d`. On the entire interval `0 <= e <= 1/72`, the components are

```
{1,3,6}: [4/9-e/3,4/9+e/6] union [5/9-e/6,5/9+e/3]
{2,3,6}: [2/9-e/6,2/9+e/3] union [7/9-e/3,7/9+e/6].
```

Thus both length germs are exactly `L(e)=e`, with `(L(0),L'(0+))=(0,1)`.
Adding speed 7 rejects both points of the first set (distance to an integer
is `1/9`) and retains both of the second (distance `4/9`). The exact interval
check also verifies that throughout `0 <= e <= 1/84` the first augmented set
is empty and the second retains precisely its displayed two components.

The whole-interval verification is not numerical sampling: all candidate
endpoints are affine rational functions of `e`. The checker enumerates every
pairwise crossing in the asserted interval, checks all boundaries, and checks
two interior points in each open ordering cell. Fixed endpoint ordering makes
each intersection endpoint affine on a cell; equality at two interior points
certifies its displayed affine formula there. There is one open cell in each
of the four cases. Singleton components are checked at the boundaries.

This replay certifies the uniform length germs and the selected continuation.
It does not claim to re-audit the source's general germ theorem or all its
directional-germ census. The controlled interpretation is that a phase-blind
length observation forgets information that a later speed constraint can test.

## Reproduction

```
python3 check_source_phase_controls.py
python3 -O check_source_phase_controls.py
```

Both produce identical JSON: **49 checks**, Hamilton counts
`[[307,777],[293,981]]`, common length germ `[0,1]`, and speed-7 survival
`[false,true]`. Checks use explicit exceptions, not optimization-removable
`assert` statements. No source files or remote state are modified.
