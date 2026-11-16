# What makes a shield reusable?

## The fixed-ratio guard has a real obstruction

Write L(w) for the vector of the two usual binary addresses and
R(w) for the vector of reversed complement addresses of a two-apex
word. The [mutual-arc contraction](apex/SHIELD.md) says that a weighted
cross term between adjacent words is det(L,R), up to a sign convention.
A family of atoms with L and R all on one common line therefore
cancels cross drift. The local endpoint term must still be checked.
That is what makes the p/q packets compose, whereas an arbitrary
signed difference of counts does not.

There is a limit on pumping a balanced atom by one constant guard.
For w_t=u c^t v, let c also denote its two-bit column and 1 the
two-coordinate all-one vector. Exactly,

    L(w_t) = 2^|v| [ 2^t(L(u)+c) − c ] + L(v),
    R(w_t) = 2^|u| [ 2^t(R(v)+1−c) − (1−c) ] + R(u).

If for two distinct positive t both vectors belong to one fixed
common line, subtracting shows that L(u)+c and R(v)+1−c lie on that
line. Their determinant is zero. By the guarded gap formula in
the same note, Δ(w_t) has the form K·2^t+C with K exactly this
determinant up to sign; thus Δ is constant for all t≥1.
A constant-letter dilation cannot produce new binary gap weights
while keeping the same exact composable response. The next candidate
is a periodically changing guard or a second seam that tracks the
cross determinant instead of forcing it to vanish.

The census also has an information ceiling. With B allowed at most
three wins, there are at most 2·2^n Σ_(j≤3) binom(n,j) words of
spine length n across both arc orientations, so across all lengths ≤n
the number of realizable counts is O(n³2^n). The rapid first-hole
growth from length 13 to 14 cannot continue at a fixed ratio above two
forever. This does not preclude eventual coverage by unbounded lengths.

## Compare what is selected, not just how many choices there are

- **Golden ratio, {2,3,11}.** The [odd-square bracket theorem](../procgen_numerology_20260923_odd_square_brackets.md)
  classifies primes p for which p and 2p fit in the same interval
  ((2h−1)²,(2h+1)²]. The interval has width 8h, whereas its lower
  endpoint grows quadratically: the ratio requirement forces a finite
  cutoff. The [golden carriers note](../collatz_golden_carriers_20261004.md)
  distinguishes that role of {2,3,11} from the halving counts of the
  known cycles and from golden denominators. The reader exponent
  1/5 comes from a different exact balance: at D≈N, the dyadic
  resolution S makes S² distinct phase pairs; the two errors are
  O(D^(1/4)S^(3/2)) and O(D^(1/2)/S), which balance at S≈D^(1/10).
  Golden and dyadic selectors share the need to keep an absolute
  lift in addition to a projective ratio. No occurrence of φ is
  needed for this exponent.

- **Fence and square packing.** The [fence and square-seam note](../inverse_fences_landing_20261116/FENCE_AND_SQUARE.md)
  charges actual lattice bonds to individual fences under its guards;
  a count of candidate contacts does not suffice unless the projection
  is injective or collisions are paid for. Its square-seam construction
  has exact landing label s−fu, hence multiplicity one when f>0.
  Compare the affine U landing of the [two-contact reader](READER.md):
  the phase-pair result alone counts addresses, but here the smooth
  derivative beats the entire discrete remainder, so every U gives
  a distinct baseline. Neither injectivity gives a prescribed landing.

- **8 and 9.** In the [archived 15-vertex chain and port repair](../coalescent_27_63_20261116/README.md),
  8 and 9 are forced ends. They have no common positive square-sum
  neighbor, because the two squares would differ by one. Three edges
  suffice by a new seam, 8–41–40–9 with colors 49,81,49.
  The lesson is a precise change in the governing equation: instead
  of a difference of two squares equal to one, solve
  2a²−b²=17 with valid distinct intermediate labels. Similarly,
  switching two contacts or passing through the second apex escapes
  the single-contact scale obstruction. Not an equivalence of the
  square-sum and tournament problems.

- **Catalan and Zsigmondy.** The [older 63/9 analysis](../collatz_mod6_20260917_zsigmondy_triad.md)
  shows exactly why the missing order 6 is a missing primitive prime:
  Φ_6(2)=3 already occurs at order 2, since 2³+1=9.
  It explains the 24 wheel; it does not forbid composite modulus 9 as
  a period-six clock. The [Catalan/cycle-gate note](../catalan_elliptic_20260921_synthesis.md)
  shows why a unit difference between pure powers is much narrower
  than an affine exact carry. The repaired tournament shell likewise
  gets +2 through cancellation of shared weighted slot profiles,
  not by making two large pure powers differ by one.
  Historical source: [Zsigmondy's 1892 paper](https://doi.org/10.1007/BF01692444).

- **The literal 6424 and the shorter named wheels.** 6424=8·11·73
  has prime coordinates of orders 10 and 9; their composite odd
  modulus 803 has order 90. That is distinct from asking for a
  *wheel whose ambient period is 6424*. The [wheel certificate](wheel/WHEEL_FOLLOWUP.md)
  rules out the latter even though abundant primitive factors exist.
  Orders 24,36,64,80,112 are covering periods, but their covering
  mechanisms differ: 24 uses the order-6 exception and its extra
  odd row; 36 uses a composite odd-part grid; 64 is dyadic-only;
  80 and 112 meet the single prime-part threshold for 5 and 7.
  The number 24 also happens to be the first supplied tournament
  order in the two-contact reader; there is no implication between
  those two occurrences.

A useful common principle: a projected family is best when its
coordinates can be **chosen from the target**, when the hidden
degrees of freedom have a native construction, and when the
landing equation comes with an explicit validity test. A label in
a dense projection is insufficient. The CRT slice, the physical
square-sum bridge, and the tournament arc contraction each keep
exactly the information used by their next operation; none alone
guarantees arbitrary Hamiltonian counts.
