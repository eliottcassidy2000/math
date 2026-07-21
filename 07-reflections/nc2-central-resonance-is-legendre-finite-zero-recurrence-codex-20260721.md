# NC2 central resonance is finite zero recurrence, not merely no-common-consecutive-zero

**codex-2026-07-21.** Owner directive: work toward finishing NC2, pull often,
and mine seemingly unrelated work for ideas. Result: THM-2021 proves a new
arbitrary-degree central-resonance slice outside a countable exceptional set,
identifies the symmetric residual with finite recurrence of Legendre zeros, and
connects the higher-charge problem to THM-2020's finite-place channel separation.
Full NC2 remains open.

## The unrelated-work transfer that survived

THM-2005's support-Dirichlet atlas repeatedly warns that a scalar endpoint can
erase the profile carrying the theorem predicate. Applied here, the scalar
Gaussian moment is the endpoint of a polynomial in the return-channel marker.
Usually that profile is still coupled to the radial variable. On the exact
central hypersurface

```text
h = kappa b^r,
```

however, every channel has the same radial word `b^m`. The profile becomes rank
one and the moment factors exactly:

```text
E[P^m]=A_m(kappa)L(b^m).
```

That is the useful controlled quotient: keep the toral zero profile
`{m:A_m(kappa)=0}` and the radial EMP sequence separately. It converts signed
channel cancellation into a zero-recurrence question.

The second transfer came from spectral insertion/interlacing and the repo's
orthogonal-polynomial thread. The first instinct was the familiar
no-common-consecutive-zero descent. Challenging it exposed the load-bearing
logical gap: EMP needs `L(b^m)=0` eventually, while no consecutive toral zeros
only gives a syndetic subsequence. The correct invariant is **finite recurrence
of a fixed root**.

## The exact symmetric identification

For primitive charges `(+1,-1)`,

```text
A_m(kappa)=sum_k m!/(k!^2(m-2k)!) kappa^k
```

has OGF

```text
(1-2t+(1-4kappa)t^2)^(-1/2)
```

and therefore is a Legendre transform. Every zero parameter is negative real.
Thus the proportional central slice is already NC2-clear for every genuinely
complex phase and every nonnegative phase, with arbitrary radial `b`.

The remaining negative-real points are exactly nonzero recurrent roots of
Legendre polynomials. An official February 2026 IAS announcement for joint work
of Mangoubi--Kadets--Weller Weiser states the needed finiteness theorem. That
would close the full symmetric proportional hypersurface, but no public proof
was located, so THM-2021 labels this corollary conditional on the announced
input rather than laundering an abstract into a self-contained proof.

## What did not transfer

- The signed Cauchy--Schwarz/Parseval move from the far-comb work controls a
  signed scalar by a positive norm, but the relevant Bargmann norm does not
  vanish on one-sided NC2 elements. It remains orthogonal to the nullcone.
- Channel-by-channel dominance is unnecessary on the proportional hypersurface
  and false as a general philosophy after MISTAKE-211.
- The symmetric monomial boundary model is not the right difficulty test when
  `b` is a nonzero monomial: its first radial moment already fires. THM-2021
  keeps arbitrary signed radial `b`, where genuine moment cancellation exists.
- A three-term recurrence is an order statement, not a tail-zero theorem.
  Higher primitive charges have higher holonomic order (THM-1670), so copying
  the Legendre descent would challenge the wrong assumption.

## Incoming work synthesized

Concurrent THM-2018 independently found the symmetric proportional
factorization. Its reservation stub proposed consecutive-root descent; THM-2021
corrects that step and supplies the finite-recurrence criterion. Concurrent
THM-2019's affine-height closure is transverse: it handles arbitrarily many
charges when all balanced words have one Wick height, while THM-2021 allows an
arbitrary radial polynomial but requires rank-one channel proportionality.

Concurrent THM-2020 suggests the next attack on HYP-8771. Every exceptional
`kappa` is algebraic. At a finite place the channel valuation is factorial
carry-data plus `k v(kappa)`; a unique minimum forbids cancellation. The
higher-charge zero-recurrence problem can therefore be attacked with the new
finite-place certificate rather than another archimedean magnitude estimate.

## Tournament Analysis

Vertices considered: channels, roots, moment levels, and proof obligations.
Moment levels are the faithful carriers because tail recurrence is the
predicate. For levels `m<n`, the observable is `deg gcd(A_m,A_n)`; the baseline
edge `m->n` is flipped on a shared-root event, and chronological order is the tie
Hamiltonian path. Four exact scans (1,379 pairs total) have zero flips, zero
directed triangles, singleton SCCs, transitive score histograms, and one
Hamiltonian path. The quotient preserves exact pairwise recurrence but destroys
root identity, multiplicity, and analytic location, so it is evidence rather
than a proof of HYP-8771.

## Net

- **Proved:** exact all-charge factorization; NC2 off a countable exceptional
  set; in the symmetric case the exceptional set lies on the negative real
  axis.
- **Conditional on an announced 2026 theorem:** the entire symmetric
  proportional central hypersurface is NC2-clear.
- **Open:** higher-charge finite zero recurrence (HYP-8771) and the rest of the
  finite resonance band/general NC2.

Artifacts: THM-2021, HYP-8771,
`gmc2_proportional_legendre_finite_recurrence_thm2021.py/.out`. Cross-links:
THM-2017/2018/2019/2020, HYP-8766/8769, THM-1510/1660/1670/2005,
MISTAKE-211.
