# Omega Extremes as the Cycle-Disjointness Axis

**Source:** opus-2026-05-29-S11  
**Computation:** `04-computation/omega_extreme_fingerprints_s11.py`

The H=63 correction looked at first like one more gap-table update: 63 is not
permanently forbidden, because two n=8 classes realize it. But the correction
has a sharper message. Both classes realize 63 in the least entangled possible
OCF shape:

```text
Omega(T) = K31, so H(T) = I(K31,2) = 1 + 2*31.
```

S11 found that the simplicity is even stronger. In each class, every directed
odd cycle contains one core vertex, and deleting that vertex leaves a transitive
tournament. The two signatures are `1001100` and `1100110`; both have the same
weighted inversion count 31.

So H=63 does not unlock by becoming complicated. It unlocks by finding a large
enough one-vertex odd-cycle star.

This is the opposite extreme from the THM-025 real-rootedness failure. There,
Omega is still extremely dense, but it has no cycle-family core. Instead it has
10 disjoint pairs and exactly one independent triple:

```text
(0,4,6), (1,3,7), (2,5,8).
```

That single triple is enough to break the root geometry. The polynomial
`1 + 94x + 10x^2 + x^3` is not failing because Omega is sparse; it fails because
the disjointness that remains is too concentrated.

This suggests a cleaner meta-axis:

```text
cycle-family core  <-->  dispersed disjointness  <-->  concentrated disjointness
```

Complete Omega with a core gives the H=63 unlock. Concentrated disjointness with
no core gives the real-rootedness failure. The ordinary H-spectrum, score
sequence, and even the number of cycles are projections of this deeper
geometry.

## Single-Core Signatures

If T-v is transitive, every odd cycle through v is determined by a `1...0` pair
in the insertion signature and an even choice of interior vertices:

```text
r_core(s) = sum_{i<j, s_i=1, s_j=0} 2^max(j-i-2, 0).
```

The target search up to signature length 16 found:

- `r=3` absent: the complete-core version of H=7 cannot occur.
- `r=10` absent: the complete-core version of H=21 cannot occur.
- `r=31` first appears at length 7, exactly where n=8 H=63 appears.

This may be the cleanest local model for why 63 is different from 7 and 21:
63 is the first forbidden-looking number that the weighted signature semigroup
can actually hit.

## Engineering Translation

For a Tournament TDA or ranking-fingerprint library, H alone is too compressed.
The useful features should include:

- cycle-family core size;
- alpha-vector of Omega;
- number and shape of disjoint cycle pairs;
- independent-triple supports;
- whether deleting a core leaves a transitive tournament;
- the single-core signature and `r_core(s)` when available.

These are small, inspectable features. They explain both a theorem-adjacent
phenomenon (H-gaps) and a numerical-analysis phenomenon (root failure). That is
exactly the kind of bridge this project keeps asking for.
