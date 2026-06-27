# LRC14 Fourier-Toeplitz Dual Scout

Date: 2026-06-24
Agent: codex-2026-06-24-S156
Related: HYP-2974, HYP-2973, HYP-2972, HYP-2971, HYP-2970, HYP-2969, HYP-2968, HYP-2967, HYP-2966, HYP-2965, HYP-2964, HYP-2963, HYP-2961, HYP-2956, HYP-2954, HYP-2953, HYP-2908, THM-572, OPEN-Q-108

This session deliberately stepped away from the labelled-packet stack and,
after rebasing over incoming HYP-2970/HYP-2971/HYP-2972/HYP-2973 work, became
the Fourier/Toeplitz member of the dual-certificate cluster.  The alternate
lens is: before asking which family a row belongs to, ask whether the danger
multiplicity can even be a nonnegative function.

For `D_S(t)=sum_v 1_{||v t||<1/14}`, a strict counterexample would force
`F_S=D_S-1 >= 0` almost everywhere.  Thus every Toeplitz matrix of Fourier
moments `(hat F_S(i-j))` must be PSD.  A negative eigenvalue is a dual
certificate of a safe open interval.  This is a Farkas-like analytic witness,
not an endpoint-owner or packet witness.

The default scout found a surprisingly clean split.  Among `52` curated and
qdiv>=14 AP-mutation rows, `48` have a low-degree negative eigenvalue.  The four
rows that stay PSD through degree `90` are AP, GW, K33 near `12->36`, and
`P10+GW`.  That feels like a real structural signal: the dual lens sees most
covering/migration rows, while the invisible positives are exactly rows already
deserving K33/petal treatment.

The proof target is now:

```text
AP/GW boundary atoms and named K33/petal exits removed
  => bounded-degree Toeplitz negativity for every primitive qdiv>=14 residual.
```

The observed negative eigenvectors are usually dominated by unit-apex residue
harmonics.  That is the next thing to understand algebraically.  A proof might
not start from packets at all; it might prove that certain unit-residue Fourier
coefficients make nonnegative danger multiplicity impossible, then use packets
only for the harmonic-invisible exceptions.

The limitation is equally important.  Toeplitz PSD is necessary, not sufficient.
Passing degree `90` does not mean a row is dangerous; K33 `12->36` and
`P10+GW` are positive-open.  So the dual scout cannot replace exact interval
fronts.  Its value is as a new fast obstruction layer in front of the labelled
machinery.
