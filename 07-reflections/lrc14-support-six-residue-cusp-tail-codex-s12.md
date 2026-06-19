# LRC14 Support-Six Residue-Cusp Tail

Codex S12, 2026-06-19.

The support-six residual changed shape today.  S10 said "execute a Minkowski
count."  HYP-2612 executed the first finite anti-coset count and showed that
span-only decay was false.  HYP-2613 measured the exact signed support-six
permanent layer and showed that the coupled absolute count was still too blunt.
S12 adds the missing coordinate:

```text
K(n_1,...,n_6,0,...,0) = C_d(n mod 7)/(n_1...n_6).
```

This is not just a computational convenience.  It changes the proof object from
a volume estimate to a residue-addressed reciprocal hyperplane sum.  The phrase
"Minkowski count" was hiding a fraction/address coordinate, and without it the
absolute envelope sees cusps that the signed kernel largely cancels.

The boundary-face ledgers make this visible.  The AP six-support shell through
`H=28` has absolute exact kernel mass about `0.920` but signed mass about
`0.0317`.  The resonant support `(1,2,3,4,5,21)` has `0.508` versus
`-0.00234`.  The sampled wide support `(2,3,4,5,6,68)` has `0.100` versus
`8.94e-5`.  On the last shell, one-coordinate boundary faces often carry most
of the relation count, but their signed totals are hundreds or thousands of
times smaller than their absolute totals.

The important correction is that the simplest possible cancellation is false.
The raw residue coefficient `C_d` does not have zero one-coordinate marginals.
So the proof cannot be a cheap "sum over a residue coordinate and vanish"
argument.  It has to know about the relation hyperplane `sum e_i n_i=0`, and it
has to delete or ledger the low-height anti-coset walls before applying a
summation-by-parts or cotangent/Dedekind-style estimate.

This gives a sharper statement of the remaining LRC(14) work:

```text
bounded finite certificate
+ finite low-height wall ledger
+ residue-cusp signed reciprocal tail
+ cluster translation quotient
< cap margin.
```

The tournament-analysis moral is also now clear.  The vertices should not be
runners or arcs in this layer.  They should be proof quotients: free envelope,
coupled absolute hyperplane, residue permanent, signed reciprocal shell,
boundary-face cancellation, wall ledger, and cluster quotient.  That quotient
destroys witness-time geometry, but it preserves the predicate that matters
here: whether the support-six correction can cross the cap margin.

LRC(14) is still open.  But the obstruction is narrower and less misty: prove
the residue-cusp tail theorem, and the previous finite machinery has something
precise to plug into.
