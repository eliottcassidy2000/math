# Joukowski-Hermite-Biehler and Perron-Frobenius Synthesis

The new shape is that both halves of the k=8 node are asking the same
question in different languages:

```text
When is a compression honest?
```

The Perron side says that compressing a covariance matrix to `1^T C 1` is
honest only if the all-ones vector is actually the positive top mode.  HYP-3202
already found a finite exact shadow of this: split consecutive's covariance
into distance layers and the ideal C6 quotient has all row sums equal, hence
`1` is exactly Perron.  But the actual empty-sector matrix has apex/boundary
deviation, so the quotient needs a boundary-aware sidecar.

The Hermite-Biehler side says that compressing the odd/even dip to root data
is honest only if the two real-rooted legs interlace, and only if the
Joukowski/self-inversive defect is carried.  The small new certificate is
surprisingly clean:

```text
E(x)=x^2+5x+4, roots -4,-1
O(x)=A_3(x)=x^2+4x+1, roots -2 +/- sqrt(3)
E O' - E' O = (x+3)^2+2 > 0.
```

So the minimal even biquadratic and odd Eulerian/Worpitzky legs are already a
strict HB pair.  That does not prove LRC14; it tells us exactly what remains:
lift this interlacing through the Joukowski circle map while measuring the
off-circle and non-self-inversive error.

This also clarifies the older compression story.  Commutativity compresses
pair mass into Pascal/binomial cap data, but it cannot see the Worpitzky dip.
Associativity-like scalarization smells like `1/7`, but HYP-3200 killed that
as a theorem.  Positive association sees the sign of pair covariance, but
HYP-3202 found `19` primitive rows with all pair covariances nonnegative.  Row
entropy is not the direction either.  The missing invariant is not a single
smaller number; it is a spectral packet:

```text
Toeplitz lambda_min
Perron alignment of 1
D1/D2/D3 covariance layers
Joukowski self-inversive defect
Hermite-Biehler interlacing
odd Worpitzky/ear sidecar
law-defect entropy meter
```

The best next proof attempt should not choose between Perron and
Hermite-Biehler.  Perron is the even covariance certificate; HB is the odd
stability certificate.  Toeplitz sits above both as the moment/stability
language, and HYP-3201 tells us when a shortcut has illegally forgotten a
coordinate.

The latest mainline fetch sharpens the packet.  HYP-3212/HYP-3213 say the
Joukowski leg is really the Chebyshev and `Q(cos(2pi/7))` arithmetic story;
HYP-3221 warns that configuration-blind algebraic positivity hits the apex-7
obstruction and that octonion/Fano structure is numerological here, not a
load-bearing LRC identity; HYP-3205 folds Perron into a certificate dictionary
and demotes it from terminal scalar to diagnostic coordinate; HYP-3204 prices
the odd central mass against `q0+q6` through the ordered-tail exchange-rate
lemma.  That makes the proof target feel less like "find the right scalar"
and more like "prove the scalar only after its sidecars have zero residual
defect."

The tournament for this session therefore uses certificates as vertices, not
runners.  The live Hamiltonian path begins with Hermite-Biehler interlacing
and Perron alignment, then Toeplitz margin and distance-layer quotient, and
ends with raw covariance, plain positive association, and row entropy.  That
ordering feels right: it keeps the proof predicate alive instead of rewarding
the most compact slogan.
