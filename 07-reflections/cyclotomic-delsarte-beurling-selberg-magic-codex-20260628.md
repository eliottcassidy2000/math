# Cyclotomic Delsarte/Beurling-Selberg magic: shell contact beats Fourier slogans

The useful surprise is that the magic function was already sitting in the
coefficient vector.  The k=8 functional

```text
10q0 + q3 + 10q6
```

is not just an arbitrary Delsarte row.  It is the shell polynomial

```text
f(n)=((n-1)(n-2)(n-4)(n-5))/4.
```

That makes the Beurling-Selberg analogy precise in the finite setting:
contact zeros at `1,2,4,5`, endpoint spikes at `0,6`, and one central
Worpitzky repair at `3`.  In centered coordinates `v=n-3`, it is
`(v^2-1)(v^2-4)/4`, the same quartic shell that keeps returning in the k=8
resolvent work.

The Delsarte face is equally clean:

```text
f(n)=10 -10 C(n,1) +10 C(n,2) -9 C(n,3) +6 C(n,4).
```

This recovers the HYP-3153 identity, but now with an interpretation: the
negative `S1,S3` coefficients are not a defect in bookkeeping; they are the
finite contact conditions that force the quartic to vanish on the four middle
nonterminal shells.  The central coefficient is the Worpitzky/odd repair that
HYP-3204 prices by loss of endpoint bimodality.

The cyclotomic face is also exact:

```text
P(z)=10+z^3+10z^6,
z^-3P(z)=10(u^3-3u)+1.
```

So the Laurent polynomial is just a shifted Chebyshev cubic under the
Joukowski map.  Reducing modulo `u^3+u^2-2u-1` gives `-10u^2-10u+11`.
This is the right kind of level-7 object: arithmetic in `Q(cos(2pi/7))`,
but still tied to the finite shell functional rather than a free-floating
root-locus slogan.

The important negative result is the PSD guardrail.  If we demand literal
nonnegativity of `10+rho z^3+10z^6` on the seventh roots, the required
central coefficient is about `18.019`, not `1`.  That means a naive
Fejer-positive or cyclic PSD magic function would destroy the tight
coefficient economics.  It would buy positivity by massively overpaying the
central shell.  This is exactly the sort of config-blind algebraic
overcertificate HYP-3221 warns against.

The rebase-added HYP-3214 makes this distinction cleaner.  The positive
sector-side cyclotomic magic function is the Fejer kernel `F_7`; this note is
about the shell `L_y` Delsarte dual.  The proof should glue those two faces
instead of forcing the shell coefficients to be the Fejer coefficients.

The later mac-mini S75 comb-overlap Gram packet is the measure-domain dual
face: a PSD kernel `K(p,q)=<1_Dp,1_Dq>` plus a single-arc peeling lemma.  That
is a better candidate for positivity than forcing the shell Laurent polynomial
itself to be PSD.

So the current proof target is not "make the magic function Fourier-positive."
It is:

```text
use the quartic shell contact as the Delsarte finite dual,
use AP support as the exposed root-locus normal,
use ordered-tail exchange to price the central q3 repair,
use Toeplitz/covariance layers to discharge trap charts,
and use Joukowski/Hermite-Biehler only as the gluing language.
```

The bounded-bank computation supports this split.  The magic functional has
no primitive beaters and only the doubled AP all-bank tie.  AP support has the
same equality set.  But the centered magic vector has cosine only `0.1827`
with the AP-support direction, so support is not the whole certificate.  It is
a coercive normal that must be coupled to the coefficient/exchange face.

This is a better invariant than "closest to cyclotomic."  Raw cyclotomic norm
failed.  Literal cyclic PSD also overpays.  The invariant that survives is a
typed packet:

```text
shell contact + Delsarte moments + Joukowski cubic
+ AP support normal + ordered-tail exchange sidecar
+ Toeplitz/covariance trap discharge.
```

The next computation I would run is a slack decomposition for each primitive
row:

```text
magic_deficit
= alpha * AP_support_deficit
  + beta * ordered_tail_exchange_slack
  + gamma * Toeplitz_or_layer_slack
  + residual
```

Then ask whether the residual has a small sign-controlled basis.  If yes,
HYP-3224's normal-fan dual becomes concrete: the shell polynomial is the
visible boundary equation, and the remaining coordinates are the sidecars
that make its slack provable.

-> HYP-3228, HYP-3214, HYP-3224, HYP-3204, HYP-3203, HYP-3153, HYP-3138, HYP-3132,
LTI-326, LTT-226, T1326, OPEN-Q-108.
