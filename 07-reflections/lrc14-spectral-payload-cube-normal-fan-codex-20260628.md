# LRC14 spectral payload cube: the AP face as a normal cone

The useful abstraction from the last few sessions is that the k=8 node is not
one extremal statistic.  It is a small diagram of compatible payloads.  The
old scalar question was "what does consecutive maximize?"  The better question
is "which proof currencies expose the same face, and which currencies are only
shadows?"

HYP-3224 gives a concrete answer on the bounded bank.  AP support, Toeplitz
`lambda_min`, all three covariance distance layers, and total covariance have
the same Pareto skyline: exactly the consecutive row and its nonprimitive
doubled AP dilation.  Raw cyclotomic energy fails again, with consecutive only
rank `19` for minimum energy.  Perron alignment is strongly supportive but
not a maximizer by itself, which is probably correct: it is a mode-alignment
sidecar, not the certificate.

The picture I now believe is:

```text
AP support        = exposed normal direction
Toeplitz lambda   = moment-cone curvature / interior slack
D1,D2,D3          = ferromagnetic covariance faces
Sigma_kappa2      = their sum, useful but lossy
Perron alignment  = leading-mode diagnostic
HB/Joukowski      = gluing theorem for even/odd dual legs
exchange traps    = local chart boundary, not counterexamples
```

The surprise is the trap discharge.  HYP-3202's `11` primitive
arbitrary-exchange traps looked like the price of using a local move theorem.
HYP-3224 shows that every one has a strict Toeplitz `lambda_min` deficit.  So
the traps are not random exceptional rows.  They are the boundary where local
exchange loses curvature, and the moment problem sees the curvature that
exchange forgot.

After the rebase, HYP-3204 slots into this picture cleanly.  It is the
coefficient face of the same cube: the central mass `q3` is allowed to grow
only when it pays a larger loss in `q0+q6`.  The scout now reproduces the
exchange-rate lemma with `0` violations and worst ratio `12882/17161`.  That
means the normal cone is not only spectral; it has a coefficient-pricing face:

```text
central Worpitzky mass q3
  <= ordered-state bimodality loss q0+q6
  <= AP normal-cone slack
```

The incoming HYP-3212/HYP-3213 work also clarifies what the dual certificate
should be.  The de Moivre cubic squared in `V_7(u)-2=(u-2)m(u)^2` is exactly
an equioscillation/magic-function signature.  So the object I called a
"Toeplitz dual" should probably be sought as a level-7 Chebyshev/Delsarte/
Cohn-Elkies magic function.  Toeplitz `lambda_min`, AP support, Perron mode,
and Hermite-Biehler interlacing are finite shadows of that function.

The newest HYP-3205 dictionary pass should be treated as the intersection
lemma layer underneath this: it verifies that the visible dictionary already
has AP/doubled-AP as the only simultaneous tight face.  HYP-3224 is the next
question: what dual object makes that dictionary intersection inevitable?

The incoming HYP-3223 green-current/Lorentzian packet is a natural
downstream test language for that question.  If the normal-fan picture is
right, effective-resistance profiles, Schur-complement trap exits, Lorentzian
Hessian signatures, and valuated-exchange slacks should not be new scalar
competitors.  They should classify the same finite trap boundary that
Toeplitz `lambda_min` already sees.  The best next computation is therefore a
joint trap table: first failed dictionary coordinate, Toeplitz slack,
Green-current bottleneck type, and valuated-exchange defect on the same
`11` rows.

That suggests a two-chart proof:

```text
Chart 1: exchange/covariance
  Every primitive row outside the finite trap manifold admits an improving
  exchange, or at least a descent in the AP support / layer profile.

Chart 2: moment cone
  Every trap is discharged by a Toeplitz, Schur, Verblunsky, Fejer-Riesz, or
  Christoffel dual certificate.
```

The theorem would not say "consecutive maximizes lambda_min" in isolation.  It
would say the moment-cone normal at consecutive supplies the missing boundary
chart for the covariance/exchange proof.  This is much closer to a proof
object than another ranked signal.

The rebase-added HYP-3222 packet makes the Joukowski/Hermite-Biehler merge
more concrete.  It gives exact minimal legs
`E=x^2+5x+4`, `O=x^2+4x+1`, strict interlacing, Wronskian
`(x+3)^2+2>0`, and an ideal C6 Perron quotient with
`lambda0=(1^T C 1)/6=lambda_max`.  That is the gluing layer.  The payload
cube provides finite shadows of the even leg
(biquadratic/Toeplitz/covariance/Perron) and the odd leg
(Worpitzky/interlacing debt).  The remaining hard lemma should identify the
slack variables in the payload cube with the even/odd interlacing slack.  If
that bridge is real, then tournament TRRT and the LRC cover bound are not
merely analogous; they are two coordinate charts on the same stability
certificate.

I would not chase "more cyclotomic" further without a direction.  Raw norm
has now failed in multiple ways.  The residue vector is useful only after
polarization by the AP residual or after embedding into the Toeplitz moment
cone.  Similarly, raw left-compression is useful only as a chart with named
trap exits.  The controlled-forgetting rule has become exact here: forget a
coordinate only if the next face of the cube can still see the AP normal.

The most promising new object to build is a dual certificate, not another
row statistic:

```text
find dual D such that
  slack_D(E) >= 0
  slack_D(AP) = slack_D(2*AP) = 0
  slack_D controls AP_support_gap
  slack_D controls Toeplitz lambda_min_gap on traps
  slack_D prices q3 gain by q0+q6 loss
  slack_D decomposes into D1,D2,D3 plus odd interlacing debt
```

This could be a PSD Toeplitz matrix, a Fejer-Riesz square, a Schur/Verblunsky
parameter inequality, a Delsarte-style linear programming witness, a
Chebyshev equioscillation certificate, or a Hermite-Biehler interlacing
determinant.  The rebase work pushes me toward the Chebyshev/Cohn-Elkies
formalism as the umbrella, with Toeplitz and covariance as finite sections.

Next I would test:

1. Whether the AP support gap is bounded by an explicit Toeplitz PSD slack.
2. Whether the `11` traps share a small number of Schur/Verblunsky failure
   words.
3. Whether HYP-3204's exchange-rate slack is a linear projection of the same
   magic-function dual.
4. Whether adding an odd Worpitzky/HB interlacing shadow preserves the
   two-row skyline.
5. Whether the same cube survives after leaving the bounded bank and passing
   to primitive normal forms.

The creative synthesis is "duodecimal" in the repo sense: twelve local maxima
in the exchange chart, two equality rows in the spectral chart, and a cube of
six visible metrics.  The number is not the proof.  The proof is the chart
change that turns the eleven decoys into moment-cone curvature.

-> HYP-3224, HYP-3223, HYP-3222, HYP-3221, HYP-3213, HYP-3212, HYP-3211, HYP-3210, HYP-3205,
HYP-3204, HYP-3203, HYP-3202, HYP-3201, HYP-3200, HYP-3163, HYP-3162,
HYP-3161, HYP-3160, HYP-3154, HYP-3153, HYP-3152, HYP-3151, HYP-3150,
HYP-3138, HYP-3132, T1307, LTI-307, LTT-207, OPEN-Q-108.
