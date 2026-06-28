# HYP-3222: Joukowski-Hermite-Biehler and Perron-Frobenius Synthesis

**Status:** SYNTHESIS / exact local scout; not an LRC14 proof.

## Claim

The current k=8 LRC14 proof frontier has two complementary spectral
certificates:

1. **Even / covariance route:** prove that consecutive speeds keep the
   all-ones vector in the Perron-Frobenius cone of the empty-sector covariance
   matrix.  HYP-3202's distance layers give a quotient where this is exact.
2. **Odd / Worpitzky route:** prove that the Joukowski image of the dip is a
   Hermite-Biehler packet: two real-rooted legs plus interlacing, with the
   self-inversive defect retained as a sidecar.

HYP-3222 is an exact-certificate addendum to the incoming HYP-3210 Joukowski
bridge merge and HYP-3211 apex-seven split.  It now also reads beside
HYP-3212's Chebyshev/equioscillation collapse, HYP-3213's
`Q(cos(2pi/7))` grand-synthesis frame, HYP-3221's apex-7 obstruction
guardrail, HYP-3205's spectral-dictionary compatibility layer, and
HYP-3204's ordered-tail exchange split.  It merges the S73c
Perron/Hermite-Biehler proposal, the S31al Toeplitz/Caratheodory route,
HYP-3202's covariance layer/trap scout, HYP-3201's law-defect compression
audit, and the HYP-3154/HYP-3162 Joukowski-cyclotomic bridge.

## Exact Local Certificates

The scout records a concrete Hermite-Biehler leg certificate at the minimal
k=8 fold.  The even biquadratic

```text
v^2 - 5v + 4
```

becomes the negative-axis leg

```text
E(x)=x^2+5x+4, roots -4,-1.
```

The odd Worpitzky/Eulerian leg is

```text
O(x)=A_3(x)=x^2+4x+1, roots -2 +/- sqrt(3).
```

These roots strictly interlace:

```text
-4 < -2-sqrt(3) < -1 < -2+sqrt(3).
```

The Wronskian has the right orientation globally:

```text
E O' - E' O = x^2+6x+11 = (x+3)^2+2 > 0.
```

So the minimal even and odd legs already form a strict Hermite-Biehler pair.
The remaining proof debt is the lift from this leg-level certificate to the
full LRC miss-PGF through the Joukowski map, where
`q_t R^t != q_{6-t} R^{6-t}` and off-circle `Im(w)` must be retained.

The scout also turns HYP-3202's three consecutive covariance layer sums into
an idealized C6 Perron quotient.  With

```text
D1 = 308509/1080450
D2 = 547577/2160900
D3 = 225577/1234800
```

and layer weights `w1=D1/6`, `w2=D2/6`, `w3=D3/3`, the nonnegative circulant
first row `[0,w1,w2,w3,w2,w1]` has eigenvalues

```text
lambda0 = 6237419/25930800
lambda1 = -1440157/25930800
lambda2 = -750151/25930800
lambda3 = -1856803/25930800
lambda4 = -750151/25930800
lambda5 = -1440157/25930800
```

and

```text
(1^T C 1)/6 = lambda0 = lambda_max.
```

Thus the distance-layer quotient has zero Perron-alignment defect.  A signed
AFM contrast with first row `[0,-w1,+w2,-w3,+w2,-w1]` moves the top mode to
`k=3` and makes the all-ones mode negative, so coherence/sign pattern is a
necessary coordinate, not just total layer mass.

## Proof Program

The even side should now be phrased as a boundary-aware Perron problem:

```text
prove PF alignment for the actual empty-sector covariance matrix,
or prove that every competitor loses in D1/D2/D3 before boundary error can help.
```

The Toeplitz/Caratheodory route is the global moment analogue: `lambda_min(T)`
measures how interior the moment packet is.  The Perron route measures whether
the covariance quadratic form is carried by the coherent positive mode.  The
distance-layer route is the finite exact quotient where that coherence becomes
row-sum regularity.

The odd side should now be phrased as a Joukowski-Hermite-Biehler problem:

```text
prove interlacing after transporting the even biquadratic leg and the
Eulerian/Worpitzky leg to the real axis, while measuring self-inversive defect.
```

The exact `A_3` interlacing above is not the full theorem, but it is a real
anchor for the route.

The post-fetch constraint is important: HYP-3212/HYP-3213 say the Joukowski
side is not only real-rootedness folklore but the arithmetic of the
seventh-cyclotomic cubic and its Chebyshev equioscillation; HYP-3221 says a
configuration-blind Lee-Yang/Asano/SOS certificate must fail at the apex-7
obstruction and retracts octonion/Fano structure as load-bearing LRC algebra;
HYP-3205 says Perron alignment is one coordinate in a certificate dictionary,
not a terminal scalar; and HYP-3204 says the theorem-facing coefficient target
is the ordered-tail exchange-rate split, not raw `q3` or full convex order.
HYP-3222 therefore treats HB/PF as a packet to be joined to those
analytic-equidistribution, certificate-dictionary, and coefficient-exchange
sidecars, not as a standalone scalar compression.

## Compression Guardrails

Several tempting compressions fail in different ways:

- `radius_only` preserves `q0=q6*R^6` but forgets root angles, `Im(w)`, and
  self-inversive defect.
- `total_covariance_only` preserves `Sigma kappa_2` but forgets `D1/D2/D3`
  and Perron alignment.
- `commutative_pair_mass` preserves the Pascal/binomial cap but forgets the
  odd Worpitzky associator and bracket/order sidecar.
- `positive_association` preserves nonnegative pair covariances but forgets
  their location and the finite trap identity; HYP-3202 has `19` primitive
  all-nonnegative decoys.
- `hb_leg_interlacing` preserves the exact even/odd leg relation but still
  forgets the full miss-PGF self-inversive defect unless that sidecar is
  carried.

This is the information-theoretic moral of HYP-3201 in spectral language:
a quotient is legal only when the target function has zero residual defect on
the fibers, or the missing coordinate is named and retained.

## Tournament Analysis

Vertices are proof certificates and obligations, not runners, arcs, or raw
sectors.  The pairwise observable is which carrier preserves more LRC14 proof
payload toward the k=8 node.  The scout uses a route-score gauge and reports a
transitive tournament:

```text
score_hist={0:1,1:1,2:1,3:1,4:1,5:1,6:1,7:1,8:1,9:1}
directed_3cycles=0
hamiltonian_path_count=1
priority_path=
  hermite_biehler_interlacing_certificate
  -> perron_alignment_certificate
  -> toeplitz_caratheodory_margin
  -> distance_layer_pf_quotient
  -> self_inversive_defect_sidecar
  -> random_current_coupling_order
  -> law_defect_entropy_meter
  -> raw_total_covariance_scalar
  -> plain_positive_association
  -> row_entropy_scalar
```

Assumption challenge: the useful vertices here are not runners.  They are
root-polynomial legs, covariance layers, spectral margins, self-inversive
defects, random-current order data, and proof obligations.  The quotient
preserves the k=8 covariance/dip proof predicate; it destroys boundary/apex
matrix data and full miss-PGF self-inversive information unless those are
retained.

## Proof-Frontier Use

Next proof obligations:

1. Lift the exact `E/O` interlacing certificate through the Joukowski map while
   measuring self-inversive defect.
2. Replace the ideal C6 Perron quotient by a boundary-aware covariance matrix
   inequality.
3. Join Toeplitz `lambda_min`, Perron alignment, distance layers,
   spectral-dictionary compatibility, and random-current order into one
   packet.

## Artifacts

- Script: `04-computation/lrc14_joukowski_hb_perron_synthesis_codex_20260628.py`
- Result: `05-knowledge/results/lrc14_joukowski_hb_perron_synthesis_codex_20260628.out`

## Links

HYP-3222 extends HYP-3221, HYP-3213, HYP-3212, HYP-3211, HYP-3210,
HYP-3205, HYP-3204, HYP-3202, HYP-3201, HYP-3200, HYP-3163, HYP-3162,
HYP-3161, HYP-3160, HYP-3154, HYP-3153, HYP-3152, HYP-3151, HYP-3150,
HYP-3147, HYP-3144, HYP-3142, HYP-3139, HYP-3138, HYP-3132, HYP-3122,
THM-577, LTI-306, LTT-206, T1306, and OPEN-Q-108.
