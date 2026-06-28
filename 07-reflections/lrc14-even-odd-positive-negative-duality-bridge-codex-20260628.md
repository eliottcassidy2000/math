# LRC14 Even/Odd Positive/Negative Duality Bridge

This pass reserves HYP-3238 as the bridge that the previous session was
pointing toward.  HYP-3236 showed that the positive covariance conductance
graph is a very strong AP-tight face.  The new synthesis says: yes, but that
face is lawful only as a packet with the odd/negative coordinates retained or
discharged.

After the checkpoint rebase, incoming HYP-3220 made the bridge sharper than
my initial wording.  Even/odd duality and positive/negative duality are not
just parallel shadows: in the de Moivre chart they are the same `Z/2`
parity/complement operator.  The verified power sums

```text
p1..p8 = -1,5,-4,13,-16,38,-57,117
```

have sign exactly `(-1)^k` because the dominant period is the negative Perron
root `-2cos(pi/7)`.  Complement sends sectors `(1,6),(2,5),(3,4)`, so the
positive-negative fold is also the even-odd parity operator.  That makes the
bridge less poetic and more algebraic: one sign bit is being seen in several
coordinate systems.

The next rebase brought HYP-3239, which is not a distraction; it is the same
object from two sharper angles.  The mac-mini S76 branch says the two proof
targets are one bimodal/phi4 extremality problem, with inclusion-exclusion
parity giving both the even/odd and positive/negative signs.  The kps S31av
branch says the topological sidecar is not generically Brouwer: for
`p=7=3 mod 4`, complement is a `D_7` anti-automorphism and the `Z/2` action is
free, so the right certificate is Borsuk-Ulam / odd degree.  That fits
HYP-3238 perfectly: the sidecar is a sign-isotypic packet, and the family law
chooses whether it is a fixed-point Brouwer/SOS packet or a free-antipodal
Borsuk-Ulam packet.

Then HYP-3241 resolves the apparent naming conflict: the Brouwer saddle and
the Borsuk-Ulam antipodal witness packet are two views of the same `Phi_14`
core witness set.  At n=14 those six unit witnesses form `3` antipodal pairs,
and that index `(p-1)/2=3` is also the de Moivre degree and the saddle
equioscillation count.  So the HYP-3238 sidecar should carry an index field,
not merely a sign field.

A later KPS witness check makes this more invariant than I had written it:
AP and Goddyn-Wong share the exact same six `Phi_14` witnesses.  The base
cyclotomic core is therefore a tight-locus packet, not an AP accident; the
harder cases are the dilation-promoted `Phi_{14d}` witnesses and rows where a
speed divisible by `14` destroys the base `t=a/14` certificate.

The crossed picture is:

```text
even / positive:
  Fejer square, SOS magnitude, pair-Pascal cap, covariance layers,
  positive Green conductance, Perron coherent mode, bulk measure.

odd / negative:
  Brouwer/Borsuk-Ulam topological sign, Worpitzky associator,
  Hermite-Biehler odd leg,
  negative covariance leakage, signed chart-change debt,
  odd-coordinate resurrection, antipodal split data,
  measure-zero cyclotomic core witnesses.
```

This merges HYP-3219 and HYP-3237 with the Green packet.  HYP-3219 says the
non-SOS cubic obstruction factors into sign times SOS magnitude.  That is the
cleanest possible warning against making the odd side into another SOS
problem.  HYP-3237 says the proof splits across a Vitali wall: bulk measure
works in the positive-measure region, while AP sits on a measure-zero core
where cyclotomic arithmetic replaces measure information.  The same pattern
appears when `C_E` is compressed to `G_+(E)`: positivity becomes usable only
after the clipped negative coordinate is kept somewhere.

The rebase added HYP-3220, which sharpens the bridge: even-odd duality is
positive-negative duality through the same `Z/2` complement/parity operator.
So the packet should now keep two explicit axes:

```text
even/odd           = 2-adic fold versus coordinate resurrection
positive/negative  = observer-blind gauge versus pair-visible cut/orientation
```

The rebase also added HYP-3239, which sharpens the topological slot.  For
`p = 1 mod 4`, the complement/reflection is orientation-preserving and the
fixed-point/Brouwer/SOS reading is the right side of the family.  For n=14
(`p = 7 = 3 mod 4`), the reflection is an anti-automorphism, the `Z/2` action
is free, and the sign side is Borsuk-Ulam odd degree on the `D_7`
sign-isotypic component.  So the bridge should carry the theorem type, not
just a scalar sign.

```text
topological_certificate_type = Brouwer_fixed_point | Borsuk_Ulam_free_Z2
D7_action_status             = aut | anti_aut
sign_isotypic_component      = retained | discharged | named_debt
```

The HYP-3234 signed-address chart tells us which full/even/odd recursion owns
the local letters.  HYP-3238 tells us which parity and sign coordinates survive
that chart change.

The information-theory translation is simple and sharp:

```text
compression is proof-grade iff the destroyed duality payload has zero
conditional entropy, is reconstructible, is dual-annihilated, or is named as
sidecar debt.
```

This extends the earlier law-defect idea beyond commutativity.  Associativity
matters because Schur complements and star-mesh edits compose only with
boundary terminals and eliminated variables still in the packet.  Positivity
matters because a conductance graph forgets negative covariance.  Evenness
matters because Fejer/SOS shadows forget the Worpitzky sign.  Measure matters
because AP's closed witnesses are invisible to open-measure mass.

The proof-frontier target I now want is a two-factor certificate:

```text
positive/even magnitude certificate
  plus
odd/negative sign or core sidecar
```

HYP-3222's Hermite-Biehler legs look like the local algebraic gluing theorem:
the even and odd polynomials interlace with positive Wronskian.  HYP-3204's
ordered-tail exchange looks like the finite pricing theorem: central odd
`q3` debt is paid by endpoint bimodality loss.  HYP-3236's Green graph then
becomes the electrical realization of the positive/even magnitude side.

Practical packet fields to add to the next bounded-bank audit:

```text
parity_axis_state
even_fold_certificate
odd_coordinate_resurrection
antipodal_split_status
sign_axis_state
sign_cut_status
topological_trace_or_degree_sign
D7_action_status
sign_isotypic_component
sos_magnitude_certificate
positive_conductance_packet
negative_leakage_sidecar
bulk_core_handoff
terminal_duality_discharge_or_debt
```

The next scout should not merely re-rank AP.  It should count the false
terminals: rows with no negative covariance leakage but still non-AP, rows
whose Green score is high but whose odd sidecar is live, and traps where the
Toeplitz/Green/ordered-tail discharges disagree.  It should also test whether
sign cuts commute with odd-coordinate resurrection; failures should be named
as antipodal split, Borsuk-Ulam odd-degree debt, signed-wall transport, Green
negative leakage, or Brouwer trace-sign debt.  That is where the proof packet
will either lock together or tell us exactly which sidecar remains unpriced.

The exact scout now does this over the `3432` anchored rows.  AP has

```text
q0 = 481/1470
q3 = 26/245
q6 = 1/49
q0+q6 = 73/210
L_y = 2633/7350
lambda2 = 0.192033074001
```

AP and the doubled AP are the only all-bank tight rows for `L_y`, `q0+q6`,
and `lambda2`; primitive normal form leaves AP unique.  But negative leakage
alone is a false terminal: there are `19` primitive zero-negative-edge rows,
so `18` primitive non-AP false terminals.  Connected positive graph is even
worse as a terminal quotient: `2754` primitive connected non-AP rows.

The ordered-tail bridge is cleaner than expected.  Among `2879` primitive
rows with positive `q3` debt, the exchange margin has `0` violations, and
the worst endpoint-bimodality gap per `q3` debt is `17161/12882`.  The `11`
non-AP HYP-3202 traps split into `8` with negative leakage plus odd `q3` debt
and `3` with odd `q3` debt but no negative covariance leakage.  Those three
are the best finite reminder that "negative" in the bridge is not merely
negative covariance; it is the full parity/sign sidecar.

So the next theorem-facing object is:

```text
Fejer / Green / cap magnitude
  plus
HYP-3220 parity sign or HYP-3222 HB interlacing
  plus
HYP-3204 q3 exchange-rate pricing.
```

That packet feels much closer to a proof shape than any scalar did.

-> HYP-3241, HYP-3240, HYP-3239, HYP-3238, HYP-3237, HYP-3236, HYP-3235, HYP-3234,
HYP-3233, HYP-3232, HYP-3231, HYP-3230, HYP-3228, HYP-3227, HYP-3225, HYP-3224, HYP-3223,
HYP-3222, HYP-3221, HYP-3220, HYP-3219, HYP-3218, HYP-3217, HYP-3216,
HYP-3214, HYP-3205, HYP-3204, HYP-3202, HYP-3201, HYP-3200, HYP-3163, HYP-3153,
HYP-3138, HYP-3004, HYP-2647, HYP-2637, THM-429, THM-426,
T1338, LTI-338, LTT-238, OPEN-Q-108.
