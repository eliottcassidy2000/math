# LRC14 Green Conductance And Algebraic Connectivity

This pass executed the electrical side of HYP-3223.  The result is stronger
than I expected: after turning positive empty-sector covariances into a
conductance graph, AP/consecutive and its doubled dilation are the only rows
maximizing algebraic connectivity and total positive conductance, and the
only rows minimizing Kirchhoff, mean, max, and distance-layer effective
resistance.  Primitive normal form leaves AP alone.

The clean conceptual object is the function chain

```text
E -> C_E -> G_+(E) -> L_E -> L_E^+ -> resistance/lambda2/current entropy.
```

That chain is a set of compressions.  The important part is not "compression
is bad"; the important part is knowing which compression has a legal sidecar.
Here the dangerous arrow is `C_E -> G_+(E)`, because it clips negative
covariance.  The scout keeps `negative_edges` and `negative_mass`; without
them, the conductance graph would falsely look like a pure positive
association theorem.  This is exactly the HYP-3201 law-defect lesson in a new
body: a quotient is allowed only when the lost coordinate is zero, recovered,
dual-annihilated, or explicitly charged.

The new coordinate has a real proof smell.  AP's Laplacian spectral gap is
`lambda2=0.192033074001`; its Kirchhoff index is `108.654718079151`; its
largest effective resistance is `9.713313375596`.  The closest non-AP Green
decoy already drops to `lambda2=0.144321509290`, raises Kirchhoff to
`136.424340938360`, and carries a negative-leakage edge.  Some decoys have no
negative edges, which is useful: it means leakage alone is not the theorem.
Connectivity and resistance are doing independent work.

The trap table is the most proof-facing output.  HYP-3202 had `11` non-AP
arbitrary-exchange local maxima.  HYP-3224 discharged all of them by strict
Toeplitz lambda-min deficit.  HYP-3236 discharges all of them by Green
resistance excess: three primarily by Kirchhoff excess and eight primarily by
max-resistance bottleneck.  That is the same finite boundary seen through a
different chart.

The right mental model is now:

```text
exchange chart:
  move until a finite trap manifold appears

moment chart:
  traps have Toeplitz / Schur / Verblunsky curvature deficit

Green chart:
  traps have effective-resistance bottlenecks or leakage

odd chart:
  Worpitzky / Hermite-Biehler sidecar carries the associator debt
```

This also clarifies how to read Perron-Frobenius.  The covariance-side Perron
story says AP aligns the all-ones mode with the positive coherent phase.  The
conductance-side Laplacian story says AP maximizes the spectral gap away from
the constants.  Those are not duplicate metrics; they are dual views of the
same current packet.  The covariance matrix wants high all-ones energy; the
conductance graph wants no cheap orthogonal bottleneck.

The incoming HYP-3214 Fejer result is the harmonic sibling of this.  It says
the 7-sector magic function is `F_7=(de Moivre cubic)^2`, positive-definite
with weights `(7-|n|)_+`, and equal to AP autocorrelation.  HYP-3236 says the
finite covariance shadow of the same AP packet becomes a best-connected
Dirichlet network.  The catch is useful: HYP-3214 separates the 7-sector Fejer
magic function from the 14-clock Johnson pair-Pascal cap.  So the Green graph
must be glued to Fejer and cap through a packet diagram, not by identifying
their scalar shadows.

The new HYP-3231 scale-recursion ledger, HYP-3232 interlocking-recursion audit,
and HYP-3216 family-law updates raise the bar again.  HYP-3230 gives the
three-gap/Stern-Brocot cap-kernel recursion; HYP-3232 says the
Mobius/Eisenstein/Legendre modes concentrate at the apex half `n/2=7`; and
HYP-3216 descends the LRC(2p) route through a moment-order ladder with a
2-adic fold.  A Green proof cannot stop at "AP has largest `lambda2`."  The
conductance packet has to be indexed by enough gap/cut data to see the same
scale-normal address, or it must prove that the missing address is dual
annihilated.  That is the information-theoretic version of the current
picture: Green resistance is a compression channel, and scale recursion is a
test for whether the channel kept the proof-relevant bits.

HYP-3217 adds another bit that is too structured to throw away: the cubic
Gaussian-period mode gives the de Moivre angles inside the subfield lattice of
`Q(zeta_7)`.  So the Green packet should not only ask whether conductance cuts
are scale-normal; it should ask whether Fiedler vectors and bottleneck currents
live in, avoid, or mix the cubic cosets `{1,6}`, `{2,5}`, `{3,4}`.  That is a
clean next measurement because it turns the abstract "Perron/HB/Joukowski"
sidecar into a concrete modal projection of the current graph.

HYP-3233 sharpens that modal projection into a polynomial grading: recursion
modes look like `(x-1)^depth * Phi_d`, so the Green-side spectral data should
ask which cyclotomic factor the conductance compression preserves.  Algebraic
connectivity is one eigenvalue; a proof packet needs the mode label too,
especially whether the `Phi_7`/de Moivre cubic hardness is present, killed, or
dual-annihilated by the Green chart.

The signed-address HYP-3234 sheaf says the same thing operationally: the
recursion formulas are local charts with slot bases and cancellation debt.  A
Schur or star-mesh Green reduction that forgets which chart it used can look
associative after scalarization while losing the actual proof transport map.

HYP-3235 and HYP-3218 are the sharpest new pressure on this Green line.  If
Fejer is now an explicit cyclotomic positive certificate, then Green resistance
has to stop being merely a beautiful diagnostic and explain its relation to
that square: factor through it, provide a Thomson dual to it, or mark the
remaining electrical residual.

Joukowski/Hermite-Biehler still has to enter because the conductance graph is
an even shadow.  It sees pair covariance and Green resistance, but it does not
see the odd Worpitzky/associator piece by itself.  The HYP-3222 exact legs
`E=x^2+5x+4`, `O=x^2+4x+1`, strict interlacing, and positive Wronskian are
still the gluing theorem candidate.  HYP-3236 is a new finite face that the
eventual magic-function dual should expose.

There is a subtle associativity/compression warning.  Schur complements and
star-mesh transformations compose cleanly only when the terminal set and
eliminated variables are still present.  If one collapses the row to
`lambda2` too early, the edit algebra loses its bracket data.  That is the
same failure pattern as subtraction or exponentiation under an associativity
quotient: the scalar shadow no longer knows which operation order produced
it.  So the proof object is not `lambda2` alone; it is a packet containing
conductance graph, Green kernel, leakage sidecar, current path, and trap edit.

The candidate theorem I would try next is a Rayleigh/Thomson separation:

```text
Every primitive non-AP k=8 row has either
  lower positive-conductance algebraic connectivity,
  larger Green resistance for a named demand,
  or a negative/odd sidecar routed to another chart.
```

This could be attacked by cut inequalities, by a Poincare/Cheeger lower-bound
comparison, or by finite Schur-complement words on the `11` traps.  The best
version would not prove the Green coordinate in isolation.  It would show that
Green slack is one projection of the HYP-3224 normal cone, alongside AP
support, Toeplitz lambda-min, covariance layers, and ordered-tail exchange.

The assumption challenge matters here.  Runners, gaps, covariance edges,
conductance graphs, current paths, Fiedler modes, Schur edits, trap rows, and
proof obligations are all plausible vertices.  For this pass, conductance
certificates win because they preserve AP extremality, distance-layer payload,
trap identity, effective-resistance bottlenecks, and algebraic-connectivity
margins.  They destroy raw runner identities and negative signs unless the
leakage sidecar is carried.  So the row statistic is not a proof; the
sidecar-aware Green packet might be.

Next tests:

1. Emit Schur-complement reduction words for the `11` non-AP traps and see
   whether the max-resistance versus Kirchhoff split corresponds to two or
   three network types.
2. Compare Green slack against HYP-3224 Toeplitz slack and AP-support slack on
   every row, not only on traps.
3. Try a cut-family certificate for `lambda2`: does AP dominate all relevant
   Fiedler cuts after leakage is named?
4. Attach Lorentzian/valuated-exchange data to the same trap rows; test
   whether Schur-complement discharge and tropical Plucker defect are the same
   finite event.
5. Search for the magic-function dual whose finite shadows include
   Toeplitz lambda-min, AP support, and this Green conductance gap.

-> HYP-3236, HYP-3235, HYP-3234, HYP-3233, HYP-3232, HYP-3231, HYP-3230, HYP-3229, HYP-3228, HYP-3227, HYP-3226, HYP-3225, HYP-3224, HYP-3223, HYP-3222, HYP-3221, HYP-3218, HYP-3217, HYP-3216, HYP-3214, HYP-3213, HYP-3212,
HYP-3211, HYP-3210, HYP-3205, HYP-3204, HYP-3203, HYP-3202, HYP-3201,
HYP-3200, HYP-3163, HYP-3162, HYP-3161, HYP-3160, HYP-3154, HYP-3153,
T1336, LTI-336, LTT-236, OPEN-Q-108.
