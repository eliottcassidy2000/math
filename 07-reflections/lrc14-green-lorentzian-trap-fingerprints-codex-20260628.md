# LRC14 Green/Lorentzian trap fingerprints: Toeplitz as chart switch

The new thing HYP-3225 adds is not another contender scalar.  It is a trap
typology.  HYP-3224 already showed that every HYP-3202 arbitrary-exchange trap
has a Toeplitz lambda-min deficit.  The question was whether the Green-current
and Lorentzian languages would explain that deficit or fight it.  The answer
looks cleaner: they classify the ways local exchange got stuck, while Toeplitz
remains the boundary chart that discharges all of them.

The surprise is that the traps are visibly heterogeneous.  Six are best read
as rank-2 pair Plucker bottlenecks, two as low-connectivity Green networks,
two as AFM/frustrated Rayleigh debt, and one as a mixed sidecar.  That means
the finite trap list should not be proved by one local exchange monotonicity
lemma.  It should be proved by a dispatcher:

```text
exchange improves generic rows
Toeplitz detects the boundary
Green / Rayleigh / Plucker labels explain the boundary type
HB/Joukowski carries the odd Worpitzky sidecar
```

This fits the recent Erdos-870/tournament lessons better than a scalar proof.
The finite hard core should have a live core, deterministic filler, and canary
fields.  Here the live core is the HYP-3202 trap manifold; the deterministic
filler is the moment-cone chart; the canaries are effective resistance,
conditional Rayleigh debt, and pair Plucker gaps.  The canary fields do not
prove extremality alone, because AP itself has raw Rayleigh negatives and a
positive pair Plucker gap.  Their job is to keep the quotient legal by naming
which coordinate was lost when local exchange stalled.

The most useful philosophical correction is this: Toeplitz lambda-min is not
the whole explanation.  It is the chart switch.  It sees all trap classes,
while `lambda2_ratio`, Plucker gap, and Rayleigh-negative count each see only
part of the obstruction.  That is why the correlations are mixed instead of
near one.  The proof should use Toeplitz for discharge and the sidecars for
case structure.

I now think the LRC14 k=8 proof target has the following shape:

```text
1. AP support / covariance layers expose the normal cone off traps.
2. Ordered-tail exchange prices central q3 mass by q0+q6 loss.
3. Toeplitz / Fejer-Riesz / Verblunsky curvature discharges trap boundary.
4. Green-current and Lorentzian sidecars make the finite boundary theorem
   small enough to prove without enumerating anonymous rows.
5. HYP-3222/HYP-3212/HYP-3213 glue the finite certificate to the
   Chebyshev-Delsarte magic-function side.
```

The dangerous shortcut would be to say "prove all pair Plucker gaps are small"
or "prove all Green graphs are well connected."  AP does not satisfy such a
clean terminal rule.  The theorem must be relative to the Toeplitz/AP normal
cone and must retain orbit, odd-side, and exchange sidecars until they are
discharged.

The next concrete move is symbolic compression of the five trap classes.  If
each class can be represented by a small exact inequality family, HYP-3225
turns the `11`-row finite check into a finite theorem schema rather than a
table.  The wild but plausible connection is that the pair-Plucker class is
the valuated-matroid face of the same Chebyshev/Delsarte dual that Toeplitz
sees as a trigonometric moment cone.

After pulling HYP-3226, the immediate coordination rule is clear: its
small-pattern atlas should record these five trap classes as typed payload
atoms.  The numbers `11` and `12` are not the content.  The content is the
legal repair data attached to them: Toeplitz chart, Green bottleneck,
Rayleigh debt, and Plucker circuit.

-> HYP-3225, HYP-3224, HYP-3223, HYP-3222, HYP-3213, HYP-3212, HYP-3205,
HYP-3204, HYP-3203, HYP-3202, HYP-3163, HYP-3132, HYP-3226, T1308,
LTI-308, LTT-208, OPEN-Q-108.
