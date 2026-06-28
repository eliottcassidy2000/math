# LRC14 Modular Magic Sidecar Audit

This pass clarifies what the modular/arithmetic prompt can contribute to the
proof frontier.

The durable core is still HYP-3214: the sector magic function is the Fejer
kernel `F_7`, equivalently the squared de-Moivre cubic in the Chebyshev
identity

```text
V_7(u) - 2 = (u-2)(u^3+u^2-2u-1)^2.
```

That is the actual Cohn-Elkies-shaped object.  It has the double zeros,
positive Fourier coefficients, and exact `7^2` discriminant.  The Johnson
`J(14,2)` cap kernel remains the separate 14-clock cap object.  The proof
target should therefore be two-kernel sharpness, not a single numerological
closed form.

The new useful object is the explicit level-7 Eisenstein sidecar

```text
E_7(tau) = (7 E_2(7 tau) - E_2(tau))/6.
```

Its q-expansion has divisor-fiber coefficients
`4*sigma_1(n)` away from multiples of 7 and
`4*sigma_1(n)-28*sigma_1(n/7)` on the 7-core.  That is exactly the kind of
payload the LRC proof needs: a coefficient source with the level-7 side channel
still attached.  The next serious move is to make this finite by producing
Toeplitz/LP certificate rows from the q-coefficients, then checking those rows
against HYP-3205/HYP-3224 packets and HYP-3227 trap weights.

The remote HYP-3215 update fits this exactly.  It names Fejer/Cohn-Elkies LP
construction as one of the main bounded gaps in the literature-facing route,
with a separate warning that the LRC induction base needs verification.  This
audit should be read as a construction attempt for the former, not as evidence
for the latter.  The Gamma0(7) coefficients are a possible finite LP basis;
they do not repair an unverified descent base.

The later mac-mini S75 update supplies the missing spatial basis.  It builds
the magic function as the comb-overlap Gram kernel

```text
K(p,q)=meas(D_p cap D_q)=<1_Dp,1_Dq>,
```

so positivity is automatic.  Its proved `K(1,q)=1/(7q)` row and single-arc
peeling recursion are the finite spatial counterpart of the Fejer Fourier
kernel.  This changes the next implementation target: do not build Gamma0(7)
coefficients in isolation.  Match them to the comb-overlap Gram rows and
leave the remaining order-3 triple-overlap constants as named debt.

The Stark/Dirichlet-L idea became more honest after computation.  The
discriminant `49` is real and important, but the direct `L(-1,chi)` values for
even primitive characters modulo 7 do not carry denominator `49`; they
contract to denominator `7`.  So the conductor-49 intuition should be stated
as "field discriminant and Bernoulli sampling grid" rather than as a proved
Dirichlet-L cap formula.  If a Stark formula exists here, it is likely a unit
or regulator statement whose rational collapse comes from the double-root
equioscillation, not the elementary `L(-1,chi)` table alone.

Beraha/Tutte and Mahler/Lehmer are real but secondary.  The Beraha constant
`B_7` satisfies

```text
B^3 - 5B^2 + 6B - 1 = 0,
```

and the Mahler measure of the de-Moivre cubic is exactly `B_7 - 1`.  This is a
good height gauge for perturbing the cubic, but it is not an LRC inequality
until it controls Fejer slack, Toeplitz margin, or Green trap conductance.

The subshift transfer-operator reading is similarly useful as a compression
model: `P(z)=1+...+z^6` has `P(z)P(z^-1)` autocorrelation coefficients
`7-|n|`, the same Fejer weights.  Consecutive/AP is the rank-one transfer
state.  A non-AP row should be compared by Perron defects or transfer-matrix
slack, not by the existence of a transfer analogy.

This also gives a better synthesis with HYP-3227.  HYP-3227 made a finite
trap-discharge graph whose Green-only version stays connected.  HYP-3229 says
the global dual above that graph should be Fejer/Gamma0(7), not raw
conductance.  The most concrete next test is:

```text
Does Fejer/Gamma0(7)/Gram-generated slack dominate every edge weight needed
to connect the non-AP traps to certificate coordinates?
```

If yes, the electrical graph becomes a finite shadow of the modular magic
certificate.  If no, the missing edge names the exact sidecar the proof still
lacks.

The tournament-analysis lesson is the same as the proof lesson.  Runners,
roots, q-coefficients, Beraha constants, Mahler measures, and L-values are all
possible vertex sets, but for this task they are shadows.  The useful vertices
are proof carriers, and the pairwise observable is how much certificate
payload survives.  Under that gauge the path is:

```text
Fejer/de-Moivre
-> comb-overlap Gram kernel
-> Johnson cap
-> Toeplitz/Green conductance bridge
-> Gamma0(7) Eisenstein sidecar
-> subshift transfer
-> Beraha/Tutte
-> Dirichlet-L/Stark
-> Mahler/Lehmer
-> raw numerology
```

That ordering is not an aesthetic preference.  It is a guardrail against
forgetting the LRC predicate.

Next work:

1. Emit a finite Gamma0(7)-coefficient / comb-overlap Gram certificate basis
   and test it on the HYP-3205/HYP-3224 bounded bank.
2. Compare Fejer/Gamma0(7)/Gram slack to HYP-3227 Green-only trap
   conductance weights and the four-trap precision-defect island.
3. Separate the S75 order-3 triple-overlap constants from the order-2 PSD
   kernel.
4. Connect the resulting finite basis to HYP-3215's Fejer/Cohn-Elkies
   LP/polyhedron-flatness Gap A while keeping the base-verification flag
   separate.
5. Add `beraha_height`, `mahler_height`, `dirichlet_l_sidecar`, and
   `subshift_perron_defect` as packet fields only after each has a named
   discharge target.

-> HYP-3229, HYP-3227, HYP-3215, HYP-3214, HYP-3213, HYP-3212, HYP-3205,
HYP-3203, HYP-3201, HYP-3162, HYP-3161, HYP-3160, THM-577, T1329,
LTI-329, LTT-229, OPEN-Q-108.
