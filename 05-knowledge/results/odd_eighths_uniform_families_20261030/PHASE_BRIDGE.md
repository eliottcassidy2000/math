# The same phase bit controls two actual responses

**PROVED elementary synthesis from E3 and the incoming Session115 harmonic
completion theorem.** Same-agent review; no new general representation
theorem or dynamical conjugacy is claimed.

Incoming [HARMONIC_COMPLETION.md](../session115_odd_eighth_classifications_20261030/lrc/HARMONIC_COMPLETION.md)
at893f84da21 proves the exact family

    V_m={1,2,3,4,5,6,8} union {7m,14m,...,49m}, m>=2,
    M(V_m)=1/8.

Its full maximizing-time census and physical mass formula are inherited.
At a maximizing point t, define the ACTUAL tail phase

    c(t)=8 frac(7m t), in {1,3,5,7}.

The inverse index of c modulo8 equals c. The contribution of that physical
point to the mass derivative is w_8(c)/(7m), where

    w_8(c)=1/c+1/(8-c)=88/105+(32/105) chi_8(c).           (1)

The incoming proof retains the original displacement index and directional
tail owners; the phase alone does not locate t or give its multiplicity.

## B1. An explicit affine response map

For our3n+1 source-character functions, define the scaled conditional energy

    K_M(c)=3^M E_[R=c mod8] |S_M+1/2|^2, M>=3.

E3 proves K_M(c)=1/4+(-1)^M(3sqrt(2)/64)chi_8(c). Combining with(1),

    w_8(c)=88/105+(1024sqrt(2)/315)(-1)^M[K_M(c)-1/4].    (2)

This is an exact affine decoder between two observables on the SAME named
phase set. It retains the clock parity M and does not transport trajectories,
positive integer seeds, or LRC endpoint locations.

Let N_m be the inherited number of actual maximizing points, and average
the conditional energies over those points with their actual multiplicity:

    Kbar_M=(1/N_m) sum_(t maximizes V_m) K_M(c(t)).

Then the entire physical terminal mass slope is

    B(V_m)=N_m/(7m) *
       [88/105+(1024sqrt(2)/315)(-1)^M(Kbar_M-1/4)].       (3)

The averaging in(3) is over ACTUAL maximizing points, not a uniform average
over four abstract phases. Replacing it by the latter deletes the incoming
signed discrepancy term. The checker recovers(3) by independently computing
literal safe intervals and the conditional Gram formula for m=2,3,6,10,
at both M parities. This checks the target physical derivative.

## B2. Complete reflection-even observable classification

Every scalar function on {1,3,5,7} invariant under c->-c has uniquely the
form alpha+beta chi_8(c): reflection has exactly the two orbits{1,7},{3,5}.
If beta is nonzero, its value determines precisely that orbit. This elementary
two-dimensional algebra explains why(1) and E3 admit an affine decoder.
They measure the same reflection-even bit, with different scales and clocks.

The full oriented phase still needs the other bit. E2's outermost complex
correlation distinguishes all four residues, whereas neither w_8 nor K_M
does. Thus the weakest additional selector after either scalar response is
a sign choosing between c and-c. An arbitrary such selector suffices on
this four-element carrier; it does not reconstruct the physical LRC point,
which also needs its owner/displacement address.

Forgetting M parity in(2) has an explicit hostile: K_M(1)=K_(M+1)(3)
after scaling as defined, but w_8(1)=8/7 and w_8(3)=8/15 differ. The clock
is therefore an essential sidecar to this cross-response decoder.

## Boundaries and the other incoming eighths

Incoming115 also classifies coherent Fano-gauge odd-eighth Walsh readings
and even-order masses of F_2^3 semidirect H for odd H<=GL_3(2). Those
fractions are probabilities of different consumers. No phase-labeled map
from those objects to c(t) or to the conditional source cylinder has been
provided, so equations(2)--(3) do not apply to them. They remain precisely
scoped neighboring classifications, not extra instances claimed by analogy.

The proved bridge here is stronger than a matching denominator because it
specifies source phase, target observable, clock, multiplicity, exact decoder
and an independently measured physical consequence. It is narrower than an
equivalence of LRC and Collatz. Both major problems remain open.
