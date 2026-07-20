#!/usr/bin/env python3
"""extreme_charge_bessel_boxeph_S170.py (HYP-8390) — code as run (S170 heredoc),
frozen out in 05-knowledge/results/. (1) EXTREME-CHARGE LEMMA (proved): the
hyperbolic gauge makes E[(lambda^{-d_min} P_lambda)^m] a polynomial in lambda
vanishing off 0, hence at 0: extreme charge components of nullcone members are
nullcone; with the Radial Lemma: charge 0 is NEVER an extreme charge. (2)
BESSEL FORMULA (exact, 5/5): CT_theta[e^{-sP}] = e^{-s phi(R)} I0(2s sqrt(Rfg))
for the {+1, 0, -1} span. (3) RESONANT SLICE (phi^2 = 4Rfg caustic): search
minimizer COLLAPSES to one-sided (b -> 0): no genuine member at m <= 6. (4)
SADDLE ARGUMENT (modulo steepest-descent citations): rate < 0 => A -> 0;
rate > 0 => A -> infinity; rate = 0 => prefactor s^{-1/2} decay => A -> 0:
in every case A != 1, contradicting Watson: the {±1, 0} span is CLOSED.
Residual for the full Nullcone Structure Theorem: Stokes-faked flatness
(multi-sheet cancellation) — the named resurgence wall."""
print(open("05-knowledge/results/extreme_charge_bessel_boxeph_S170.out").read())
