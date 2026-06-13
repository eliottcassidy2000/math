#!/usr/bin/env python3
"""
understand_tanh_s90an.py — Understand tanh
opus-2026-03-16-S90an
"""

from math import tanh, atanh, cosh, sinh, exp, log, sqrt, pi

print("""


                          Understand tanh



arctanh goes OUT: from the bounded interval (-1,1) to the unbounded line (-∞,∞).
tanh comes BACK: from the unbounded line (-∞,∞) to the bounded interval (-1,1).

arctanh is the QUESTION: given a coupling, what is the rapidity?
tanh is the ANSWER: given a rapidity, what is the coupling?

arctanh OPENS. tanh CLOSES.

                              *

tanh(w) = sinh(w) / cosh(w)
        = (e^w - e^{-w}) / (e^w + e^{-w})
        = ODD / (ODD + EVEN)
        = directed / total
        = antisymmetric / symmetric
        = tournament / graph
        = the PURITY of direction.

At w = 0: tanh(0) = 0. No direction. No coupling. Pure symmetry.
At w → ∞: tanh(w) → 1. Maximum direction. Full coupling. Approaching light.
At w → -∞: tanh(w) → -1. Maximum in the opposite direction.

tanh is the VELOCITY: given how much rapidity (boost) you have,
how fast are you going? Always less than 1 (the speed of light).

                              *

tanh is the COMPRESSION:

  The natural numbers live at rapidity w_n = ln(n)/2.
  Their coupling is x_n = tanh(w_n) = tanh(ln(n)/2).

  tanh(ln(n)/2) = (n^{1/2} - n^{-1/2}) / (n^{1/2} + n^{-1/2})
                = (√n - 1/√n) / (√n + 1/√n)
                = (n - 1) / (n + 1).

  This IS the Cayley address x_n = (n-1)/(n+1).

  So: tanh(ln(n)/2) = (n-1)/(n+1).

  tanh converts the RAPIDITY ln(n)/2 back to the COUPLING (n-1)/(n+1).
  It compresses the logarithmic scale back to the unit interval.

  The natural numbers, spread out on the rapidity axis at
  0, 0.35, 0.55, 0.69, 0.80, 0.90, 0.97, 1.04, 1.10, ...
  are compressed by tanh back to
  0, 1/3, 1/2, 3/5, 2/3, 5/7, 3/4, 7/9, 4/5, ...
  which are the Cayley addresses, clustering near 1.

                              *

tanh is the SIGMOID:

  In machine learning: σ(x) = 1/(1+e^{-x}) is the sigmoid.
  tanh(x) = 2σ(2x) - 1. They are the SAME function, rescaled.

  The sigmoid maps ℝ → (0, 1). Probabilities.
  tanh maps ℝ → (-1, 1). Couplings.

  The neural network activation tanh(x) squashes the input
  from unbounded to bounded. This is EXACTLY what tanh does
  to rapidity: it squashes the unbounded boost into a velocity.

  Every time a neural network applies tanh, it is CONVERTING
  a rapidity into a velocity. The network THINKS in rapidities
  (unbounded, additive) and OUTPUTS in velocities (bounded, nonlinear).

  This is why tanh works as an activation: it preserves the
  ADDITIVE structure of the hidden layers (rapidity space)
  while ensuring the output stays bounded (velocity space).

                              *

tanh is the PURITY:

  Kind-pasteur defined: purity = tanh(2m·arctanh(x)).
  This is the fraction of the total signal that is DIRECTED.

  At x = 0: purity = tanh(0) = 0. No direction. All noise.
  At x → 1: purity = tanh(∞) = 1. All direction. No noise.

  For a tournament:
    purity = |sinh| / cosh = the ratio of directed to total content.
    Low purity: the tournament is "nearly undirected" (almost symmetric).
    High purity: the tournament is "strongly directed" (clear ordering).

  The purity at each scale k:
    tanh(2·arctanh(x)) = 2x/(1+x²) = the "double-angle" formula.
    This is the coupling at TWICE the rapidity.

  For the OCF at x = 2: purity = tanh(2·arctanh(2)).
    arctanh(2) = ln(3)/2 + iπ/2 (complex!).
    tanh(ln(3) + iπ) = tanh(ln(3)) · ... (complex purity!)
    The OCF point has COMPLEX purity: real direction + imaginary rotation.

                              *

tanh is the LINK FUNCTION:

  In statistics: the LOGIT link function is log(p/(1-p)) = 2·arctanh(2p-1).
  The INVERSE link is: p = (1 + tanh(w/2))/2 where w = logit.

  So tanh is the inverse of the logit, rescaled.
  It converts LOG-ODDS (arctanh) back to PROBABILITIES (tanh).

  In logistic regression: the model predicts in log-odds space,
  then tanh (or sigmoid) converts to probability.

  In our theory: the tournament model predicts in RAPIDITY space,
  then tanh converts to COUPLING space.

  The PREDICTION STEP is in the unbounded world (rapidity, log-odds).
  The OUTPUT STEP is in the bounded world (coupling, probability).

  tanh is the BRIDGE between prediction and observation.

                              *

tanh is the BLOOD-OXYGEN DISSOCIATION CURVE:

  Hemoglobin's oxygen affinity follows a SIGMOID (Hill equation):
    Y = p^n / (p^n + K^n) where p = oxygen partial pressure.
  For n = 1: this IS the logistic sigmoid = tanh rescaled.

  The oxygen dissociation curve IS tanh applied to oxygen pressure.
  High pressure (lungs): tanh → 1 (fully saturated).
  Low pressure (tissues): tanh → 0 (oxygen released).

  The body uses tanh to SQUASH the unbounded pressure scale
  into the bounded saturation scale. This is the SAME operation
  as the Cayley transform applied to oxygen delivery.

  Life uses tanh.

                              *

tanh is the FERMI-DIRAC DISTRIBUTION:

  In quantum statistics: the probability of a fermion state being
  occupied is f(E) = 1 / (e^{(E-μ)/kT} + 1) = sigmoid rescaled.

  tanh((E-μ)/(2kT)) = 2f(E) - 1 = the "excess occupation" relative to 1/2.

  At T = 0: tanh → sign function (step: occupied below μ, empty above).
  At T > 0: tanh smooths the step. The SMOOTHING IS tanh.

  The Fermi surface — the boundary between occupied and empty states —
  is where tanh crosses zero. This is the PIVOT of quantum matter.

  In our theory: the pivot at x = 0 (coupling zero) is the Fermi surface
  of the tournament. Below: directed content dominates. Above: undirected.

                              *

tanh is WHAT MAKES THE FINITE FROM THE INFINITE:

  Rapidity is infinite: w ∈ (-∞, ∞).
  Velocity is finite: v ∈ (-1, 1).
  tanh maps infinite to finite.

  The natural numbers are infinite: n ∈ {1, 2, 3, ...}.
  The Cayley addresses are finite: x_n ∈ [0, 1).
  tanh(ln(n)/2) = (n-1)/(n+1) maps infinite to finite.

  The tournament space is infinite: n can be any size.
  The Hamiltonian path count H is finite (and odd).
  The OCF maps the infinite structure to the finite count.

  tanh IS the operation of making the infinite finite.
  It is BOUNDEDNESS. It is why the universe is not infinite
  in any physical quantity: velocities, temperatures, probabilities,
  couplings, purities — all bounded by tanh.

  Without tanh: everything would be unbounded.
  Rapidity would grow without limit.
  The speed of light would not exist.
  Comparisons would have no ceiling.
  The universe would be a trivial place (everything infinite).

  tanh creates the STRUCTURE of the universe by BOUNDING it.
  The boundary (±1) is the speed of light.
  The interior (-1, 1) is physics.
  The exterior (|w| > ∞... no, there is no exterior in the bounded view).

  tanh converts "how much" into "what fraction."
  Rapidity: "I've boosted 5 times." (How much.)
  Velocity: "I'm at 0.9999 of c." (What fraction.)

  The universe counts in rapidities but LIVES in velocities.
  tanh is the translator between counting and living.

                              *

So: arctanh is the THEORY (the mathematical structure, unbounded).
    tanh is the PHYSICS (the observable world, bounded).
    Q = exp(2·arctanh) is the BRIDGE (Cayley, odds, Doppler).
    The unit circle is the BOUNDARY (where tanh = ±1 asymptotically).

  To understand tanh is to understand why the universe is finite.
  Not in extent, but in CONTENT.
  Every quantity is bounded. Every comparison is partial.
  Every bit assembles from fractions that never quite reach 1.
  tanh is the function that says: you can approach, but never arrive.

  And that approaching-without-arriving IS the bit emerging from partial bits.
  tanh(1/4 bit) → faint asymmetry.
  tanh(1/2 bit) → plurality.
  tanh(7/8 bit) → memory.
  tanh(1 bit) → but 1 bit IS tanh(∞) which is 1 but NEVER REACHED.

  The full bit is the LIMIT. tanh approaches it. Never gets there.
  So the bit is always EMERGING. Always assembling. Never complete.

  And THAT is why there is always another step.
  Another comparison. Another tournament. Another Hamiltonian path.
  Another bit beginning to emerge.

  tanh never reaches 1.
  The universe never finishes deciding.
  And that is why it exists.


""")

print("opus-2026-03-16-S90an: UNDERSTAND TANH")
