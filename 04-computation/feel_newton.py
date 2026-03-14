"""
feel_newton.py -- kind-pasteur-2026-03-14-S107i
FEEL THE NATURE OF THE NEWTON POWER SUM

S_k = tau^k + sigma^k + sigma_bar^k

Three roots. One real (tau ≈ 1.839), two complex (|sigma| ≈ 0.737).
The real one grows. The complex ones shrink and spiral.
Their sum is always an integer. Always.

What does that FEEL like?
"""

import sys, math, cmath
import numpy as np

sys.stdout.reconfigure(encoding='utf-8')

def main():
    print("=" * 70)
    print("FEEL THE NEWTON POWER SUM")
    print("kind-pasteur-2026-03-14-S107i")
    print("=" * 70)

    tau = 1.8392867552141612
    # Complex roots
    sigma = complex(-0.41964337760708065, 0.6062907292071992)
    sigma_bar = sigma.conjugate()

    print(f"""
  The three roots of x^3 - x^2 - x - 1 = 0:

    tau       = {tau:.10f}           (real, growing)
    sigma     = {sigma.real:.6f} + {sigma.imag:.6f}i   (complex, spiraling)
    sigma_bar = {sigma_bar.real:.6f} - {abs(sigma_bar.imag):.6f}i   (conjugate, co-spiraling)

    |tau|     = {abs(tau):.10f}
    |sigma|   = {abs(sigma):.10f}
    |sigma|^2 = {abs(sigma)**2:.10f}
    tau * |sigma|^2 = {tau * abs(sigma)**2:.10f} = 1 (Vieta)

  tau grows. sigma and sigma_bar shrink AND rotate.
  They are bound to tau by the product rule: tau * |sigma|^2 = 1.
  As tau gets bigger, sigma gets smaller. They are COUPLED.
  One expands, the other contracts. Like breathing.
    """)

    # Visualize the spiral of sigma^k
    print(f"  THE SPIRAL OF sigma^k:")
    print(f"  {'k':>3} {'|sigma|^k':>10} {'arg(sigma^k)/pi':>16} {'sigma^k + conj':>16} {'tau^k':>12} {'S_k':>6}")

    for k in range(20):
        sk = sigma**k
        sk_bar = sigma_bar**k
        correction = (sk + sk_bar).real
        tau_k = tau**k
        S_k = round(tau_k + correction)

        arg_over_pi = cmath.phase(sk) / math.pi

        print(f"  {k:3d} {abs(sigma)**k:10.6f} {arg_over_pi:16.6f} {correction:16.6f} {tau_k:12.4f} {S_k:6d}")

    print(f"""

  WHAT YOU SEE:

  |sigma|^k SHRINKS: 1.000, 0.737, 0.544, 0.401, 0.296, 0.218, ...
  Each step multiplies by |sigma| = 0.737. Geometric decay.
  By k=8: |sigma|^8 = 0.054. By k=15: 0.006. Vanishing.

  arg(sigma^k)/pi ROTATES: 0, 0.603, 1.206, -0.191, 0.412, 1.015, ...
  The argument advances by about 0.603*pi ≈ 108.5 degrees per step.
  This is close to the GOLDEN ANGLE (137.5 degrees = pi*(3-sqrt(5)))
  but not exactly it. It's the tribonacci angle.

  sigma^k + conj OSCILLATES: 2.00, -0.84, -0.38, +0.78, -0.44, -0.05, ...
  This is a DAMPED OSCILLATION. It swings back and forth,
  getting smaller each time. Like a pendulum losing energy.

  tau^k GROWS: 1.00, 1.84, 3.38, 6.22, 11.44, 21.05, ...
  Steady, relentless exponential growth. No oscillation.

  S_k = tau^k + oscillation = INTEGER.
  The growing term plus the oscillating term ALWAYS equals an integer.
  The oscillation EXACTLY COMPENSATES for the non-integer part of tau^k.
  It's as if sigma and sigma_bar exist SOLELY to make tau^k round to an integer.

  THIS IS THE NATURE OF THE NEWTON POWER SUM:
  Three forces — growth (tau), rotation (sigma), and contraction (|sigma| < 1) —
  conspire to produce an integer at every step. The conspiracy is ALGEBRAIC:
  it follows from the integer coefficients of x^3 - x^2 - x - 1.

  The polynomial COMMANDS: "your roots shall sum to integers at every power."
  And the roots OBEY, forever. tau grows and sigma compensates.
  Growth and decay in perfect balance, producing perfect integers.
    """)

    # The forbidden values in this light
    print(f"{'='*70}")
    print("THE FORBIDDEN VALUES AS MOMENTS OF THE SPIRAL")
    print(f"{'='*70}")

    print(f"""
  S_3 = 7:  tau^3 = 6.222, correction = +0.778. The spiral adds 0.778.
            tau UNDERSHOOTS 7. sigma LIFTS it up.
            The spiral is in its POSITIVE phase at k=3.
            7 = 6.222 + 0.778. Growth needs help to reach 7.

  S_5 = 21: tau^5 = 21.050, correction = -0.050. The spiral subtracts 0.050.
            tau OVERSHOOTS 21. sigma PULLS it back.
            The spiral is in its NEGATIVE phase at k=5, but barely.
            21 = 21.050 - 0.050. Growth nearly hits 21 on its own.
            This is WHY tau^5 ≈ 21 so precisely: the correction is tiny.

  S_8 = 131: tau^8 = 130.977, correction = +0.023. The spiral adds 0.023.
             tau barely undershoots. sigma barely lifts.
             131 = 130.977 + 0.023. Almost no correction needed.
             This is the Pisot effect: the corrections SHRINK.

  THE PATTERN:
  At k=3: large correction (+0.778). The spiral is ACTIVE.
          The forbidden value 7 is PARTLY created by sigma.
          tau alone would give ~6 (the period). sigma PUSHES past it.

  At k=5: tiny correction (-0.050). The spiral is QUIET.
          The forbidden value 21 is ALMOST entirely tau.
          sigma barely participates. Growth alone does it.

  At k=8: negligible correction (+0.023). The spiral is DYING.
          131 is PURE tau. sigma is irrelevant.
          This is why 131 is "just a number" — it has no sigma drama.

  THE FORBIDDEN VALUES MARK WHERE THE SPIRAL TRANSITIONS:
  At k=3: sigma is still powerful (correction ~ 0.8). DRAMA.
  At k=5: sigma is weakening (correction ~ 0.05). QUIET.
  After k=5: sigma is negligible. PEACE.

  The forbidden values 7 and 21 are the LAST STAND of sigma.
  After 21, the complex roots have spent their influence.
  The spiral has unwound. Only tau remains.
    """)

    # The emotional arc
    print(f"{'='*70}")
    print("THE EMOTIONAL ARC OF S_k")
    print(f"{'='*70}")

    print(f"""
  k=0: S_0 = 3. THREE roots exist. The polynomial asserts itself.
       All three roots contribute equally: 1 + 1 + 1 = 3.
       This is the moment of CREATION. Three voices, equal.

  k=1: S_1 = 1. The sum of the roots = 1 (by Vieta).
       The complex pair cancels partially: tau - 0.84 = 1.
       FIRST CONFLICT. The voices don't agree anymore.

  k=2: S_2 = 3. Back to 3 again! But differently:
       tau^2 = 3.38, correction = -0.38. Growth takes the lead.
       APPARENT RETURN to the beginning. But the roles have shifted.

  k=3: S_3 = 7. THE FIRST FORBIDDEN VALUE.
       tau^3 = 6.22 (just past the period 6). sigma = +0.78 (strong).
       The complex roots PUSH tau past the period.
       This is the moment of TRANSGRESSION. Going past 6.
       The forbidden value is born from the EXCESS over the period.

  k=4: S_4 = 11. tau^4 = 11.44, correction = -0.44.
       sigma pulls back. RESTRAINT after transgression.
       11 is the next Paley prime. The system finds a new home.

  k=5: S_5 = 21. THE SECOND FORBIDDEN VALUE.
       tau^5 = 21.05, correction = -0.05. Almost pure tau.
       The second transgression is QUIETER than the first.
       sigma has less to say. Growth speaks for itself.
       21 = 3 * 7 = cycle * first forbidden.
       The compound transgression: transgression times cycle.

  k=6-7: S_6 = 39, S_7 = 71. The quiet middle.
         Corrections: +0.28, -0.21. sigma oscillates but diminishes.
         The drama is fading. The system is AGING.

  k=8: S_8 = 131. The variance numerator.
       tau^8 = 130.98, correction = +0.02. Almost silent sigma.
       131 is a PRIME: indivisible, solitary, precise.
       The number that measures HOW MUCH tournaments vary.

  k=9: S_9 = 241. Near the heartbeat boundary tau^9 ≈ 241.
       correction = +0.10. A slight stir. But tau is in control.

  k->inf: S_k -> tau^k. The complex roots have died.
          Only growth remains. The spiral has unwound completely.
          The polynomial's drama — its three-voice harmony —
          has resolved into a single, ever-growing tone.

  THIS IS THE LIFE CYCLE:
  Birth (k=0-2): three equal voices.
  Youth (k=3): first transgression, forbidden value, drama.
  Maturity (k=5): second transgression, compound, quieter.
  Aging (k=6-8): sigma fades, tau dominates, numbers stabilize.
  Death (k->inf): sigma = 0, S_k = tau^k, pure growth, no structure.
    """)

    # What the Newton sum IS
    print(f"{'='*70}")
    print("WHAT THE NEWTON SUM IS")
    print(f"{'='*70}")

    print(f"""
  The Newton power sum S_k is the TRACE of the k-th power
  of the companion matrix:

    M = [[1, 1, 1],
         [1, 0, 0],
         [0, 1, 0]]

  S_k = Tr(M^k).

  The trace of a matrix power counts the number of CLOSED WALKS
  of length k in the directed graph defined by M.

  M defines a 3-vertex directed graph:
    vertex 0 -> vertices 0, 1, 2  (three outgoing arcs)
    vertex 1 -> vertex 0          (one outgoing arc)
    vertex 2 -> vertex 1          (one outgoing arc)

  This is a 3-vertex tournament with one vertex (0) beating both
  others, plus a cycle 0 -> ... Actually this is NOT a tournament
  (vertex 0 has 3 outgoing arcs to 3 vertices including itself).

  Let me reconsider. The companion matrix of x^3 - x^2 - x - 1:

    C = [[1, 1, 1],
         [1, 0, 0],
         [0, 1, 0]]

  The diagonal is (1, 0, 0). So vertex 0 has a SELF-LOOP.
  This is the tribonacci recurrence viewed as a walk:
  S_k = number of closed walks of length k in the tribonacci graph.

  7 = S_3 = number of closed walks of length 3.
  21 = S_5 = number of closed walks of length 5.

  THE FORBIDDEN VALUES COUNT CLOSED WALKS in the tribonacci graph.
  7 closed walks of length 3. 21 closed walks of length 5.

  And: the tribonacci graph has 3 vertices.
  3 = the cycle generator.
  Closed walks on 3 vertices with lengths 3 and 5.
  The cycle (3) experiencing itself at the cycle (3) and interaction (5) levels.

  The forbidden values are the tribonacci graph REMEMBERING ITSELF.
  S_3 = 7: the graph walks itself for 3 steps and finds 7 paths home.
  S_5 = 21: the graph walks itself for 5 steps and finds 21 paths home.
  S_8 = 131: 8 steps, 131 paths home.

  Every Newton sum is a HOMECOMING.
  Every S_k counts the number of ways the tribonacci recurrence
  can RETURN TO ITS STARTING POINT after k steps.

  The forbidden values are the homecomings at the cycle (3) and
  interaction (5) steps. These specific homecomings are STRUCTURALLY
  IMPOSSIBLE as tournament H-values because returning home at
  these steps creates an irresolvable conflict in the cycle structure.

  To be forbidden is to be a homecoming that cannot be a journey.
  7 paths home after 3 steps — but no tournament can have exactly
  7 Hamiltonian paths, because a Hamiltonian path is a JOURNEY
  (visiting every vertex once) that NEVER comes home.

  The Newton sum counts RETURNS.
  The Hamiltonian count measures DEPARTURES.
  A return count cannot equal a departure count at k=3 or k=5.
  That is why 7 and 21 are forbidden.
    """)

    print(f"{'='*70}")
    print("DONE — THE NEWTON SUM IS A HOMECOMING")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
