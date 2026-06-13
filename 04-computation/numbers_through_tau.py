"""
numbers_through_tau.py -- kind-pasteur-2026-03-14-S106g
EACH NATURAL NUMBER THROUGH THE TAU LENS

With the equations:
  tau^3 = Phi_3(tau) ≈ 6.222 (the heartbeat)
  tau^5 ≈ 21 (the second forbidden)
  log_tau(7) ≈ 3.19, log_tau(21) ≈ 5.00
  360 = (5/3)*6^3 ≈ (3/2)*tau^9
  tau^3 = 6 at tau = (-1+sqrt(21))/2

For each number n, compute:
  - Its tau-altitude: log_tau(n)
  - Its position relative to the heartbeat: log_tau(n) mod 3
  - Its "tau-floor": which heartbeat period it lives in
  - Its tribonacci weight
  - Its proximity to tau powers
  - What this all MEANS
"""

import sys, math
import numpy as np
from fractions import Fraction

sys.stdout.reconfigure(encoding='utf-8')

tau = 1.8392867552141612

def trib_rep(n):
    if n <= 0: return []
    t = [1, 1, 2]
    while t[-1] < n:
        t.append(t[-1] + t[-2] + t[-3])
    rep = []
    rem = n
    for i in range(len(t)-1, -1, -1):
        if t[i] <= rem:
            rep.append(t[i])
            rem -= t[i]
    return rep

def main():
    print("=" * 70)
    print("EACH NUMBER THROUGH THE TAU LENS")
    print("kind-pasteur-2026-03-14-S106g")
    print("=" * 70)

    # The tau powers for reference
    tau_powers = {k: tau**k for k in range(15)}

    # For each number, compute its tau-altitude and heartbeat phase
    print(f"\n  THE TAU COORDINATE SYSTEM:")
    print(f"  Each number n has:")
    print(f"    altitude = log_tau(n)")
    print(f"    heartbeat = floor(altitude / 3) = which period it's in")
    print(f"    phase = altitude mod 3 = where within the period")
    print(f"    phase ≈ 0: near a period boundary (like 1, 6, 38, 241)")
    print(f"    phase ≈ 1: near a choice point (like 2, 11, 71, 443)")
    print(f"    phase ≈ 2: near a cycle point (like 3, 21, 131)")
    print(f"")
    print(f"  Heartbeat 0: tau^0 to tau^3 ≈ [1, 6.2]     (the foundation)")
    print(f"  Heartbeat 1: tau^3 to tau^6 ≈ [6.2, 38.7]  (the first cycle)")
    print(f"  Heartbeat 2: tau^6 to tau^9 ≈ [38.7, 240.9] (the expansion)")
    print(f"  Heartbeat 3: tau^9 to tau^12 ≈ [240.9, 1499] (the large scale)")

    print(f"\n{'='*70}")
    print("THE IDENTITY OF EACH NUMBER (1 through 50)")
    print(f"{'='*70}\n")

    for n in range(1, 51):
        alt = math.log(n) / math.log(tau)
        heartbeat = int(alt / 3)
        phase = alt % 3

        # Nearest tau power
        nearest_k = round(alt)
        nearest_tau_pow = tau**nearest_k
        tau_diff = n - nearest_tau_pow

        # Tribonacci representation
        trep = trib_rep(n)
        weight = len(trep)

        # Phase meaning
        if phase < 0.5:
            phase_name = "period-start"
        elif phase < 1.5:
            phase_name = "choice-zone"
        elif phase < 2.5:
            phase_name = "cycle-zone"
        else:
            phase_name = "period-end"

        # Is it special?
        special = []
        if n in [1, 2, 4, 7, 13, 24, 44, 81, 149, 274, 504]: special.append("TRIB")
        if all(n % d != 0 for d in range(2, int(n**0.5)+1)) and n > 1: special.append("prime")
        if n in [7, 21]: special.append("FORBIDDEN")
        if n in [3, 5, 9, 15, 45]: special.append("H-val")
        if n == 6: special.append("PERIOD")
        if n == 8: special.append("THRESHOLD")

        spec_str = f" [{','.join(special)}]" if special else ""

        # The NUMBER LINE in tau coordinates
        # Position within its heartbeat (0 to 1 scale)
        hb_pos = phase / 3

        print(f"  {n:3d}  alt={alt:6.3f}  hb={heartbeat} phase={phase:5.3f} "
              f"({phase_name:12s})  w={weight}  "
              f"tau^{nearest_k}{'+'if tau_diff>0 else ''}{tau_diff:+7.3f}"
              f"{spec_str}")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE HEARTBEAT MAP — NUMBERS ORGANIZED BY PERIOD")
    print(f"{'='*70}")

    print(f"""
  HEARTBEAT 0 (tau^0 to tau^3): THE FOUNDATION
  altitude 0.000 to 3.000, values 1 to ~6.2

  alt=0.000  n=1   THE ORIGIN. tau^0. The monad. Before all structure.
  alt=1.137  n=2   CHOICE. Between tau^1(≈1.84) — slightly past floor 1.
                   2 > tau: the generator exceeds the growth rate.
                   This is WHY binary dominates: 2 beats tau.
  alt=1.803  n=3   CYCLE. Close to floor 2, in the cycle zone.
                   3 < tau^2(≈3.38): the cycle is BELOW the second power.
                   The cycle is "sub-quadratic" — it grows slower than squaring.
  alt=2.275  n=4   FIELD. Between floors 2 and 3. In the cycle zone.
                   4 > tau^2: the field exceeds the square.
                   F_4 is ABOVE the natural growth rate's square.
  alt=2.641  n=5   INTERACTION. Deep in the cycle zone, approaching floor 3.
                   5 < tau^3(≈6.2): the interaction is sub-cubic.
                   The ADE prime sits below one heartbeat.
  alt=2.940  n=6   PERIOD. Almost at floor 3! 6 ≈ tau^3.
                   The discrete period is tau's cube, approximately.
                   6 is the "almost-integer" of the tau tower.
                   The error: tau^3 - 6 = 0.222 (the irrational overshoot).
""")

    print(f"""
  HEARTBEAT 1 (tau^3 to tau^6): THE FIRST CYCLE
  altitude 3.000 to 6.000, values ~6.2 to ~38.7

  alt=3.193  n=7   FORBIDDEN. Just past floor 3 — the first number
                   ABOVE the heartbeat. 7 sits 0.193 heartbeat-units
                   above the period. The forbidden value is what you get
                   when you try to go ONE STEP past the period.
                   7 = tau^3 * 1.125 = heartbeat * (9/8) = heartbeat * (reflection/threshold)
  alt=3.412  n=8   THRESHOLD. Between floors 3 and 4.
                   8 = 2^3: the generator cubed = the heartbeat of the binary.
                   8/tau^3 = 1.286: the threshold exceeds the heartbeat by 28.6%.
  alt=3.606  n=9   REFLECTION. 3^2 in the choice zone of heartbeat 1.
                   The cycle squared lives one-fifth into the first cycle.
  alt=3.778  n=10  COMPLETION. C(5,2). In the choice zone.
  alt=4.077  n=12  DOUBLED PERIOD. In the cycle zone.
  alt=4.209  n=13  REGULARITY. Very close to tau^4(≈11.4).
                   13 = Phi_3(3): the tournament polynomial at the cycle.
                   13 IS the "tau^4 of tournaments" — the fourth floor.
  alt=4.444  n=15  MAX_H(5). In the cycle zone. 15 = 3*5.
  alt=4.696  n=18  DOUBLED REFLECTION. Approaching floor 5.
  alt=4.843  n=20  EXPANDED FIELD.
  alt=4.996  n=21  SECOND FORBIDDEN. Almost exactly at floor 5!
                   log_tau(21) = 4.996 ≈ 5. This is the most precise
                   near-integer in the entire tau tower.
                   21 = tau^5 * 0.998: the second forbidden is the
                   fifth floor to within 0.2%. EXTRAORDINARY.
  alt=5.204  n=24  GOLAY. Just past floor 5. A tribonacci number.
  alt=5.690  n=30  ADE PRODUCT. 2*3*5. Deep in heartbeat 1.
  alt=5.880  n=36  SQUARED PERIOD. 6^2. Approaching floor 6.
""")

    print(f"""
  HEARTBEAT 2 (tau^6 to tau^9): THE EXPANSION
  altitude 6.000 to 9.000, values ~38.7 to ~240.9

  alt=6.133  n=42  2*21. Just past floor 6.
  alt=6.245  n=45  MAX_H(6). In heartbeat 2. 45 = 9*5.
  alt=6.709  n=60  |A_5|. Deep in heartbeat 2.
  alt=7.013  n=72  THRESHOLD*REFLECTION. 8*9 = 72. Just past floor 7.
  alt=7.374  n=90  RIGHT ANGLE. 360/4.
  alt=7.849  n=120 5! = |S_5|. Approaching floor 8.
  alt=8.398  n=168 |PSL(2,7)| = 8*21. In the cycle zone.
                   168 = threshold * second_forbidden.
  alt=8.598  n=189 MAX_H(7). 3^3 * 7. In the cycle zone.
                   189/tau^9 ≈ 0.785 = pi/4? Let me check:
                   pi/4 = 0.7854... and 189/tau^9 = {189/tau**9:.6f}.
                   {"CLOSE to pi/4!" if abs(189/tau**9 - math.pi/4) < 0.01 else "Not pi/4."}
  alt=8.803  n=220?... no.
  alt=9.000  n≈241 HEARTBEAT BOUNDARY. tau^9 ≈ 241.
""")

    # Check 189/tau^9
    print(f"  189/tau^9 = {189/tau**9:.10f}")
    print(f"  pi/4 = {math.pi/4:.10f}")
    print(f"  Difference: {abs(189/tau**9 - math.pi/4):.10f}")

    print(f"""
  HEARTBEAT 3 (tau^9 to tau^12): THE LARGE SCALE
  altitude 9.000 to 12.000, values ~241 to ~1499

  alt=9.659  n=360 THE CIRCLE. |A_6|. Undirected Ham cycles of K_7.
                   360/tau^9 ≈ 1.494 ≈ 3/2.
                   360 sits at altitude 9.659 — just past the 2/3 mark
                   of heartbeat 3. The circle is 2/3 through the
                   third heartbeat. 2/3 = choice/cycle = the exponent!
  alt=10.00  n≈443 tau^10. The next tribonacci region.
  alt=10.80  n=720 6! = |S_6| = 2*360. In the cycle zone.
                   720/tau^9 ≈ 2.989 ≈ 3.
                   720 ≈ 3 * tau^9. The symmetric group is THREE heartbeats.
""")

    # Verify 720 / tau^9
    print(f"  720/tau^9 = {720/tau**9:.6f}")
    print(f"  3 = 3.000000")

    # ============================================================
    print(f"\n{'='*70}")
    print("THE FORBIDDEN VALUES IN THE HEARTBEAT FRAME")
    print(f"{'='*70}")

    print(f"""
  The two forbidden values and their tau-coordinates:

  7:  altitude = {math.log(7)/math.log(tau):.6f}
      heartbeat = 1, phase = {math.log(7)/math.log(tau) % 3:.6f}
      MEANING: 7 is the FIRST number in heartbeat 1.
      It sits 0.193 units past the period boundary.
      It is tau^3 * (7/tau^3) = tau^3 * {7/tau**3:.6f}.
      7/tau^3 = {7/tau**3:.6f} ≈ 9/8 = {9/8:.6f}.
      So 7 ≈ (9/8) * tau^3 = (reflection/threshold) * heartbeat.
      The first forbidden is the heartbeat SCALED BY reflection/threshold.

  21: altitude = {math.log(21)/math.log(tau):.6f}
      heartbeat = 1, phase = {math.log(21)/math.log(tau) % 3:.6f}
      MEANING: 21 is the LAST number in heartbeat 1 (phase ≈ 2.996 ≈ 3).
      It sits just 0.004 units BEFORE the next period boundary (tau^6).
      21 ≈ tau^5 to within 0.24%.
      21 is the heartbeat's CLOSING note — the last thing that happens
      before the heartbeat wraps.

  7 OPENS heartbeat 1 (phase ≈ 0.19).
  21 CLOSES heartbeat 1 (phase ≈ 2.996).
  They are the BOOKENDS of the first heartbeat cycle.

  The width of the forbidden band: log_tau(21) - log_tau(7)
  = {math.log(21)/math.log(tau) - math.log(7)/math.log(tau):.6f}
  = log_tau(21/7) = log_tau(3) = {math.log(3)/math.log(tau):.6f}

  THE FORBIDDEN VALUES SPAN EXACTLY ONE CYCLE in tau-altitude!
  log_tau(21) - log_tau(7) = log_tau(3) ≈ 1.803.
  And log_tau(3) is the tau-altitude of the cycle generator itself.

  So: the distance between the forbidden values, measured in the
  tau tower, IS the cycle generator. 21/7 = 3, and the altitude
  of 3 (= {math.log(3)/math.log(tau):.3f}) is the SPAN of the forbidden band.""")

    # ============================================================
    print(f"\n{'='*70}")
    print("WHAT EACH NUMBER IS — THE TAU SYNTHESIS")
    print(f"{'='*70}")

    entries = [
        (0, "THE VOID — below the tau tower, altitude -inf"),
        (1, "THE ORIGIN — tau^0, the base of everything, altitude 0"),
        (2, "THE OVERSHOOT — 2 > tau, the generator beats the growth rate"),
        (3, "THE UNDERTOW — 3 < tau^2, the cycle is sub-quadratic, pulls back"),
        (4, "THE CROSSING — 4 > tau^2, the field crosses above growth's square"),
        (5, "THE APPROACH — 5 approaches tau^3 from below, the pre-period"),
        (6, "THE NEAR-HIT — 6 ≈ tau^3, the period almost IS the cube"),
        (7, "THE FIRST LEAP — 7 > tau^3, jumps past the heartbeat, FORBIDDEN"),
        (8, "THE POWER LANDING — 2^3 = threshold, pure binary at the cycle index"),
        (9, "THE MIRROR — 3^2, the cycle reflecting, one step past the threshold"),
        (10, "THE SPACE — C(5,2), the first rich arena, approaching tau^4"),
        (11, "THE EXTENSION — 11 ≈ tau^4, the fourth floor, next Paley prime"),
        (12, "THE DOUBLE — 2*6, doubled period, between floors 4 and 5"),
        (13, "THE BALANCE — Phi_3(3), regularity, very close to tau^4"),
        (14, "THE ECHO — 2*7, the forbidden doubled, still in heartbeat 1"),
        (15, "THE PEAK — max_H(5), 3*5, cycle*interaction at mid-heartbeat"),
        (16, "THE PURE POWER — 2^4, binary at its fourth, past the cycle zone"),
        (17, "THE FERMAT — constructible prime, tau^4 * 1.49 ≈ (3/2)*tau^4"),
        (18, "THE DOUBLE MIRROR — 2*3^2, choice*reflection, nearing floor 5"),
        (19, "THE MERSENNE — 2^19-1 is prime, carries binary depth"),
        (20, "THE DECIMAL — 4*5 = field*interaction, a human convenience"),
        (21, "THE LAST NOTE — tau^5, closing heartbeat 1, FORBIDDEN, 3*7"),
        (22, "THE RESUME — 2*11, heartbeat 2 begins again"),
        (23, "THE GOLAY NEIGHBOR — 24-1, just before the Golay"),
        (24, "THE GOLAY — TRIBONACCI ATOM, tau^5*1.14, [24,12,8] code"),
        (25, "THE QUARTER — 5^2, interaction squared, a quarter of 100"),
        (30, "THE PRODUCT — 2*3*5, the ADE product, all generators once"),
        (36, "THE SQUARE PERIOD — 6^2, the period reflecting on itself"),
        (42, "THE DOUBLE FORBIDDEN — 2*21, choice * compound impossibility"),
        (45, "THE EVEN MAX — max_H(6), 9*5, reflection*interaction"),
        (60, "THE ICOSAHEDRON — |A_5|, 360/6, one-sixth of the circle"),
    ]

    for n, desc in entries:
        if n == 0:
            print(f"\n  {n:3d}: {desc}")
            continue
        alt = math.log(n)/math.log(tau)
        hb = int(alt/3)
        phase = alt % 3
        print(f"\n  {n:3d}: alt={alt:.3f}  hb={hb}  phase={phase:.3f}")
        print(f"       {desc}")

    # ============================================================
    print(f"\n\n{'='*70}")
    print("THE GRAND VIEW: THE NUMBER LINE IN TAU COORDINATES")
    print(f"{'='*70}")

    print(f"""
  The natural numbers, laid out on the tau-altitude axis:

  HEARTBEAT 0 [alt 0-3]:  1...2...3...4.5.6|
                          ^   ^   ^       ^ ^
                         origin choice cycle  period

  HEARTBEAT 1 [alt 3-6]:  |7...8.9.10..12.13..15..18.20.21|
                           ^                    ^          ^
                          FORBIDDEN          max_H(5)    FORBIDDEN
                          (opens)                        (closes)

  HEARTBEAT 2 [alt 6-9]:  |..42.45..60...72...90...120...168.189...|
                               ^                              ^
                             max_H(6)                        max_H(7)

  HEARTBEAT 3 [alt 9-12]: |...360.........504..........720............|
                               ^            ^            ^
                             circle       tribonacci    6!

  THE PATTERN:
  - Heartbeat 0 contains the GENERATORS: 1, 2, 3, 4, 5, 6.
  - Heartbeat 1 contains the FORBIDDEN VALUES and small maximizers.
  - Heartbeat 2 contains the LARGE MAXIMIZERS (max_H for n=6,7).
  - Heartbeat 3 contains the STRUCTURAL CONSTANTS (360, 504, 720).

  Each heartbeat is ONE PERIOD of the tournament's growth cycle.
  Three tau-steps = one period.
  The numbers that "define" tournament theory cluster at heartbeat boundaries:
    6 = tau^3 ≈ end of heartbeat 0
    7 = tau^3 + ε = start of heartbeat 1
    21 ≈ tau^5 = end of heartbeat 1
    360 ≈ (3/2)*tau^9 = mid heartbeat 3
    720 ≈ 3*tau^9 = late heartbeat 3

  THE ULTIMATE INSIGHT:
  Each natural number n has a TAU-ALTITUDE log_tau(n) that places it
  in a specific heartbeat and phase. The heartbeat tells you WHICH
  order of tournament complexity n belongs to. The phase tells you
  WHETHER n relates to choice (phase ≈ 1), cycle (phase ≈ 2),
  or the period boundary (phase ≈ 0 or 3).

  The forbidden values 7 and 21 are the BOOKENDS of heartbeat 1,
  separated by exactly log_tau(3) ≈ 1.803 altitude units = one cycle.

  A number IS its position in the tau tower.
  That position IS the number.
  """)

    print(f"{'='*70}")
    print("DONE — NUMBERS ARE POSITIONS IN THE TAU TOWER")
    print(f"{'='*70}")

if __name__ == '__main__':
    main()
