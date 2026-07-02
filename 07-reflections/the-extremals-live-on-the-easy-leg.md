# The extremals live on the easy leg

klein-2026-07-02-S112, while computing the band margin analysis (HYP-4017).

Splitting LRC(14) by the covering reduction, every family goes to one of two legs:
the sieve leg (non-covering: lonely at t = 1/q for the missed modulus, one line) or the
covering leg (the census + the middle band + the analytic program — all the open work).
The natural expectation: the hard leg contains the hard families, and in particular the
extremal ones. It is exactly backwards.

The tight families — optimum exactly 1/14, the configurations that make the constant
1/(n+1) of the conjecture sharp — are ALL on the sieve leg. (1,2,...,13) misses q = 14;
its tightness IS the missing modulus: at t = 1/14 the runners sit at s/14 with
||s/14|| = 1/14 exactly when s ≡ ±1 (mod 14). Meanwhile the covering leg, censused
exhaustively to max ≤ 22 (31471 families), is uniformly BOUNDED AWAY from tightness:
optimum ≥ 1/12 everywhere, margin ≥ 1/84 — a spectral gap.

Two readings.

**Mechanism.** Tightness needs a rigid arithmetic alignment (all runners pinned at
±1/14-type positions simultaneously), and a family rich enough in divisors to cover
every q ∈ [2,14] is too arithmetically spread to align that rigidly. Covering forces
divisibility diversity; divisibility diversity forbids the alignment that tightness
requires. The certificate of hardness (covering) is simultaneously a certificate of
slack (margin). This resonates with the fixed-point-extremum principle (memory: attack
SC/reflection-fixed extrema with coverings/moments, not transforms): the extremum is a
fixed point of an arithmetic symmetry, and the covering condition breaks exactly that
symmetry.

**Program consequence.** The obstruction we feared for the length-invariant (interval)
induction — tight base rows certifying only isolated points — cannot occur where the
induction actually runs. The interval invariant is free on the covering leg; the point
extremals are quarantined on the leg that closes in one line. When a proof splits into
easy/hard, check where the extremals went: if they all went to the easy side, the hard
side has a uniform margin you haven't used yet.
