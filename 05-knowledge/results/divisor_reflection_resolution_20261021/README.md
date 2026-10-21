# Reflection obstruction and minimal source-preserving resolution

See [the proof](REFLECTION_RESOLUTION.md). The ten-divisor action cannot map equivariantly to graph complement plus path reversal; a minimal twelve-state resolution does. All-exponent orbit capacities are explicit.

Run `python3 verify_reflections.py > /tmp/reflections.json` and compare with `verification.json`; repeat with `python3 -O`. The checker reports 68,861 live checks. No full metagraph, divisibility-order, or Collatz equivalence is asserted.
