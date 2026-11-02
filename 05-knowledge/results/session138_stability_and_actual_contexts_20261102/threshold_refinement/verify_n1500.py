#!/usr/bin/env python3
"""Exact endpoints for PROOF_N1500.md; no tournament enumeration."""
from fractions import Fraction as F
from math import factorial
import hashlib
from pathlib import Path

checks = {}
checks['exp4_partial_above_50'] = sum(F(4)**j/factorial(j) for j in range(8)) > 50
checks['kappa_at_1500_below_1_over_187'] = F(8,1499) < F(1,187)
checks['kappa_below_eta'] = F(1,187) < F(1,44)
checks['variance_upper_n_over_7'] = F(1,8)+F(4,1499) < F(1,7)
heavy = F(49*1499,1600)-F(16*1500,1499)-F(1500,88)-F(1500,120)
checks['heavy_mass_endpoint'] = heavy == F(9253539,26382400) and heavy > 0
checks['heavy_mass_positive_slope'] = F(49,1600)-F(1,88)-F(1,120) == F(577,52800) > 0
checks['epsilon_endpoint'] = 1500000*3**50 < 4**50
checks['epsilon_log_slope_bound'] = F(1,1500) < F(1,120)
checks['light_absorption_slack'] = F(49,200)-F(1,11)-F(1,1750)-F(3,20) == F(271,77000) > 0
checks['contradiction_linear_term'] = F(44,187) == F(4,17)
checks['contradiction_squared_slack'] = F(13,17)**2-F(320,561) == F(137,9537) > 0
for name, passed in checks.items():
    print(f'{name}: {"PASS" if passed else "FAIL"}')
failed = [name for name, passed in checks.items() if not passed]
if failed:
    raise RuntimeError("Exact endpoint verification failed: " + ", ".join(failed))
proof = Path(__file__).with_name('PROOF_N1500.md')
print(f'proof_sha256: {hashlib.sha256(proof.read_bytes()).hexdigest()}')
print(f'exact_endpoint_checks: {len(checks)} / {len(checks)} PASS')
