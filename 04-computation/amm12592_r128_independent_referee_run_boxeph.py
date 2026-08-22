"""R = 128 gamma* floor witness through the INDEPENDENT hostile referee.

Loads amm12592_independent_witness_referee_boxeph.py (which re-implements
everything from the theorem statements: exact Z[x] arithmetic, Lucas-box
admissibility, and the EXACT floor d = floor(m gamma*) via integer sign
tests 5^d <= phi^{2m} < 5^{d+1} using phi^{2m} = (L_{2m} + F_{2m} sqrt5)/2 --
no floating point, no rational gamma* proxy) and runs the R = 128 checks:
  1. witness profile == true gamma* floor profile (Fib/Lucas integer test);
  2. verify_epoch: admissibility of all 128 blocks + exact epoch identity
     sum_i x^i Delta_i == (1-x)^127 in Z[x].
"""
import importlib.util, io, contextlib, json, os, sys

HERE = os.path.dirname(os.path.abspath(__file__))
spec = importlib.util.spec_from_file_location(
    "iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
mod = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):   # its own R<=64 self-run
    spec.loader.exec_module(mod)

w = json.load(open(os.path.join(HERE, "amm12592_floor_witness_R128_direct.json")))
R, prof, blocks = w["R"], w["profile"], w["blocks"]
true_prof = [mod.floor_gamma_star(R + i) for i in range(R)]
p_ok = (prof == true_prof)
print(f"R128 profile == exact gamma* floor profile (Fib/Lucas): {p_ok}")
res = mod.verify_epoch("R128 direct (beam1000 row119 + 8-row banded repair)",
                       R, prof, blocks)
print(f"profile_ok={res['profile_ok']} adm_ok={res['adm_ok']} "
      f"adm_bad={res['adm_bad']} identity_ok={res['identity_ok']} "
      f"unit_ok={res['unit_ok']} eff_rate={res['eff_rate']}")
print(f"signword: {res['signword']}")
ok = p_ok and res["profile_ok"] and res["adm_ok"] and res["identity_ok"]
print("R128 INDEPENDENT REFEREE: " + ("ALL CHECKS PASSED" if ok else "FAIL"))
sys.exit(0 if ok else 1)
