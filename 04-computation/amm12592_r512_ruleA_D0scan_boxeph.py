"""R=512 rule A slack scan: D0 = 2, 3, 4 until closure. Exact floors."""
import json, time, io, contextlib, importlib.util, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)
from amm12592_allR_greedy_attractor_rule_boxeph import adm_clamp, Em
from amm12592_allR_family_toolbox_boxeph import *

R = 512
floor_prof = [iref.floor_gamma_star(R+i) for i in range(R)]

def attempt(D0):
    pr = [d + D0 for d in floor_prof]
    sig = qpow(R-1)
    blocks = []
    t0 = time.time()
    for i in range(R):
        d = pr[i]
        ideal = psub(sig, pshift(Em(R-2-i), 1))
        delta = adm_clamp(ideal[:d+1], d, 'tozero')
        D = bern_to_poly(delta, d)
        t = psub(sig, D)
        if t and t[0] != 0:
            print(f"D0={D0}: DIE row {i} const 2^{abs(t[0]).bit_length()} ({time.time()-t0:.0f}s)", flush=True)
            return None
        sig = t[1:] if t else []
        blocks.append(delta)
        if i % 128 == 0: print(f"  D0={D0} row {i} ({time.time()-t0:.0f}s)", flush=True)
    if sig != []:
        print(f"D0={D0}: residual nonzero", flush=True)
        return None
    adm = all(admissible(blocks[i], pr[i]) for i in range(R))
    idt = epoch_sum(R, blocks, pr) == qpow(R-1)
    print(f"D0={D0}: CLOSED adm={adm} identity={idt} ({time.time()-t0:.0f}s)", flush=True)
    if adm and idt:
        json.dump({"R": R, "profile": pr, "floor_profile": floor_prof, "D0": D0,
                   "blocks": blocks, "verified": True,
                   "source": f"greedy attractor rule A tozero, d_i=exact_floor+{D0}"},
                  open(f"amm12592_witness_R512_ruleA_D0_{D0}_boxeph.json", "w"))
        print(f"wrote amm12592_witness_R512_ruleA_D0_{D0}_boxeph.json", flush=True)
    return blocks

for D0 in (2, 3, 4):
    if attempt(D0) is not None:
        break
