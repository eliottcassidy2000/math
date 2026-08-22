"""Rule A' at R=512: truncation cut at d_i - s (cut below capacity) x slack D0.
Tests the cut-point-junk theory: the discontinuity coefficient binom(R-1,cut)
drives the junk; lowering the cut should move/kill the bottom-edge death."""
import json, time, sys, os, io, contextlib, importlib.util
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

def attempt(D0, s, tag):
    pr = [d + D0 for d in floor_prof]
    sig = qpow(R-1); blocks = []
    t0 = time.time()
    for i in range(R):
        d = pr[i]
        ideal = psub(sig, pshift(Em(R-2-i), 1))
        cut = d - s if len(ideal)-1 > d else d   # only lower the cut while overlong
        delta = adm_clamp(ideal[:cut+1], d, 'tozero')
        t = psub(sig, bern_to_poly(delta, d))
        if t and t[0] != 0:
            print(f"{tag}: DIE row {i} const 2^{abs(t[0]).bit_length()} ({time.time()-t0:.0f}s)", flush=True)
            return None
        sig = t[1:] if t else []
        blocks.append(delta)
    if sig != []:
        print(f"{tag}: residual nonzero", flush=True); return None
    adm = all(admissible(blocks[i], pr[i]) for i in range(R))
    idt = epoch_sum(R, blocks, pr) == qpow(R-1)
    print(f"{tag}: CLOSED adm={adm} identity={idt} ({time.time()-t0:.0f}s)", flush=True)
    if adm and idt:
        json.dump({"R": R, "profile": pr, "floor_profile": floor_prof, "D0": D0, "cut_offset": s,
                   "blocks": blocks, "verified": True,
                   "source": f"rule A' tozero, cut=d_i-{s} while residual overlong, d_i=exact_floor+{D0}"},
                  open(f"amm12592_witness_R512_ruleAcut{s}_D0_{D0}_boxeph.json", "w"))
        print(f"wrote amm12592_witness_R512_ruleAcut{s}_D0_{D0}_boxeph.json", flush=True)
    return blocks

for D0, s in [(0, 4), (0, 8), (0, 16), (1, 8), (1, 16), (0, 32)]:
    if attempt(D0, s, f"D0={D0} cut-{s}") is not None:
        break
