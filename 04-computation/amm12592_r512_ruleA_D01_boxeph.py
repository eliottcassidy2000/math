"""R=512 epoch, rule A with slack D0=1, EXACT Fib/Lucas floor profile."""
import json, time, io, contextlib, importlib.util, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
spec = importlib.util.spec_from_file_location("iref", os.path.join(HERE, "amm12592_independent_witness_referee_boxeph.py"))
iref = importlib.util.module_from_spec(spec)
with contextlib.redirect_stdout(io.StringIO()):
    spec.loader.exec_module(iref)
from amm12592_allR_greedy_attractor_rule_boxeph import adm_clamp, Em
from amm12592_allR_family_toolbox_boxeph import *

R, D0 = 512, 1
t0 = time.time()
floor_prof = [iref.floor_gamma_star(R+i) for i in range(R)]
pr = [d + D0 for d in floor_prof]
print(f"exact floor profile computed ({time.time()-t0:.0f}s); d_0={pr[0]} d_last={pr[-1]}", flush=True)
GS = (5979874356654402, 10**16)
proxy = [GS[0]*(R+i)//GS[1] for i in range(R)]
print(f"proxy==exact on m in [512,1023]: {proxy == floor_prof}", flush=True)

sig = qpow(R-1)
blocks = []
ok_run = True
for i in range(R):
    d = pr[i]
    ideal = psub(sig, pshift(Em(R-2-i), 1))
    delta = adm_clamp(ideal[:d+1], d, 'tozero')
    D = bern_to_poly(delta, d)
    t = psub(sig, D)
    if t and t[0] != 0:
        print(f"DIE row {i}: const 2^{abs(t[0]).bit_length()}", flush=True)
        ok_run = False
        break
    sig = t[1:] if t else []
    blocks.append(delta)
    if i % 64 == 0:
        print(f"  row {i} done ({time.time()-t0:.0f}s)", flush=True)
        json.dump({"R": R, "checkpoint_row": i}, open("amm12592_r512_ruleA_ckpt_boxeph.json", "w"))
if ok_run and sig == []:
    print(f"R=512 D0=1: ALL ROWS EMITTED, residual exhausted ({time.time()-t0:.0f}s)", flush=True)
    adm = all(admissible(blocks[i], pr[i]) for i in range(R))
    idt = epoch_sum(R, blocks, pr) == qpow(R-1)
    print(f"exact re-verification: admissible={adm} identity={idt}", flush=True)
    if adm and idt:
        json.dump({"R": R, "profile": pr, "floor_profile": floor_prof, "D0": D0,
                   "blocks": blocks, "verified": True,
                   "source": "greedy attractor rule A tozero, d_i=exact_floor(gamma*(R+i))+1"},
                  open("amm12592_witness_R512_ruleA_D0_1_boxeph.json", "w"))
        print("wrote amm12592_witness_R512_ruleA_D0_1_boxeph.json", flush=True)
elif ok_run:
    print(f"R=512 D0=1: final residual nonzero L1~2^{sum(abs(v) for v in sig).bit_length()}", flush=True)
