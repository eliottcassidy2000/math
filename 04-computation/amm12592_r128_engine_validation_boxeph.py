"""Validation for the R=128 engine before the long hunt:
V1: optimized step()/final_decode() agree with the original R=64 solver's
    on randomized inputs (semantics preserved).
V2: the saved R=64 floor witness's residual trajectory passes the new exact
    capacity + eval-1 prune at every row (prune soundness on a real winner);
    ditto R=8/16/32 witnesses.
V3: timing of one R=128 beam row expansion (calibration).
Exact int arithmetic only."""
import json, os, random, sys, time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r64_floor_solver_boxeph as base
import amm12592_r128_floor_solver_boxeph as eng

# --- V1: randomized agreement of step / final_decode -----------------------
rng = random.Random(20260803)
agree = 0
for trial in range(4000):
    d = rng.randrange(3, 40)
    ln = rng.randrange(1, d + 8)
    sig = [rng.choice([1, -1])] + [rng.randrange(-50, 51)
                                   for _ in range(ln - 1)]
    ctrl = rng.randrange(0, 4)
    tg = tuple(rng.choice([None, -3, -2, -1, 0, 1, 2, 3])
               for _ in range(ctrl))
    r1 = base.step(list(sig), d, tg)
    r2 = eng.step(list(sig), d, tg)
    assert r1 == r2, (sig, d, tg, r1, r2)
    f1 = base.final_decode(list(sig), d)
    f2 = eng.final_decode(list(sig), d)
    assert f1 == f2, (sig, d, f1, f2)
    agree += 1
print(f"V1 PASS: step/final_decode agree on {agree} randomized cases")

# --- V2: prune soundness on saved winners ----------------------------------
def residual_trajectory(R, blocks, d):
    """sigma_i from the recursion p sigma_i = sigma_{i-1} - Delta_i."""
    sig = eng.qpow(R - 1)
    out = []
    for i in range(R):
        de = blocks[i]
        res = list(sig) + [0] * max(0, d[i] + 2 - len(sig))
        for k, v in enumerate(de):
            if v:
                t = eng.qk(k)
                off = d[i] - k
                for s in range(k + 1):
                    res[off + s] -= v * t[s]
        assert res[0] == 0, f"row {i}: residual not divisible by p"
        sig = eng.trim(res[1:])
        out.append(sig)
    assert sig == [], "final residual nonzero"
    return out

def check_witness_prune(R, blocks, d, label):
    caps = eng.coeff_caps(R, d)
    traj = residual_trajectory(R, blocks, d)
    for i in range(R - 1):
        ok = eng.sig_ok(traj[i], caps[i], R - 1 - i)
        assert ok, f"{label}: winner residual at row {i} FAILS prune (BUG)"
    print(f"V2 PASS: {label} winner trajectory passes caps+eval1 prune "
          f"at all {R-1} rows")

with open(os.path.join(HERE, "amm12592_floor_witness_R64.json")) as f:
    w64 = json.load(f)
check_witness_prune(64, w64["blocks"], w64["profile"], "R=64 slim")
with open(os.path.join(HERE, "amm12592_floor_witnesses_R8_R16_R32.json")) as f:
    wsm = json.load(f)
for key, w in (wsm.items() if isinstance(wsm, dict) else []):
    if isinstance(w, dict) and "blocks" in w and "profile" in w:
        check_witness_prune(w["R"], w["blocks"], w["profile"], f"{key}")
with open(os.path.join(HERE, "amm12592_floor_witness_R64_boxeph.json")) as f:
    w64f = json.load(f)
if "blocks" in w64f:
    check_witness_prune(64, w64f["blocks"], w64f["profile"], "R=64 fat")

# --- V3: timing one R=128 row ----------------------------------------------
R = 128
d = eng.prof(R)
print(f"R=128 profile: d[0]={d[0]} d[-1]={d[-1]}")
caps = eng.coeff_caps(R, d)
t0 = time.time()
states, snaps, hist = None, None, None
# expand rows 0..2 with beam 400 to time it
from itertools import product
opts = [None] + list(range(-2, 3))
grids = [tg for tg in product(opts, repeat=2) if tg[0] in (1, -1)]
st = [([], eng.qpow(R - 1))]
for i in range(3):
    nxt = []
    for acc, sig in st:
        for tg in grids:
            r = eng.step(sig, d[i], tg)
            if r is None:
                continue
            de, ns = r
            if not ns or abs(ns[0]) != 1:
                continue
            if not eng.sig_ok(ns, caps[i], R - 1 - i):
                continue
            nxt.append((acc + [de], ns))
    nxt.sort(key=lambda s2: (len(s2[1]), eng.l1(s2[1])))
    st = nxt[:400]
    print(f"V3 row {i}: children={len(nxt)} kept={len(st)} "
          f"t={time.time()-t0:.2f}s")
print(f"V3: 3 rows in {time.time()-t0:.2f}s "
      f"(row cost grows with beam fill; extrapolate x126)")
print("ALL VALIDATION PASSED")
