"""Extract the late-row residual attractor from the R=64 floor witnesses:
exact residual trajectory sigma_i for rows ~45..63 of both the slim (direct
beam) and fat (doubling) witnesses; also test, for each row i and each
attractor shape tau, whether Delta = sigma_{i-1} - p*tau is Bernstein-d_i
admissible ("hook" feasibility on the actual winner) and whether the ride
Delta = E-step decodes at the profile degrees of R=128's tail."""
import json, os, sys
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import amm12592_r128_floor_solver_boxeph as eng

def residual_traj(R, blocks, d):
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
        assert res[0] == 0
        sig = eng.trim(res[1:])
        out.append(sig)
    assert sig == []
    return out

for tag, path in (("slim", "amm12592_floor_witness_R64.json"),
                  ("fat", "amm12592_floor_witness_R64_boxeph.json")):
    with open(os.path.join(HERE, path)) as f:
        w = json.load(f)
    R, d, blocks = w["R"], w["profile"], w["blocks"]
    traj = residual_traj(R, blocks, d)
    print(f"== {tag} witness late residuals ==")
    for i in range(48, R):
        sig = traj[i]
        if i >= R - 1:
            print(f"  row {i}: sigma = {sig}")
        elif len(sig) <= 14:
            print(f"  row {i}: deg={len(sig)-1} sigma={sig}")
        else:
            print(f"  row {i}: deg={len(sig)-1} L1={eng.l1(sig)} "
                  f"head={sig[:6]} tail={sig[-6:]}")

# The ride blocks: if sigma_{i-1} = tau_{m+1} and sigma_i = tau_m, the row
# block is Delta_i = tau_{m+1} - p tau_m.  Print those for the slim witness
# to identify the traveling shape.
with open(os.path.join(HERE, "amm12592_floor_witness_R64.json")) as f:
    w = json.load(f)
R, d, blocks = w["R"], w["profile"], w["blocks"]
traj = residual_traj(R, blocks, d)
print("== slim ride blocks Delta_i = sigma_{i-1} - p sigma_i (rows 54..63) ==")
for i in range(54, 64):
    prev = traj[i - 1] if i > 0 else eng.qpow(R - 1)
    cur = traj[i] if i < R else []
    de = blocks[i]
    # reconstruct Delta_i as polynomial
    poly = [0] * (d[i] + 1)
    for k, v in enumerate(de):
        if v:
            t = eng.qk(k)
            off = d[i] - k
            for s in range(k + 1):
                poly[off + s] += v * t[s]
    print(f"  row {i}: d={d[i]} Delta={eng.trim(poly)} "
          f"nonzero cells={[(k, v) for k, v in enumerate(de) if v]}")
