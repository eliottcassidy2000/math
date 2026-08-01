"""Systematic search over BOUNDED offset profiles e_i at R=128, gamma=3/5.
d_i = floor(3(R+i)/5) + e_i;  T(n)=n+1+floor(3n/5)+e(n) with e bounded gives
C = 8/5 for ANY such e.  The fixed-D0 family is the single point e = const."""
import random, time, json
from perrow import run_e, verify_e
import amm12592_gamma35_beam_deathstar as beam

R = 128
best = (1e9, None)
rng = random.Random(11)
cands = []
# structured families
for base in (2, 3):
    for seg in (8, 16, 32):
        for bump in (1, 2):
            e = [base] * R
            for j in range(R - seg, R): e[j] = base + bump
            for j in range(R - seg // 2, R): e[j] = base + bump + 1
            cands.append((f"base{base} seg{seg} bump{bump} staircase", e))
# ramps
for base in (2, 3):
    for top in (4, 5, 6):
        e = [base + round((top - base) * i / (R - 1)) for i in range(R)]
        cands.append((f"ramp {base}->{top}", e))
# random bounded profiles
for t in range(8):
    e = [rng.randint(2, 5) for _ in range(R)]
    e.sort()
    cands.append((f"random sorted #{t}", e))

for name, e in cands:
    t0 = time.time()
    sol, msg, d = run_e(R, e, bw=2000, tail=48)
    ok = verify_e(R, sol, d) if sol else False
    m = None
    if "best |want|/cap" in msg:
        try: m = float(msg.split("=")[-1])
        except Exception: pass
    if m is not None and m < best[0]: best = (m, name)
    print(f"{name:32s}: {msg}  verify={ok}  ({time.time()-t0:.0f}s)", flush=True)
    if ok:
        print(f"*** R=128 CLOSES with bounded offset '{name}' -> C = 8/5 for n <= 255 ***")
        json.dump({"R": R, "e": e, "sol": sol}, open("sol128.json", "w")); break
print(f"\nBEST MISS over {len(cands)} bounded-e profiles: {best[0]} at '{best[1]}'")
