#!/usr/bin/env python3
"""Audit B (independent) Monte Carlo for HYP-9241 (Karlin-McGregor / vicious-walker law for translation partners).

For random y (BITS-bit, top bit set; parities exact for TMAX << BITS Terras steps) follow the Terras orbits of
y, y+1, y+2, y+3 up to TMAX and record, for each of the six pairs, the first EQUAL-time coincidence (merge).
A global value -> (class, time) dictionary also detects every UNEQUAL-time coincidence between any two orbits.
Outputs q_1, q_2, q_3 (no merge among y..y+R by time T), every pair survival q_ij, the ratios
q_2/q_1^3, q_3/q_1^6 and the pair-independence ratios q_2/(q01 q12 q02), q_3/prod(6 pairs), local slopes,
first-merge statistics for the triple, the number of meetings (equal odd counts) before a merge,
and the parity-difference frequencies of the three walkers while all three are apart.
usage: km_audit.py N BITS TMAX SEED
"""
import random, math, sys, time, json

N = int(sys.argv[1]) if len(sys.argv) > 1 else 5000
BITS = int(sys.argv[2]) if len(sys.argv) > 2 else 2048
TMAX = int(sys.argv[3]) if len(sys.argv) > 3 else 1024
SEED = int(sys.argv[4]) if len(sys.argv) > 4 else 20261007
assert TMAX + 16 <= BITS          # parities exact (Haar) and values stay >= 2^(BITS-1-TMAX)
rnd = random.Random(SEED)
PAIRS = [(0, 1), (1, 2), (2, 3), (0, 2), (1, 3), (0, 3)]
NEVER = TMAX + 1
merge_times = []          # per sample: dict pair -> first equal-time merge time (NEVER if none by TMAX)
unequal = 0               # unequal-time coincidences found
first_triple = {'01': 0, '12': 0, '02': 0, 'simul': 0, 'none': 0}
meet_before_merge = []    # pair (0,1): number of arrivals at equal odd-step count before its merge
par_counts = [0]*8        # parity patterns (b0,b1,b2) while 0,1,2 pairwise apart, t >= 64
t0 = time.time()
for it in range(N):
    y = rnd.getrandbits(BITS - 1) | (1 << (BITS - 1))
    x = [y, y + 1, y + 2, y + 3]
    cls = [0, 1, 2, 3]                     # class label of each orbit
    seen = {}
    for c in range(4): seen[x[c]] = (c, 0)
    mt = {pq: NEVER for pq in PAIRS}
    odd = [0, 0, 0, 0]
    meets01, prev_eq01 = 0, True
    for t in range(1, TMAX + 1):
        reps = sorted(set(cls))
        if len(reps) == 1: break
        if t >= 64 and cls[0] != cls[1] and cls[1] != cls[2] and cls[0] != cls[2]:
            par_counts[(x[0] & 1) | ((x[1] & 1) << 1) | ((x[2] & 1) << 2)] += 1
        newval = {}
        for c in reps:
            v = x[cls.index(c)]
            if v & 1:
                v = (3*v + 1) >> 1
                for i in range(4):
                    if cls[i] == c: odd[i] += 1
            else:
                v >>= 1
            newval[c] = v
        for i in range(4): x[i] = newval[cls[i]]
        # coincidence detection
        merged_now = []
        for c in reps:
            v = newval[c]
            if v in seen:
                c2, s = seen[v]
                if s == t:
                    merged_now.append((c2, c))
                else:
                    unequal += 1
            else:
                seen[v] = (c, t)
        for (a, b) in merged_now:          # a = first class to reach the value at time t (never relabelled)
            for i in range(4):
                if cls[i] == b: cls[i] = a
        for (i, j) in PAIRS:
            if mt[(i, j)] == NEVER and cls[i] == cls[j]: mt[(i, j)] = t
        # meetings of the (0,1) walkers (equal odd counts) while unmerged
        if cls[0] != cls[1]:
            eq = (odd[0] == odd[1])
            if eq and not prev_eq01: meets01 += 1
            prev_eq01 = eq
    merge_times.append(mt)
    if mt[(0, 1)] <= TMAX: meet_before_merge.append(meets01)
    a, b, c = mt[(0, 1)], mt[(1, 2)], mt[(0, 2)]
    m = min(a, b, c)
    if m == NEVER: first_triple['none'] += 1
    elif [a, b, c].count(m) > 1: first_triple['simul'] += 1
    elif m == a: first_triple['01'] += 1
    elif m == b: first_triple['12'] += 1
    else: first_triple['02'] += 1
elapsed = time.time() - t0

checkpoints = [c for c in (16, 32, 64, 128, 256, 512, 1024, 2048, 4096) if c <= TMAX]
def surv(pairs, T):
    return sum(1 for mt in merge_times if all(mt[pq] > T for pq in pairs)) / N
G = {1: [(0, 1)], 2: [(0, 1), (1, 2), (0, 2)], 3: PAIRS}
lines = [f"km_audit: N={N} BITS={BITS} TMAX={TMAX} SEED={SEED} ({elapsed:.0f}s)"]
lines.append(f"unequal-time coincidences between any two of the 4 orbits: {unequal}")
tab = []
for T in checkpoints:
    q = {R: surv(G[R], T) for R in (1, 2, 3)}
    pq = {pq: surv([pq], T) for pq in PAIRS}
    se = {R: math.sqrt(q[R]*(1 - q[R])/N) for R in (1, 2, 3)}
    r2 = q[2]/q[1]**3 if q[1] > 0 else float('nan')
    r3 = q[3]/q[1]**6 if q[1] > 0 else float('nan')
    ind2 = q[2]/(pq[(0, 1)]*pq[(1, 2)]*pq[(0, 2)]) if q[2] > 0 else float('nan')
    ind3 = q[3]/math.prod(pq[x] for x in PAIRS) if q[3] > 0 else float('nan')
    tab.append(dict(T=T, q1=q[1], q2=q[2], q3=q[3], se1=se[1], se2=se[2], se3=se[3], ev2=round(q[2]*N), ev3=round(q[3]*N),
                    q01=pq[(0, 1)], q12=pq[(1, 2)], q23=pq[(2, 3)], q02=pq[(0, 2)], q13=pq[(1, 3)], q03=pq[(0, 3)],
                    r2=r2, r3=r3, ind2=ind2, ind3=ind3))
lines.append("T | q1 q2 q3 (events q2,q3) | q2/q1^3 [rel.err] | q3/q1^6 [rel.err] | q2/(q01 q12 q02) | q3/prod6 | lag-1 q01,q12,q23 | lag-2 q02,q13 | lag-3 q03")
for r in tab:
    e2 = math.sqrt((r['se2']/r['q2'])**2 + 9*(r['se1']/r['q1'])**2) if r['q2'] > 0 else float('nan')
    e3 = math.sqrt((r['se3']/r['q3'])**2 + 36*(r['se1']/r['q1'])**2) if r['q3'] > 0 else float('nan')
    lines.append(f"{r['T']:5d} | {r['q1']:.4f} {r['q2']:.4f} {r['q3']:.5f} ({r['ev2']},{r['ev3']}) | {r['r2']:.3f} [{e2:.3f}] | {r['r3']:.3f} [{e3:.3f}]"
                 f" | {r['ind2']:.3f} | {r['ind3']:.3f} | {r['q01']:.4f},{r['q12']:.4f},{r['q23']:.4f} | {r['q02']:.4f},{r['q13']:.4f} | {r['q03']:.4f}")
for R in (1, 2, 3):
    sl = []
    for a, b in zip(tab, tab[1:]):
        if b[f'q{R}'] > 0 and b[f'q{R}']*N >= 20:
            sl.append(round(-math.log(b[f'q{R}']/a[f'q{R}'])/math.log(b['T']/a['T']), 3))
    lines.append(f"R={R} local slopes (between consecutive checkpoints, >=20 events): {sl}")
# global slope ratios over the window [16, T*] where q3 has >= 50 events
good = [r for r in tab if r['q3']*N >= 50]
if len(good) >= 2:
    a, b = good[0], good[-1]
    s = {R: -math.log(b[f'q{R}']/a[f'q{R}'])/math.log(b['T']/a['T']) for R in (1, 2, 3)}
    lines.append(f"mean slopes on [{a['T']},{b['T']}]: {[round(s[R], 3) for R in (1, 2, 3)]}; ratio 1 : {s[2]/s[1]:.2f} : {s[3]/s[1]:.2f}")
lines.append(f"first merge within the triple (y,y+1,y+2): {first_triple}  (02 = the NON-adjacent pair merges strictly first)")
if meet_before_merge:
    lines.append(f"pair (y,y+1): mean number of returns to equal odd counts before the merge = "
                 f"{sum(meet_before_merge)/len(meet_before_merge):.2f} over {len(meet_before_merge)} merged samples")
tot = sum(par_counts)
if tot:
    P = lambda f: sum(par_counts[m] for m in range(8) if f(m & 1, (m >> 1) & 1, (m >> 2) & 1))/tot
    f01, f02, f12 = P(lambda a, b, c: a != b), P(lambda a, b, c: a != c), P(lambda a, b, c: b != c)
    g = P(lambda a, b, c: a != b and a != c)
    # wedge angles of the non-collision chambers after whitening the increment covariance of (k01, k02)
    import numpy as np
    S = np.array([[f01, g], [g, f02]]); Si = np.linalg.inv(S)
    dirs = [np.array([0.0, 1.0]), np.array([1.0, 0.0]), np.array([1.0, 1.0])]
    W = np.linalg.cholesky(Si).T
    angs = sorted((math.atan2(*(W @ d)[::-1]) % math.pi) for d in dirs)
    wedges = [angs[1] - angs[0], angs[2] - angs[1], math.pi - (angs[2] - angs[0])]
    lines.append(f"walker parity statistics while (y,y+1,y+2) pairwise apart, t>=64 ({tot} steps): P(flip 01)={f01:.3f}, "
                 f"P(flip 02)={f02:.3f}, P(flip 12)={f12:.3f}, P(01 and 02 flip)={g:.3f}; whitened chamber angles/pi = "
                 f"{[round(w/math.pi, 3) for w in wedges]} (vicious walkers: 1/3 each -> exponent 3/2)")
out = "\n".join(lines)
print(out)
with open(f"km_audit_N{N}_B{BITS}_T{TMAX}_s{SEED}.out", "w") as f: f.write(out + "\n")
with open(f"km_audit_N{N}_B{BITS}_T{TMAX}_s{SEED}.json", "w") as f: json.dump(tab, f, indent=0)
