"""Audit THM-4824 (E): parity-coset shell lemma on Z[w] (basis 1, w; Q = a^2 - ab + b^2).
Direct slack from THM-4820 (B): slack = beta_+ - floor((r + s + beta)/3), r = Q(y) mod 3, s = Q(z) mod 3.
Compares, for every x with Q(x) <= QMAX, the set of tight nontrivial splits found by brute force over y with the set
predicted by the shell criterion; checks Q(v) = Q(x) - 2 beta and the bijection; runs the 11 and 51 seed tests."""
import sys
Q = lambda a, b: a*a - a*b + b*b
QMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 700
R = int(2 * (QMAX + 3) ** 0.5) + 3
pts = [(a, b) for a in range(-R, R + 1) for b in range(-R, R + 1)]
byQ = {}
for p in pts: byQ.setdefault(Q(*p), []).append(p)
mism = 0; nx = 0; seeds = []
for x in pts:
    qx = Q(*x)
    if qx > QMAX or qx < 2: continue
    nx += 1
    brute = set()
    for y in pts:
        z = (x[0] - y[0], x[1] - y[1])
        if y == (0, 0) or z == (0, 0): continue
        qy, qz = Q(*y), Q(*z)
        if qy + qz > qx + 2: continue           # tight needs beta >= -2
        beta = qx - qy - qz
        slack = max(beta, 0) - (qy % 3 + qz % 3 + beta) // 3
        assert slack >= 0
        if slack == 0: brute.add((y, z))
    shells = [qx, qx + 2] if qx % 3 == 1 else [qx - 2, qx, qx + 2, qx + 4]
    pred = set()
    for qv in shells:
        for v in byQ.get(qv, []):
            if (v[0] - x[0]) % 2 or (v[1] - x[1]) % 2: continue
            if v == x or v == (-x[0], -x[1]): continue
            y = ((x[0] + v[0]) // 2, (x[1] + v[1]) // 2); z = ((x[0] - v[0]) // 2, (x[1] - v[1]) // 2)
            beta = qx - Q(*y) - Q(*z)
            assert Q(*v) == qx - 2 * beta
            pred.add((y, z))
    if brute != pred: mism += 1
    if not brute: seeds.append(x)
print(f"x with 2 <= Q(x) <= {QMAX}: {nx}; mismatches between brute-force tight splits and the shell criterion: {mism}")
rat = sorted({abs(a) for (a, b) in seeds if b == 0})
print("rational-integer seeds (no tight split) with n^2 <= QMAX:", rat)
print("seed norms <= QMAX:", sorted({Q(*s) for s in seeds}))
def seedtest(n):
    qx = n*n
    shells = [qx, qx + 2] if qx % 3 == 1 else [qx - 2, qx, qx + 2, qx + 4]
    out = {}
    Rn = 2*n + 4
    for qv in shells:
        out[qv] = sorted((a, b) for a in range(-Rn, Rn + 1) for b in range(-Rn, Rn + 1) if Q(a, b) == qv and (a - n) % 2 == 0 and b % 2 == 0)
    return out
for n in (11, 51):
    print(f"seed test for {n}: coset points on shells:", seedtest(n))
