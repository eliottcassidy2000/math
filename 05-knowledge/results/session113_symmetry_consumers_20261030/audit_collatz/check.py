#!/usr/bin/env python3
# Verification gates are explicit and remain active under python -O.
"""Independent scalar-integer check; does not import producer or NumPy auditor."""
import argparse, struct, hashlib, json, math, time
from pathlib import Path
p = argparse.ArgumentParser()
p.add_argument('E', type=int)
p.add_argument('vector')
p.add_argument('--clock', type=int, default=18)
p.add_argument('--out')
args = p.parse_args()
start = time.time()
E = args.E
q = 5 ** E
D = 8 * q
N = 5 * q - 5
cap = 64
raw = Path(args.vector).read_bytes()
if not len(raw) == 16 * q:
    raise AssertionError()
V = struct.unpack('<' + 'q' * (2 * q), raw)

def height(n):
    if n == 0:
        return E
    ans = 0
    while ans < E and n % 5 == 0:
        ans += 1
        n //= 5
    return ans
hs = [height(z) for z in range(q)]
inv = [pow(2 ** a, -1, q) for a in range(cap + 1)]
ks = [height(2 ** a - 2) for a in range(cap + 1)]
minimum_slack = None
for a in range(1, cap + 1):
    ia = inv[a]
    for flag in range(2):
        for z in range(q):
            f = (flag + a + 1) % 2
            destz = (3 * z + 2 ** a - 2) * ia % q
            s = flag * q + z
            t = f * q + destz
            if a == 1:
                reward_slope = sum((D // (2 * 5 ** (e - 1)) for e in range(1, hs[z] + 1))) if flag == 0 else 0
                cost = reward_slope - N
            else:
                cost = (D * max(hs[destz], ks[a]) if f == 0 else 0) - N * a
            slack = V[t] - V[s] - cost
            if not slack >= 0:
                raise AssertionError((a, flag, z, slack))
            minimum_slack = slack if minimum_slack is None else min(minimum_slack, slack)
            if a > 1 and s == 0:
                if not V[t] >= cost:
                    raise AssertionError(('entry', a))
mn, mx = (min(V), max(V))
if not mx <= 0:
    raise AssertionError()
if not mx + D * E - N * (cap + 1) <= mn:
    raise AssertionError()
if not D * E - N * (cap + 1) <= mn:
    raise AssertionError()
count = 0
largest_excess = None
stack = [(1, 1, 0, 0, 0, 0, 0, False)]
while stack:
    P, Q, B, A, S, C, flag, entered = stack.pop()
    for a in range(1, args.clock - A + 1):
        pp, qq, bb = (3 * P, 2 ** a * Q, 3 * B + Q)
        aa = A + a
        ff = (flag + a + 1) % 2
        gg = math.gcd(abs(pp - qq), bb)
        ss = S + height(gg)
        if not D * ss <= N * aa:
            raise AssertionError(('true bound', pp, qq, bb, aa))
        zz = (bb + qq - pp) * pow(qq, -1, q) % q
        ent = entered or a > 1
        if a > 1:
            cc = C + (D * max(height(bb + qq - pp), height(2 ** a - 2)) if ff == 0 else 0) - N * a
        elif entered:
            cc = C + (sum((D // (2 * 5 ** (e - 1)) for e in range(1, height(bb + qq - pp) + 1))) if ff == 0 else 0) - N
        else:
            cc = 0
        if not D * ss - N * aa <= cc:
            raise AssertionError(('envelope', aa, ss, cc))
        if ent:
            if not cc <= V[ff * q + zz]:
                raise AssertionError(('path', aa, cc))
        count += 1
        if aa < args.clock:
            stack.append((pp, qq, bb, aa, ss, cc, ff, ent))
r = dict(E=E, states=2 * q, edges=2 * q * cap, minimum_slack=minimum_slack, minimum_V=mn, maximum_V=mx, entry_tail_margin=mn - (D * E - N * (cap + 1)), edge_tail_margin=mn - (mx + D * E - N * (cap + 1)), exact_positive_words=count, maximum_word_clock=args.clock, sha256=hashlib.sha256(raw).hexdigest(), seconds=time.time() - start)
s = json.dumps(r, indent=2)
print(s)
if args.out:
    Path(args.out).write_text(s + '\n')
