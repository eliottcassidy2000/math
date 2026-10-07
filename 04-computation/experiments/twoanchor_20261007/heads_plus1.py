#!/usr/bin/env python3
"""Complete list of head collisions at the anchor +1 with debt 3 (mac-mini-2026-10-07-twoanchor):
words u, v with |v| = |u| + 3, sum(u) = sum(v) + 2, F_v(1) = 4 F_u(1) + 1 (THM-4601 (iv)), for sum(u) <= SMAX,
u not starting with 2 (a leading 2 belongs to the two-run), u and v of the same run type (first letter 1 iff r = 1).
F_w(1) = num / 2^sum(w), with num_{wa} = 3 num_w + 2^{sum w}. The source class of a pair (exact letters of u from X)
has measure 2^(r - sum u) relative to the run type {v2(X-1) = r}; the child's word is then forced (THM-4601 (iv)).
Classes of different u are disjoint unless one u is a prefix of the other; the union counts prefix-minimal pairs."""
import sys
from collections import defaultdict
SMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 22
words = defaultdict(list)          # (length, total, num) -> list of words
def gen(w, S, num):
    if w: words[(len(w), S, num)].append(w)
    for a in range(1, SMAX - 1 - S + 1):
        gen(w + (a,), S + a, 3*num + (1 << S))
gen((), 0, 1)
pairs = []
for (L, S, num), us in list(words.items()):
    if S < 3: continue
    # target: F_v(1) = 4 num/2^S + 1 = (4 num + 2^S)/2^S ; v has total S-2: num_v / 2^(S-2) => num_v = (4num + 2^S)/4 = num + 2^(S-2)
    tnum = num + (1 << (S - 2))
    for v in words.get((L + 3, S - 2, tnum), []):
        for u in us:
            if u[0] == 2: continue
            if (u[0] == 1) != (v[0] == 1): continue
            pairs.append((S, u, v))
pairs.sort()
minimal = []
for S, u, v in pairs:
    if not any(u[:len(m[1])] == m[1] for m in minimal): minimal.append((S, u, v))
print(f"head collisions at +1 with debt 3, sum(u) <= {SMAX}: {len(pairs)} pairs, {len(minimal)} prefix-minimal sources")
cov = {1: 0.0, 2: 0.0}
for S, u, v in minimal:
    r = 1 if u[0] == 1 else 2
    cov[r] += 2.0 ** (r - S)
for S, u, v in pairs[:25]:
    print(f"  sum(u)={S:2d}  u={u}  v={v}")
print(f"union measure of source classes (relative to run type): r=1 {cov[1]:.6f}, r=2 {cov[2]:.6f}")
