"""
audit_A_lattice.py -- independent recount (audit A) of the material-signature lattice (section 9 of the note)
and of the per-file pawn-move lemma (Lemma 3.1) with an own exhaustive model; the 118 formula.
Run: python3 audit_A_lattice.py
"""
import itertools, collections, functools

# ---- lattice without promotion ----
S0 = [(p, n, bl, bd, r, q) for p in range(9) for n in range(3) for bl in range(2) for bd in range(2) for r in range(3) for q in range(2)]
inv0 = {s[1:] for s in S0}
cap0 = 2 * sum(sum(1 for c in s if c > 0) for s in S0) * len(S0)
print(f'no promotion: piece inventories {len(inv0)}; per-side signatures {len(S0)}; joint {len(S0)**2}; capture edges {cap0}')

# ---- with promotion ----
def promoted(s):
    p, n, bl, bd, r, q = s
    return max(0, n - 2) + max(0, bl - 1) + max(0, bd - 1) + max(0, r - 2) + max(0, q - 1)

S1 = [s for s in itertools.product(range(9), range(11), range(10), range(10), range(11), range(10)) if promoted(s) <= 8 - s[0]]
assert all(sum(s) <= 15 for s in S1)
byp = collections.Counter(s[0] for s in S1)
S1set = set(S1)
cap1 = 2 * sum(sum(1 for c in s if c > 0) for s in S1) * len(S1)
promo = 0
for s in S1:
    p = s[0]
    if p == 0:
        continue
    for t in range(1, 6):
        u = list(s); u[0] -= 1; u[t] += 1
        if tuple(u) in S1set:
            promo += 1
print(f'with promotion: per-side {len(S1)}; by pawn count {dict(sorted(byp.items()))}; pawnless inventories {byp[0]}; '
      f'joint {len(S1)**2}; capture edges {cap1}; promotion edges {2 * promo * len(S1)}')
print('expected from note: 72 / 648 / 419,904 / 3,125,952 ; 26,298 / 9,160 / 691,584,804 / 6,129,222,264 / 4,506,951,240')

# ---- per-file lemma, own model: longest-path DP over (white rank or None, black rank or None) ----
# white pawn starts on rank 2 moving up (1 step, or 2 from rank 2), black starts on rank 7 moving down;
# neither pawn captures (so both stay on the file); they never pass or share a square; at any moment
# either pawn (or both, in sequence) may be removed by an outside capture.
@functools.lru_cache(None)
def longest(w, b, captures_allowed):
    """max number of further pawn moves on this file from state (w, b); captures_allowed = how many
    outside captures may still happen (0, 1 or 2)."""
    best = 0
    if w is not None and w < 8:
        for st in ([1, 2] if w == 2 else [1]):
            nw = w + st
            if b is not None and nw >= b:
                continue
            best = max(best, 1 + longest(nw, b, captures_allowed))
    if b is not None and b > 1:
        for st in ([1, 2] if b == 7 else [1]):
            nb = b - st
            if w is not None and nb <= w:
                continue
            best = max(best, 1 + longest(w, nb, captures_allowed))
    if captures_allowed > 0:
        if w is not None:
            best = max(best, longest(None, b, captures_allowed - 1))
        if b is not None:
            best = max(best, longest(w, None, captures_allowed - 1))
    return best

print(f'per-file lemma: max pawn moves with no capture of either pawn: {longest(2, 7, 0)}; '
      f'with at most one pawn captured: {longest(2, 7, 1)}; with both allowed to be captured: {longest(2, 7, 2)}')
# the 118 formula
vals = [10 * (8 - k) + 12 * k + (30 - k) for k in range(9)]
print('110 + k for k = 0..8:', vals, '; max', max(vals))
# cross-check against the classical count 16*6 + 30 - 8
print('classical count 16*6 + 30 - 8 =', 16 * 6 + 30 - 8)
