"""MISTAKE-683 side claims: (a) the order-6 strong class 111011111111101 (gentourng upper-triangle code) is rigid,
not self-converse, scores 1,1,2,3,4,4, L = 1; (b) S_4 => S_5' (S_5' a rigid prime 5-tournament) is rigid,
non-self-converse and has no eligible vertex for delta0 (THM-4818 (C2): hence a delta0 failure)."""
from itertools import permutations, product

def from_code(code, n):
    arcs = set(); b = 0
    for i in range(n):
        for j in range(i + 1, n):
            arcs.add((i, j) if code[b] == '1' else (j, i)); b += 1
    return arcs

def autanti(arcs, n):
    aut = anti = 0
    for p in permutations(range(n)):
        if all(((p[i], p[j]) in arcs) == ((i, j) in arcs) for i in range(n) for j in range(n) if i != j): aut += 1
        if all(((p[j], p[i]) in arcs) == ((i, j) in arcs) for i in range(n) for j in range(n) if i != j): anti += 1
    return aut, anti

def scores(arcs, n): return [sum((i, j) in arcs for j in range(n)) for i in range(n)]

def is_prime(arcs, n):
    out = [{j for j in range(n) if (i, j) in arcs} for i in range(n)]
    for u in range(n):
        for v in range(u + 1, n):
            M = {u, v}
            while True:
                add = {w for w in range(n) if w not in M and 0 < len(out[w] & M) < len(M)}
                if not add: break
                M |= add
            if len(M) < n: return False
    return True

A = from_code('111011111111101', 6)
s = scores(A, 6); L = sum(abs(a - i) for i, a in enumerate(sorted(s))) / 2
print("(a) order-6 class: scores", sorted(s), "(|Aut|, #anti) =", autanti(A, 6), "L =", L)

# rigid prime 5-tournaments
rig5 = []
seen = set()
for bits in product('01', repeat=10):
    T = from_code(''.join(bits), 5)
    if not is_prime(T, 5): continue
    a, an = autanti(T, 5)
    key = min(tuple(sorted((p[i], p[j]) for (i, j) in T)) for p in permutations(range(5)))
    if key in seen: continue
    seen.add(key)
    if a == 1: rig5.append(T)
S4 = {(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)}
for T5 in rig5:
    arcs = set(S4) | {(a + 4, b + 4) for (a, b) in T5} | {(a, b) for a in range(4) for b in range(4, 9)}
    n = 9
    sc = scores(arcs, n)
    d0 = min(min(x, n - 1 - x) for x in sc)
    top, bot = range(0, 4), range(4, 9)          # initial (top) and terminal (bottom) strong components
    eligible = [w for w in top if sc[w] == d0] + [w for w in bot if n - 1 - sc[w] == d0]
    aut, anti = autanti(arcs, n)
    print("(b) S_4 => S_5': S_5' scores", sorted(scores(T5, 5)), "| scores", sc, "delta0 =", d0, "eligible:", eligible, "|Aut| =", aut, "#anti =", anti)
