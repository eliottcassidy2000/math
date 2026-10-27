"""audit_named_classes.py -- brute-force (no nauty) checks of the named classes in THM-4815 (D) and of the
order-5/6/7 examples found by the audit.  Pure Python, permutations for isomorphism.
Run: python3 -I audit_named_classes.py
"""
from itertools import permutations, combinations


def decode(n, s):
    A = [[0] * n for _ in range(n)]
    k = 0
    for i in range(n):
        for j in range(i + 1, n):
            if s[k] == '1':
                A[i][j] = 1
            else:
                A[j][i] = 1
            k += 1
    return A


def scores(A):
    return [sum(r) for r in A]


def ham_paths(A):
    n = len(A)
    res = []

    def rec(path, used):
        if len(path) == n:
            res.append(tuple(path))
            return
        u = path[-1]
        for v in range(n):
            if not used >> v & 1 and A[u][v]:
                path.append(v)
                rec(path, used | 1 << v)
                path.pop()
    for s in range(n):
        rec([s], 1 << s)
    return res


def is_iso(A, B):
    n = len(A)
    sa, sb = scores(A), scores(B)
    if sorted(sa) != sorted(sb):
        return False
    for p in permutations(range(n)):
        if any(sa[i] != sb[p[i]] for i in range(n)):
            continue
        if all(A[i][j] == B[p[i]][p[j]] for i in range(n) for j in range(n) if i != j):
            return True
    return False


def aut_count(A):
    n = len(A)
    sa = scores(A)
    c = 0
    for p in permutations(range(n)):
        if any(sa[i] != sa[p[i]] for i in range(n)):
            continue
        if all(A[i][j] == A[p[i]][p[j]] for i in range(n) for j in range(n) if i != j):
            c += 1
    return c


def converse(A):
    n = len(A)
    return [[A[j][i] if i != j else 0 for j in range(n)] for i in range(n)]


def strong(A, vs):
    vs = list(vs)
    if len(vs) <= 1:
        return True
    if len(vs) == 2:
        return False
    S = set(vs)
    for s in vs:
        seen = {s}
        st = [s]
        while st:
            u = st.pop()
            for v in S:
                if A[u][v] and v not in seen:
                    seen.add(v)
                    st.append(v)
        if seen != S:
            return False
    return True


def kappa(A):
    n = len(A)
    if not strong(A, range(n)):
        return 0
    for k in range(1, n - 1):
        for S in combinations(range(n), k):
            if not strong(A, [v for v in range(n) if v not in S]):
                return k
    return n - 1


def report(n, s, label):
    A = decode(n, s)
    sc = scores(A)
    hp = ham_paths(A)
    mx, mn = max(sc), min(sc)
    ends = sorted(set((p[0], p[-1]) for p in hp))
    minmax = [(a, b) for (a, b) in ends if sc[a] == mn and sc[b] == mx]
    t = sorted(sc)
    L = sum(abs(t[i] - i) for i in range(n)) // 2
    print(f"{label}: {s} n={n} scores={sc} sorted={t} L={L} |Aut|={aut_count(A)} self-converse={is_iso(A, converse(A))} "
          f"strong={strong(A, range(n))} kappa={kappa(A)} #HP={len(hp)}")
    print(f"   Hamiltonian-path endpoint pairs: {ends}")
    print(f"   min->max endpoint pairs: {minmax}")
    return A, hp


if __name__ == '__main__':
    A, hp = report(6, '111011111111101', 'named order-6 class (THM-4815 D)')
    print('   out-neighbourhoods:', {v: [w for w in range(6) if A[v][w]] for v in range(6)})
    with_534 = [p for p in hp if any(p[i:i + 3] == (5, 3, 4) for i in range(4))]
    print(f"   paths containing 5->3->4: {len(with_534)} of {len(hp)}; paths without it: {[p for p in hp if p not in with_534]}")
    B, _ = report(6, '101011111111110', 'second order-6 class with scores 1,1,2,3,4,4, L=1, d=3')
    print('   is it the converse of the named class?', is_iso(B, converse(A)))
    report(5, '1111100111', 'order-5 rigid L-failure (K1 => S4)')
    report(5, '1101111111', 'order-5 rigid L-failure (S4 => K1)')
    report(7, '110011111011111100110', 'order-7 strong 2-connected L-failure')
    report(7, '111110110011101111111', 'order-7 black colour obstruction (scores 1,3^5,5)')
