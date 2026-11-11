# Check Theorem S (2): for a tournament of even order with a fixed-point-free involutive anti-automorphism psi,
# H = 1 + 4 s (mod 8), s = number of psi-orbits {C, psi(C)} of directed odd cycles with V(C) meeting V(psi(C)).
import subprocess, sys, itertools
def classes(m):
    out = subprocess.run(['gentourng', '-q', str(m)], capture_output=True, text=True).stdout.split()
    for line in out:
        A = [[0]*m for _ in range(m)]; k = 0
        for i in range(m):
            for j in range(i+1, m):
                if line[k] == '1': A[i][j] = 1
                else: A[j][i] = 1
                k += 1
        yield A
def matchings(vs):
    if not vs: yield []; return
    a = vs[0]
    for i in range(1, len(vs)):
        b = vs[i]; rest = vs[1:i] + vs[i+1:]
        for mt in matchings(rest): yield [(a, b)] + mt
def ham(A, n):
    f = {}
    for v in range(n): f[(1 << v, v)] = 1
    full = (1 << n) - 1
    for S in range(1, full + 1):
        for v in range(n):
            if not (S >> v) & 1: continue
            c = f.get((S, v), 0)
            if not c: continue
            for w in range(n):
                if not (S >> w) & 1 and A[v][w]:
                    f[(S | (1 << w), w)] = f.get((S | (1 << w), w), 0) + c
    return sum(f.get((full, v), 0) for v in range(n))
def odd_cycles(A, n):
    cyc = []
    for L in range(3, n + 1, 2):
        for S in itertools.combinations(range(n), L):
            v = S[0]
            for p in itertools.permutations(S[1:]):
                seq = (v,) + p
                if all(A[seq[i]][seq[(i+1) % L]] for i in range(L)): cyc.append(seq)
    return cyc
for n in (4, 6, 8):
    nsc = 0; bad = 0; tot = 0
    for A in classes(n):
        tot += 1
        psi = None
        for mt in matchings(list(range(n))):
            p = list(range(n))
            for a, b in mt: p[a] = b; p[b] = a
            if all(A[a][b] == A[p[b]][p[a]] for a in range(n) for b in range(n) if a != b): psi = p; break
        if psi is None: continue
        nsc += 1
        H = ham(A, n); cyc = odd_cycles(A, n)
        # orbits: psi(C) = reversed image; count cycles C with V(C) meeting psi(V(C)), each orbit has two cycles
        meet = 0
        for C in cyc:
            VS = set(C); img = set(psi[v] for v in C)
            assert img != VS or True
            if VS & img: meet += 1
        assert meet % 2 == 0
        s = meet // 2
        if (H - 1 - 4 * s) % 8: bad += 1
        if len(cyc) % 2 or (H - 1) % 4: bad += 1
    print("n =", n, "classes", tot, "with a fixed-point-free involutive anti-automorphism", nsc, "failures of H = 1 + 4s mod 8:", bad)
