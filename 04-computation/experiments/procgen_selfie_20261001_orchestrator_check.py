"""procgen_selfie_20261001_orchestrator_check.py -- the orchestrator's independent audit of the selfie lane, 2026-10-01.

Part I (C engine procgen_selfie_20261001_orchestrator_check.c, compiled into a temp dir; needs gentourng):
  labeled N = 3..7: HP-covered counts, all-odd/all-even counts, min/max #odd arcs, dead-arc formula,
  "every arc on some HP and some arc on every HP" impossible, sum_L H(switch_L T) = 2 N!,
  sum_L (-1)^|L| H(switch_L T) = 0, the loop-Walsh degree bound on every tiling, constant-H classes,
  the x = 0 fixed-point OCF identity, selfie Redei (N <= 6);
  classes N = 3..9 (gentourng): the census table; beta = Hall bound for every class N <= 7;
  Redei shaving (single-arc deletions with c(e) even, down to H = 1) for every class N <= 7;
  beta = N - 1 for every regular class N = 3, 5, 7;
  QR_q and QR_q minus a vertex (q = 7, 11, 19 exact; 23 mod 2), circulants Z_3..Z_13, Cayley Z_3 x Z_3,
  the two all-odd N = 10 classes, the four 5-odd-arc N = 10 witnesses, the N = 7 sigma counterexample.

Part II: Theorem E1 (OPEN-Q-060), as follows.

Claim: A049313(n) = number of unlabeled Euler graphs F on n vertices such that every automorphism g of F
reverses an even number of edges of a fixed orientation O (O: i -> j for i < j).

Checks (written from the definitions; the lane's code was not read):
 (1) per-permutation identity  #Fix(g on labeled switching classes of tournaments) = sum_{F Euler, gF=F} eps_F(g),
     eps from the ORIENTATION definition; every permutation for n <= 5, one per cycle type for n = 6, 7;
 (2) Burnside: (1/n!) sum_g #Fix(g) = A049313(n) for n <= 7 (OEIS values);
 (3) direct orbit count of 'untwisted' Euler graphs (Aut computed by brute force) for n <= 6;
 (4) cycle-form of eps (antipodal pairs of even cycles) agrees with the orientation form on every
     (g, F) checked; Burnside with the cycle form for n = 8, 9 via F_2 linear algebra vs OEIS.
"""
import itertools, math, sys
from collections import Counter

A049313 = [1, 1, 1, 2, 2, 6, 12, 79, 792, 19576, 886288, 75369960]   # n = 1..12 (OEIS)
A002854 = [1, 1, 2, 3, 7, 16, 54, 243, 2038, 33120]                     # n = 1..10

def pairs(n):
    return [(i, j) for i in range(n) for j in range(i + 1, n)]

def euler_graphs(n):
    P = pairs(n); out = []
    for m in range(1 << len(P)):
        deg = [0] * n
        for k, (i, j) in enumerate(P):
            if m >> k & 1:
                deg[i] ^= 1; deg[j] ^= 1
        if not any(deg):
            out.append(frozenset(P[k] for k in range(len(P)) if m >> k & 1))
    return out

def eps_orient(F, g):
    """(-1)^#{edges {i<j} of F with g(i) > g(j)}  (g carries O-orientation against O)."""
    s = 0
    for (i, j) in F:
        if g[i] > g[j]:
            s ^= 1
    return -1 if s else 1

def eps_cycle(F, g):
    """(-1)^#{even cycles of g (length 2k) whose antipodal pairs {x, g^k x} are edges of F}."""
    n = len(g); seen = [False] * n; s = 0
    for x in range(n):
        if seen[x]:
            continue
        cyc = []; y = x
        while not seen[y]:
            seen[y] = True; cyc.append(y); y = g[y]
        L = len(cyc)
        if L % 2 == 0:
            a, b = cyc[0], cyc[L // 2]
            if (min(a, b), max(a, b)) in F:
                s ^= 1
    return -1 if s else 1

def apply_graph(F, g):
    return frozenset((min(g[i], g[j]), max(g[i], g[j])) for (i, j) in F)

# --- switching classes of tournaments: tournament = dict pair -> True if i->j (i<j)
def tiling_rep(T, n):
    """unique switching-equivalent tournament with base path i+1 -> i for all i (THM-474 gauge)."""
    inL = [0] * n
    for i in range(n - 1):
        # pair (i, i+1): T[(i,i+1)] True means i -> i+1, which must be reversed
        inL[i + 1] = inL[i] ^ (1 if T[(i, i + 1)] else 0)
    return tuple(T[(i, j)] ^ bool(inL[i] ^ inL[j]) for (i, j) in pairs(n))

def act_tour(T, g, n):
    R = {}
    for (i, j), fwd in T.items():
        a, b = (g[i], g[j]) if fwd else (g[j], g[i])     # arc a -> b
        if a < b:
            R[(a, b)] = True
        else:
            R[(b, a)] = False
    return R

def fixed_classes(n, g, reps):
    P = pairs(n); cnt = 0
    for t in reps:
        T = dict(zip(P, t))
        if tiling_rep(act_tour(T, g, n), n) == t:
            cnt += 1
    return cnt

def tiling_reps(n):
    P = pairs(n); res = []
    base = [k for k, (i, j) in enumerate(P) if j == i + 1]
    for m in range(1 << len(P)):
        if any(m >> k & 1 for k in base):
            continue
        res.append(tuple(bool(m >> k & 1) for k in range(len(P))))
    return res

def cycle_type_reps(n):
    seen = {}
    for part in partitions(n):
        g = []; start = 0
        for L in part:
            g += [start + (k + 1) % L for k in range(L)]
            start += L
        seen[part] = tuple(g)
    return seen

def partitions(n, maxp=None):
    if maxp is None:
        maxp = n
    if n == 0:
        yield ()
        return
    for p in range(min(n, maxp), 0, -1):
        for rest in partitions(n - p, p):
            yield (p,) + rest

def class_size(part, n):
    c = Counter(part); d = 1
    for L, m in c.items():
        d *= L ** m * math.factorial(m)
    return math.factorial(n) // d

def euler_main():
    ok = True
    for n in range(1, 8):
        E = euler_graphs(n)
        assert len(E) == 2 ** ((n - 1) * (n - 2) // 2)
        reps = tiling_reps(n)
        assert len(reps) == 2 ** ((n - 1) * (n - 2) // 2)
        perms = list(itertools.permutations(range(n))) if n <= 5 else list(cycle_type_reps(n).values())
        mism = 0; cyc_mism = 0
        for g in perms:
            fx = fixed_classes(n, g, reps)
            rhs = 0
            for F in E:
                if apply_graph(F, g) == F:
                    e1 = eps_orient(F, g); e2 = eps_cycle(F, g)
                    if e1 != e2:
                        cyc_mism += 1
                    rhs += e1
            if fx != rhs:
                mism += 1
        # Burnside total
        tot = 0
        for part, g in cycle_type_reps(n).items():
            tot += class_size(part, n) * fixed_classes(n, g, reps)
        orbits = tot // math.factorial(n)
        good = (mism == 0 and cyc_mism == 0 and tot % math.factorial(n) == 0 and orbits == A049313[n - 1])
        ok &= good
        print(f"n={n}: perms checked={len(perms)} identity mismatches={mism} cycle-form mismatches={cyc_mism} "
              f"Burnside orbits={orbits} A049313={A049313[n - 1]} {'OK' if good else 'FAIL'}", flush=True)
        if n <= 6:
            # direct: orbits of untwisted Euler graphs
            allp = list(itertools.permutations(range(n)))
            canon_seen = set(); untw = 0; allorb = 0
            for F in E:
                c = min(tuple(sorted(apply_graph(F, g))) for g in allp)
                if c in canon_seen:
                    continue
                canon_seen.add(c); allorb += 1
                aut = [g for g in allp if apply_graph(F, g) == F]
                if all(eps_orient(F, g) == 1 for g in aut):
                    untw += 1
            good2 = (untw == A049313[n - 1] and allorb == A002854[n - 1])
            ok &= good2
            print(f"   direct: Euler graph orbits={allorb} (A002854={A002854[n - 1]}), untwisted={untw} {'OK' if good2 else 'FAIL'}", flush=True)
    # (4) cycle-form Burnside by F_2 linear algebra, n = 8, 9, 10
    for n in (8, 9, 10):
        tot = 0; totE = 0
        for part, g in cycle_type_reps(n).items():
            # edge orbits of g
            P = pairs(n); idx = {p: k for k, p in enumerate(P)}; seen = [False] * len(P); orbs = []
            for k, (i, j) in enumerate(P):
                if seen[k]:
                    continue
                orb = []; a, b = i, j
                while True:
                    kk = idx[(min(a, b), max(a, b))]
                    if seen[kk]:
                        break
                    seen[kk] = True; orb.append((min(a, b), max(a, b))); a, b = g[a], g[b]
                orbs.append(orb)
            # constraint rows: vertex degree parity, as bitmasks over orbits
            rows = []
            for v in range(n):
                r = 0
                for o, orb in enumerate(orbs):
                    if sum(1 for (a, b) in orb if v in (a, b)) % 2:
                        r |= 1 << o
                rows.append(r)
            # antipodal functional
            ell = 0
            ginv = g
            for o, orb in enumerate(orbs):
                a, b = orb[0]
                # orbit is antipodal iff g^{|orb|} swaps a and b
                x, y = a, b
                for _ in range(len(orb)):
                    x, y = g[x], g[y]
                if (x, y) == (b, a):
                    ell |= 1 << o
            def rank(vs):
                basis = []
                for v in vs:
                    for bvec in basis:
                        v = min(v, v ^ bvec)
                    if v:
                        basis.append(v)
                return len(basis)
            r0 = rank(rows); r1 = rank(rows + [ell])
            dimW = len(orbs) - r0
            cs = class_size(part, n)
            totE += cs * 2 ** dimW
            if r1 == r0:
                tot += cs * 2 ** dimW
        f = math.factorial(n)
        good = (tot % f == 0 and tot // f == A049313[n - 1] and totE // f == A002854[n - 1])
        ok &= good
        print(f"n={n}: cycle-form Burnside untwisted orbits={tot // f} (A049313={A049313[n - 1]}), Euler orbits={totE // f} (A002854={A002854[n - 1]}) {'OK' if good else 'FAIL'}", flush=True)
    return ok



# ------------------------------------------------------------------ Part I driver
import os, re, subprocess, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))

def run(binary, args, stdin_cmd=None):
    if stdin_cmd:
        p1 = subprocess.Popen(stdin_cmd, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
        r = subprocess.run([binary] + args, stdin=p1.stdout, capture_output=True, text=True)
        p1.wait()
    else:
        r = subprocess.run([binary] + args, capture_output=True, text=True)
    return r.stdout

def expect(label, text, pattern, ok_list):
    good = re.search(pattern, text) is not None
    ok_list.append(good)
    print(f"[{'OK' if good else 'FAIL'}] {label}", flush=True)

def c_main():
    oks = []
    with tempfile.TemporaryDirectory() as d:
        b = os.path.join(d, 'sa')
        subprocess.run(['cc', '-O2', '-o', b, os.path.join(HERE, 'procgen_selfie_20261001_orchestrator_check.c')], check=True)
        lab = {3: (2, 0, 2, 0, 2), 4: (40, 0, 0, 3, 5), 5: (664, 0, 184, 0, 8), 6: (26048, 240, 0, 5, 15), 7: (1934528, 0, 96608, 0, 18)}
        for n, (cov, aodd, aeven, mn, mx) in lab.items():
            out = run(b, ['labeled', str(n)]); print(out.strip())
            expect(f'labeled N={n} counts', out, rf'cover={cov} allodd={aodd} alleven={aeven} minodd={mn} maxodd={mx} bad\(cover&on-all\)=0 zformula_fail=0 oddcount_parity_fail=0 x0_fail=0 selfieRedei_fail=0', oks)
            expect(f'labeled N={n} switching sums', out, r"2N!: 0 ; sum_L \(-1\)\^\|L\| H != 0: 0", oks)
            expect(f'labeled N={n} Walsh degree bound', out, r'walsh_degree_fail=0', oks)
        expect('constant-H switching classes only at N=4 (two, H=3)', ''.join(run(b, ['labeled', '4'])), r'constant_H_classes=2', oks)
        cls = {3: (2, 1, 0, 1, 0, 2, 0), 4: (4, 3, 0, 0, 3, 5, 0), 5: (12, 7, 0, 3, 0, 8, 1), 6: (56, 44, 1, 0, 5, 15, 2),
               7: (456, 412, 0, 28, 0, 18, 9), 8: (6880, 6674, 0, 0, 7, 27, 33), 9: (191536, 189992, 0, 899, 0, 32, 167)}
        for n, (nc, cov, aodd, aeven, mn, mx, sd) in cls.items():
            args = ['classes', 'beta', 'shave'] if n <= 7 else ['classes']
            out = run(b, args, ['gentourng', '-q', str(n)]); print(out.strip())
            expect(f'classes N={n} census', out, rf'n={nc} cover={cov} allodd={aodd} alleven={aeven} minodd={mn} maxodd={mx} strong_with_dead_arc={sd} zformula_fail=0', oks)
            if n <= 7:
                expect(f'classes N={n} beta = Hall bound', out, r'beta!=hall=0', oks)
                ns = 1 if n == 6 else 0
                expect(f'classes N={n} Redei shaving to one HP (H odd at every step) fails for exactly {ns} class(es)', out, rf'no_redei_shaving={ns}\b', oks)
                if n == 6:
                    expect('N=6: the class with no Redei shaving is the all-odd one (H=45)', out, r'NO_REDEI_SHAVING 110101101110111 H=45', oks)
                if n % 2:
                    regs = re.findall(r'REGULAR \S+ beta=(\d+)', out)
                    expect(f'classes N={n}: every regular class has beta = N-1 ({len(regs)} classes)', out if regs and all(int(x) == n - 1 for x in regs) else '', r'.', oks)
        out = run(b, ['beta', '111000101111111111101']); print(out.strip())
        expect('N=7 sigma counterexample: beta = hall = 3 < sigma = 4', out, r'beta=3 hall=3 sigma=4', oks)
        for q, m, pat in [(7, 0, r'#odd arcs=0 of 21'), (7, 1, r'H=45 #odd arcs=15 of 15'), (11, 0, r'#odd arcs=0 of 55'),
                          (11, 1, r'H=15745 #odd arcs=45 of 45'), (19, 0, r'#odd arcs = 0 of 171'), (19, 1, r'#odd arcs=153 of 153'),
                          (23, 1, r'#odd arcs = 231 of 231')]:
            out = run(b, ['paley', str(q), str(m)]); print(out.strip().splitlines()[0])
            expect(f'QR_{q}{" minus 0" if m else ""}', out, pat, oks)
        for n in (3, 5, 7, 9, 11, 13):
            out = run(b, ['circ', str(n)]); print(out.strip().splitlines()[0])
            expect(f'circulants Z_{n} all-even', out, rf'{2 ** ((n - 1) // 2)} tournaments, all-even: {2 ** ((n - 1) // 2)}', oks)
        out = run(b, ['z3z3']); print(out.strip().splitlines()[0])
        expect('Cayley Z3xZ3 all-even', out, r'16 tournaments, all-even: 16', oks)
        for t, pat in [('111111001110111111111101110111111111110111110', r'H=3929 #odd arcs=45 of 45'),
                       ('101010110111010101100101111100110111101110111', r'H=15745 #odd arcs=45 of 45'),
                       ('110011110111110011110001111111111001111111110', r'#odd arcs=5 of 45'),
                       ('111001110101110111111111111101101111101111101', r'#odd arcs=5 of 45'),
                       ('111000101101101101111011111101110111111111101', r'#odd arcs=5 of 45'),
                       ('101011111111011011111011110111110001111111101', r'#odd arcs=5 of 45')]:
            out = run(b, ['one', t]); print(out.strip().splitlines()[0][:90])
            expect(f'N=10 class {t[:12]}...', out, pat, oks)
    return all(oks)

if __name__ == "__main__":
    import time
    t0 = time.time()
    print("==== Part I: C engine ====", flush=True)
    ok1 = c_main()
    print("==== Part II: Theorem E1 (OPEN-Q-060) ====", flush=True)
    ok2 = euler_main()
    print(f"elapsed {time.time() - t0:.1f} s")
    print("ALL CHECKS PASSED" if (ok1 and ok2) else "SOME CHECK FAILED")
