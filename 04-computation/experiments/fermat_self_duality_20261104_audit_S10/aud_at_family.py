"""aud_at_family.py -- independent audit checks (fermat_self_duality_20261104_audit_S10) of the PROVED AT_n claims in
fermat_self_duality_fewcycle_coincidences_20261104.md.  Written from scratch; no code from the session's scripts.

Vertices are 0..n-1 standing for the labels 1..n.  AT_n: every arc i -> j (i < j) except the extreme pair, n-1 -> 0.
Hamiltonian paths are enumerated two independent ways: all n! permutations (n <= 9) and an explicit-stack bitmask DFS
(n <= 16).  Checks:
 (a) substitution form C3[1, TT_(n-2), 1];  (b) every induced subtournament is AT_k or TT_k (all subsets, n <= 12), and
 AT_k, TT_k are self-converse;  (c) H(AT_n) = 2^(n-2) + 1;  (d) the mirror law;  (e) parity H = Fix (mod 2) under an
 involutive anti-automorphism (random self-converse tournaments built from an involution);  (f) the x^(1/x) pair
 (x, x^2);  (g) end_(m-1) + end_m;  (h) W4 = AT_4 and S75's block list (brute force, orders 3..6) plus the block
 conditions for (AT_m, {m-1, m});  (i) the generating-function identity;  (j) |Aut(AT_n)|;  (k) C3 modules at
 orders 7 and 8."""
import itertools, random, sys
from math import comb
import networkx as nx
from networkx.algorithms.isomorphism import DiGraphMatcher
import sympy as sp

random.seed(20261104)
FAIL = []

def check(cond, msg):
    if not cond:
        FAIL.append(msg)
        print("  FAIL:", msg)
    return cond

def at(n):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if (i, j) == (0, n - 1):
                out[j] |= 1 << i
            else:
                out[i] |= 1 << j
    return out

def tt(n):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            out[i] |= 1 << j
    return out

def arc(out, u, v):
    return (out[u] >> v) & 1

def paths_perm(n, out):
    res = []
    for p in itertools.permutations(range(n)):
        if all(arc(out, p[i], p[i + 1]) for i in range(n - 1)):
            res.append(p)
    return res

def paths_dfs(n, out):
    res = []
    full = (1 << n) - 1
    for s in range(n):
        stack = [(s, 1 << s, (s,))]
        while stack:
            v, used, path = stack.pop()
            if used == full:
                res.append(path)
                continue
            nb = out[v] & ~used & full
            while nb:
                w = (nb & -nb).bit_length() - 1
                nb &= nb - 1
                stack.append((w, used | (1 << w), path + (w,)))
    return res

def H_dp(n, out):
    full = (1 << n) - 1
    dp = [[0] * n for _ in range(full + 1)]
    for v in range(n):
        dp[1 << v][v] = 1
    for S in range(1, full + 1):
        for v in range(n):
            c = dp[S][v]
            if c:
                nb = out[v] & ~S
                while nb:
                    w = (nb & -nb).bit_length() - 1
                    nb &= nb - 1
                    dp[S | (1 << w)][w] += c
    return sum(dp[full]), dp[full]

def to_nx(n, out):
    G = nx.DiGraph()
    G.add_nodes_from(range(n))
    for u in range(n):
        for v in range(n):
            if arc(out, u, v):
                G.add_edge(u, v)
    return G

def converse(n, out):
    op = [0] * n
    for u in range(n):
        for v in range(n):
            if arc(out, u, v):
                op[v] |= 1 << u
    return op

def induced(out, S):
    S = sorted(S)
    k = len(S)
    sub = [0] * k
    for a in range(k):
        for b in range(k):
            if a != b and arc(out, S[a], S[b]):
                sub[a] |= 1 << b
    return sub

def strong(n, out):
    G = to_nx(n, out)
    return nx.is_strongly_connected(G)

print("== (a) substitution form AT_n = C3[1, TT_(n-2), 1], n = 3..20")
for n in range(3, 21):
    out = at(n)
    I = list(range(1, n - 1))
    mod = all(len({arc(out, x, m) for m in I}) == 1 for x in (0, n - 1))
    trans = induced(out, I) == tt(len(I))
    quot = all(arc(out, 0, m) for m in I) and all(arc(out, m, n - 1) for m in I) and arc(out, n - 1, 0)
    check(mod and trans and quot, f"substitution form fails at n={n}")
print("  I = {2..n-1} is a transitive module and the quotient on {1, I, n} is 1 -> I -> n -> 1: verified for n = 3..20")

print("== (b) induced subtournaments of AT_n are AT_k (if they contain 1 and n) or TT_k, all subsets, n = 3..12")
for n in range(3, 13):
    out = at(n)
    for r in range(2, n + 1):
        for S in itertools.combinations(range(n), r):
            sub = induced(out, S)
            want = at(r) if (0 in S and n - 1 in S and r >= 2) else tt(r)
            if r == 2 and 0 in S and n - 1 in S:
                want = [0, 1]          # the single arc n -> 1
            check(sub == want, f"heredity fails n={n} S={S}")
selfconv = []
for k in range(1, 17):
    for name, out in (("AT", at(k) if k >= 3 else tt(k)), ("TT", tt(k))):
        sig = lambda v: k - 1 - v
        ok = all(arc(out, u, v) == arc(out, sig(v), sig(u)) for u in range(k) for v in range(k) if u != v)
        selfconv.append(ok)
check(all(selfconv), "AT_k or TT_k not self-converse under i -> k+1-i")
print("  verified for every subset at n = 3..12; AT_k and TT_k are self-converse via i -> k+1-i (k <= 16)")

print("== (c), (d) path count and the mirror law (two independent path enumerators)")
for n in range(3, 17):
    out = at(n)
    m = n - 2
    P = paths_dfs(n, out)
    if n <= 9:
        Pp = paths_perm(n, out)
        check(sorted(Pp) == sorted(P), f"perm and DFS path sets differ at n={n}")
    if n <= 14:
        h, _ = H_dp(n, out)
        check(h == len(P), f"DP count differs at n={n}")
    check(len(P) == 2 ** m + 1 and len(set(P)) == len(P), f"H(AT_{n}) != 2^(n-2)+1")
    sig = lambda v: n - 1 - v
    check(all(arc(out, u, v) == arc(out, sig(v), sig(u)) for u in range(n) for v in range(n) if u != v),
          f"sigma not an anti-automorphism at n={n}")
    S = set(P)
    iota = {p: tuple(sig(v) for v in reversed(p)) for p in P}
    check(all(q in S for q in iota.values()), f"iota leaves the path set at n={n}")
    check(all(iota[iota[p]] == p for p in P), f"iota not an involution at n={n}")
    fix = sum(1 for p in P if iota[p] == p)
    exp = 1 + 2 ** (m // 2) if m % 2 == 0 else 1
    check(fix == exp, f"Fix mismatch at n={n}: {fix} vs {exp}")
    # explicit form A n 1 B, binomial row and the mirror j -> m - j
    dist = [0] * (m + 1)
    ident = tuple(range(n))
    check(ident in S, f"transitive path missing at n={n}")
    check(iota[ident] == ident, f"transitive path not fixed at n={n}")
    for p in P:
        if p == ident:
            continue
        t = p.index(n - 1)
        check(t + 1 < n and p[t + 1] == 0, f"non-transitive path without n -> 1 at n={n}: {p}")
        A, B = p[:t], p[t + 2:]
        check(list(A) == sorted(A) and list(B) == sorted(B) and set(A) | set(B) == set(range(1, n - 1)),
              f"path not of the form A n 1 B at n={n}")
        dist[len(A)] += 1
        q = iota[p]
        sB = tuple(sorted(sig(v) for v in B)); sA = tuple(sorted(sig(v) for v in A))
        check(q == sB + (n - 1, 0) + sA, f"iota(A n 1 B) != sigma(B) n 1 sigma(A) at n={n}")
        check(len(q[:q.index(n - 1)]) == m - len(A), f"mirror j -> m-j fails at n={n}")
    check(dist == [comb(m, j) for j in range(m + 1)], f"binomial row fails at n={n}")
    line = f"  n={n:2d}: H={len(P):6d} = 2^{m}+1, Fix(iota)={fix:4d}"
    if m % 2 == 0:
        k = m // 2
        Hk = 2 ** k + 1
        check(fix == Hk, f"Fix(AT_{n}) != H(AT_{k + 2})")
        line += f" = H(AT_{k + 2}) = {Hk}"
    if m >= 1 and (m & (m - 1)) == 0:
        i = m.bit_length() - 1
        line += f"   [n = 2^{i} + 2: H = F_{i} = {2 ** m + 1}, Fix = " + (f"F_{i - 1} = {fix}]" if i >= 1 else f"{fix}]")
    print(line)

print("== (e) parity: H = Fix(iota) (mod 2) for tournaments with an involutive anti-automorphism (random, n = 4..9)")
def random_with_involutive_anti(n):
    verts = list(range(n))
    random.shuffle(verts)
    sig = list(range(n))
    for a in range(0, n - 1, 2):
        x, y = verts[a], verts[a + 1]
        sig[x], sig[y] = y, x            # at most one fixed vertex (n odd)
    out = [0] * n
    done = set()
    for u in range(n):
        for v in range(u + 1, n):
            if (u, v) in done:
                continue
            # orbit of the unordered pair {u, v} under sigma; arc u -> v forces sig(v) -> sig(u)
            if random.random() < 0.5:
                a, b = u, v
            else:
                a, b = v, u
            a2, b2 = sig[b], sig[a]
            if {a2, b2} == {a, b} and (a2, b2) != (a, b):
                # pair fixed by sigma with reversal: sigma swaps u and v, consistent (arc a->b maps to b->a? no)
                pass
            out[a] |= 1 << b
            out[a2] |= 1 << b2
            done.add((min(a, b), max(a, b))); done.add((min(a2, b2), max(a2, b2)))
    return out, sig

tested = 0
for n in range(4, 10):
    for _ in range(40):
        out, sig = random_with_involutive_anti(n)
        # validity: tournament and sigma anti-automorphism
        okT = all(arc(out, u, v) + arc(out, v, u) == 1 for u in range(n) for v in range(u + 1, n))
        okS = okT and all(arc(out, u, v) == arc(out, sig[v], sig[u]) for u in range(n) for v in range(n) if u != v)
        if not okS:
            continue
        P = paths_dfs(n, out)
        fix = sum(1 for p in P if tuple(sig[v] for v in reversed(p)) == p)
        check(len(P) % 2 == fix % 2 and fix % 2 == 1, f"parity reading fails n={n}")
        tested += 1
print(f"  {tested} random tournaments with an involutive anti-automorphism: H = Fix (mod 2) and Fix odd in all")

print("== (f) (x, x^2) as an x^(1/x) mirror pair")
x = sp.symbols('x', positive=True)
sols = sp.solve(sp.Eq(x ** 2 * sp.log(x), 2 * x * sp.log(x)), x)
print(f"  x^(x^2) = (x^2)^x, x > 0  <=>  x^2 ln x = 2x ln x: solutions {sorted(sols)}  (x = 1 gives the degenerate pair (1, 1))")

print("== (g) end counts of AT_m: end_(m-1) + end_m = 2^(m-3) + 1 = H(AT_(m-1))")
for m in range(4, 17):
    out = at(m)
    ends = [0] * m
    for p in paths_dfs(m, out):
        ends[p[-1]] += 1
    ok = ends[m - 2] == 2 ** (m - 3) and ends[m - 1] == 1
    check(ok and ends[m - 2] + ends[m - 1] == 2 ** (m - 3) + 1, f"end counts fail at m={m}: {ends}")
print("  end_(m-1) = 2^(m-3) and end_m = 1 for m = 4..16")

print("== (h) W4 = AT_4; S75's blocks by brute force (orders 3..6); (AT_m, {m-1, m}) as a block")
w4 = [0] * 4
for (u, v) in [(0, 1), (1, 2), (2, 3), (3, 0), (0, 2), (1, 3)]:     # a->b->c->d->a, a->c, b->d
    w4[u] |= 1 << v
check(w4 == at(4), "W4 (a,b,c,d) = (1,2,3,4) is not AT_4")
strong4 = []
for code in range(64):
    out = [0] * 4
    bit = 0
    for j in range(1, 4):
        for i in range(j):
            if (code >> bit) & 1:
                out[i] |= 1 << j
            else:
                out[j] |= 1 << i
            bit += 1
    if strong(4, out):
        strong4.append(out)
check(all(nx.is_isomorphic(to_nx(4, o), to_nx(4, at(4))) for o in strong4), "a strong 4-tournament is not AT_4")
_, endsW4 = H_dp(4, at(4))
print(f"  W4 = AT_4 exactly; {len(strong4)} labelled strong 4-tournaments, all isomorphic to AT_4; "
      f"end counts (a,b,c,d) = {tuple(endsW4)}; sole exits = {[v + 1 for v in range(4) if bin(at(4)[v]).count('1') == 1]}")

def block_ok(n, out, A, ends):
    k = sum(ends[a] for a in range(n) if (A >> a) & 1)
    if any(bin(out[v]).count("1") < 2 for v in range(n) if not (A >> v) & 1):
        return None
    full = (1 << n) - 1
    for Z in range(1, full + 1):
        dz = sum(bin(out[z] & ~Z & full).count("1") for z in range(n) if (Z >> z) & 1)
        if dz + bin(Z & A).count("1") == 1:
            return None
    return k

def canon_pair(n, out, A):
    best = None
    for p in itertools.permutations(range(n)):
        code = tuple(arc(out, p.index(i), p.index(j)) for j in range(n) for i in range(j))
        am = tuple(1 if (A >> p.index(i)) & 1 else 0 for i in range(n))
        key = (code, am)
        if best is None or key < best:
            best = key
    return best

blocks = {}
for n in range(3, 7):
    E = n * (n - 1) // 2
    seen = set()
    for code in range(1 << E):
        out = [0] * n
        bit = 0
        for j in range(1, n):
            for i in range(j):
                if (code >> bit) & 1:
                    out[i] |= 1 << j
                else:
                    out[j] |= 1 << i
                bit += 1
        if not strong(n, out):
            continue
        _, ends = H_dp(n, out)
        for A in range(1, 1 << n):
            if bin(A).count("1") < 2:
                continue
            if sum(ends[a] for a in range(n) if (A >> a) & 1) > 5:
                continue
            k = block_ok(n, out, A, ends)
            if k is None:
                continue
            key = canon_pair(n, out, A)
            if key in seen:
                continue
            seen.add(key)
            hc_uh = None
            blocks.setdefault(k, []).append((n, out, A))
for k in sorted(blocks):
    desc = []
    for (n, out, A) in blocks[k]:
        iso_at = nx.is_isomorphic(to_nx(n, out), to_nx(n, at(n))) if n >= 3 else False
        tag = f"order {n}" + (" (AT_%d)" % n if iso_at else "") + f", |A| = {bin(A).count('1')}"
        if iso_at:
            # is A the set of the two sole exits?
            se = [v for v in range(n) if bin(out[v]).count("1") == 1]
            tag += " A = all" if A == (1 << n) - 1 else (" A = the two sole exits" if sorted(se) == [v for v in range(n) if (A >> v) & 1] else "")
        desc.append(tag)
    print(f"  blocks of type {k} (orders 3..6): {len(blocks[k])}: " + "; ".join(desc))
check(len(blocks.get(2, [])) == 0, "a type-2 block exists")
check(len(blocks.get(3, [])) == 2, "type-3 block count != 2")
# (AT_5, {4, 5}) among the type-5 blocks
at5_found = any(n == 5 and nx.is_isomorphic(to_nx(5, out), to_nx(5, at(5))) and bin(A).count("1") == 2 for (n, out, A) in blocks.get(5, []))
check(at5_found, "(AT_5, sole exits) not among the type-5 blocks")
print(f"  (AT_5, its two sole exits) is one of the order-5 type-5 blocks: {at5_found}")
for m in range(4, 13):
    out = at(m)
    _, ends = H_dp(m, out)
    A = (1 << (m - 2)) | (1 << (m - 1))
    k = block_ok(m, out, A, ends)
    check(k == 2 ** (m - 3) + 1, f"(AT_{m}, {{m-1, m}}) fails the block conditions")
print("  (AT_m, {m-1, m}) satisfies all of S75's block conditions with type 2^(m-3) + 1, for every m = 4..12")

print("== (i) u_3(n) = F(2n - 9) = R(n - 2) and U_3(x) = x^2 (R(x) - x^2 - x^3)")
X = sp.symbols('X')
R = X ** 2 * (1 - 2 * X) / (1 - 3 * X + X ** 2)
U3 = X ** 6 * (2 - X) / (1 - 3 * X + X ** 2)
check(sp.simplify(X ** 2 * (R - X ** 2 - X ** 3) - U3) == 0, "GF identity fails")
sR = sp.series(R, X, 0, 45).removeO(); sU = sp.series(U3, X, 0, 45).removeO()
F = lambda k: sp.fibonacci(k) if k >= 0 else (-1) ** (k + 1) * sp.fibonacci(-k)
okR = all(sR.coeff(X, mm) == F(2 * mm - 5) for mm in range(2, 43))
okU = all(sU.coeff(X, nn) == (F(2 * nn - 9) if nn >= 6 else 0) for nn in range(0, 45))
okRel = all(sU.coeff(X, nn) == sR.coeff(X, nn - 2) for nn in range(6, 45))
check(okR and okU and okRel, "coefficient identities fail")
print(f"  symbolic identity holds; R(m) = F(2m-5) for m = 2..42 (R(2) = F(-1) = 1, R(3) = 1); u_3(n) = F(2n-9) = R(n-2) for n = 6..44; "
      f"u_3(4) = u_3(5) = 0 (the subtracted x^4, x^5)")

print("== (j) automorphism groups of AT_n")
auts = []
for n in range(3, 13):
    G = to_nx(n, at(n))
    a = sum(1 for _ in DiGraphMatcher(G, G).isomorphisms_iter())
    auts.append(a)
print(f"  |Aut(AT_n)|, n = 3..12: {auts}")
check(auts[0] == 3 and all(a == 1 for a in auts[1:]), "AT_n not rigid for some n >= 4")

print("== (k) C3 modules in {-3}-self-dual tournaments")
def tour(n, arcs):
    out = [0] * n
    for (u, v) in arcs:
        out[u] |= 1 << v
    return out
c3 = [(0, 1), (1, 2), (2, 0)]
x_dom = tour(4, c3 + [(3, 0), (3, 1), (3, 2)])
c3_dom = tour(4, c3 + [(0, 3), (1, 3), (2, 3)])
iso = lambda o1, o2, n: nx.is_isomorphic(to_nx(n, o1), to_nx(n, o2))
check(iso(converse(4, x_dom), c3_dom, 4) and not iso(x_dom, c3_dom, 4), "order-7 lemma fails")
print("  x => C3 and C3 => x are each other's converse and not isomorphic, so neither is self-converse (order-7 lemma)")
# order 8: M = C3 module plus two outside vertices x, y
res = []
for xs in (0, 1):            # x dominates M (1) or is dominated (0)
    for ys in (0, 1):
        for xy in (0, 1):
            arcs = list(c3)
            for (w, s) in ((3, xs), (4, ys)):
                arcs += [(w, m) for m in range(3)] if s else [(m, w) for m in range(3)]
            arcs.append((3, 4) if xy else (4, 3))
            o = tour(5, arcs)
            res.append((xs, ys, xy, iso(o, converse(5, o), 5)))
same_side_sc = [r for r in res if r[0] == r[1] and r[3]]
opp_side_sc = [r for r in res if r[0] != r[1] and r[3]]
print(f"  M (C3 module) + x + y: self-converse configurations: {[(r[0], r[1], r[2]) for r in res if r[3]]} "
      f"(1 = dominates M); same-side configurations self-converse: {len(same_side_sc)}")
check(len(same_side_sc) == 0 and len(opp_side_sc) == 4, "order-8 configuration table unexpected")
print("  => at order 8 a C3 module leaves 5 outside vertices, two of them on the same side; deleting the other three\n"
      "     gives a non-self-converse 5-tournament.  So C3 modules are excluded at order 8 by proof as well.")

print("== (l) section 1 Fact 4: AT_n uniquely Hamiltonian, non-critical run = interior, 3 cut arcs, sole exits {n-1, n} (n = 4..12)")
for n in range(4, 13):
    out = at(n)
    G0 = to_nx(n, out)
    cyc = [c for c in nx.simple_cycles(G0) if len(c) == n]
    nc = [v for v in range(n) if nx.is_strongly_connected(to_nx(n - 1, induced(out, [u for u in range(n) if u != v])))]
    C = cyc[0] if cyc else []
    cyc_arcs = [(C[i], C[(i + 1) % n]) for i in range(n)] if C else []
    cut = [e for e in cyc_arcs if not nx.is_strongly_connected(nx.restricted_view(G0, [], [e]))]
    se = [v for v in range(n) if bin(out[v]).count("1") == 1]
    sn = [v for v in range(n) if n - 1 - bin(out[v]).count("1") == 1]
    ok = len(cyc) == 1 and nc == list(range(1, n - 1)) and len(cut) == 3 and se == [n - 2, n - 1] and sn == [0, 1]
    check(ok, f"Fact 4 fails at n={n}: hc={len(cyc)} nc={nc} cut={cut} sole exits={se}")
print("  hc = 1, non-critical vertices = {2..n-1} (run of length n-2), 3 cut arcs, sole exits {n-1, n}, sole entries {1, 2}: n = 4..12"
      " (n = 3 is C3: three sole exits, no non-critical vertex; THM-4868 assumes n >= 4)")

print("== (m) the arithmetic of section 4")
ogg = [2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 41, 47, 59, 71]
fermat = [3, 5, 17, 257, 65537]
print(f"  Fermat primes & Ogg primes = {sorted(set(ogg) & set(fermat))}; primes <= 31 all supersingular: "
      f"{all(p in ogg for p in sp.primerange(2, 32))}")
fi22, B_ = {2, 3, 5, 7, 11, 13}, {2, 3, 5, 7, 11, 13, 17, 19, 23, 31, 47}
print(f"  3, 5 in the Fi22 layer: {3 in fi22 and 5 in fi22}; 17 in B but not Fi22: {17 in B_ and 17 not in fi22}")
plat = [(p, q) for p in range(3, 10) for q in range(3, 10) if sp.Rational(1, p) + sp.Rational(1, q) > sp.Rational(1, 2)]
print(f"  Schlafli pairs with 1/p + 1/q > 1/2: {plat}; face sizes {sorted({p for p, q in plat})}")
print(f"  2^m + 1 prime for m <= 64 only at m = {[m for m in range(1, 65) if sp.isprime(2 ** m + 1)]} (powers of 2)")

print()
print("ALL CHECKS PASSED" if not FAIL else f"{len(FAIL)} FAILURES")
