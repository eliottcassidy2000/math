# Audit S83 / THM-4886, step 1: identities, written independently of the session's scripts.
# (i)   H(w) from an explicit tournament (Hamiltonian-path DP over subsets) equals g(A(w)) for ALL words of length <= 9.
# (ii)  Concatenation laws (1) and (2) of Session 193, EXHAUSTIVELY for |u|,|v|,|w| <= 4 (|v| >= 1), plus random long words.
# (iii) Split facts used in the proof: with p(u), eps(w) fixed, A(u) = 2^(a-1) + p + 2j (0 <= j < 2^(a-2)),
#       C(w) = eps + 2i (0 <= i < 2^(l-1)); A and C are bijections {0,1}^L -> [0, 2^L); number of ordered pairs
#       (v, v') with C(v) - C(v') = d (and with A(v) - A(v') = d) is 2^L - |d|.
# (iv)  The transfer x = x' 1^t used in (C), H odd, and g(v) <= v + 1.
import itertools, random, sys

def ham_paths(nv, adj):
    # adj[i] = bitmask of out-neighbours of i; count directed Hamiltonian paths
    full = (1 << nv) - 1
    dp = [[0] * nv for _ in range(1 << nv)]
    for i in range(nv):
        dp[1 << i][i] = 1
    for S in range(1, 1 << nv):
        row = dp[S]
        for i in range(nv):
            c = row[i]
            if not c:
                continue
            nb = adj[i] & ~S
            while nb:
                b = nb & -nb
                j = b.bit_length() - 1
                dp[S | b][j] += c
                nb ^= b
    return sum(dp[full])

def H_tournament(w, chain_up=True):
    # vertices 0..m-1 chain, vertex m apex; chain i -> j iff i < j (chain_up) ; apex -> i iff w_i = 1
    m = len(w)
    adj = [0] * (m + 1)
    for i in range(m):
        for j in range(m):
            if i != j and ((i < j) == chain_up):
                adj[i] |= 1 << j
    for i in range(m):
        if w[i] == 1:
            adj[m] |= 1 << i
        else:
            adj[i] |= 1 << m
    return ham_paths(m + 1, adj)

def g(v):
    if v == 0:
        return 1
    return g(v >> 1) if v & 1 else g(v >> 1) + (v >> 1) + ((v >> 1) & 1)

def A(w):
    a = 0
    for b in w:
        a = 2 * a + b
    return a

def C(w):
    return sum((1 - b) << j for j, b in enumerate(w))

def p(w):
    return w[-1] if w else 0

def eps(w):
    return 1 if (w and w[0] == 0) else 0

def H(w):
    return g(A(w))

out = []
# (i)
for chain_up in (True, False):
    bad = 0; tot = 0
    for m in range(0, 10):
        for w in itertools.product((0, 1), repeat=m):
            tot += 1
            if H_tournament(list(w), chain_up) != H(list(w)):
                bad += 1
    out.append(f"(i) chain orientation i->j iff i<j is {chain_up}: mismatches H_tournament vs g(A) over all {tot} words of length 0..9: {bad}")

# (ii) exhaustive
bad1 = bad2 = cnt = 0
words = [list(w) for m in range(0, 5) for w in itertools.product((0, 1), repeat=m)]
for u in words:
    for v in words:
        if not v:
            continue
        if H(u + v) != H(u) + H(v) - 1 + A(u) * C(v) + p(u) * eps(v):
            bad1 += 1
        for w in words:
            cnt += 1
            L = len(v)
            rhs = H(u) + H(v) + H(w) - 2 + A(u) * C(v) + ((1 << L) * A(u) + A(v)) * C(w) + p(u) * eps(v) + p(v) * eps(w)
            if H(u + v + w) != rhs:
                bad2 += 1
out.append(f"(ii) exhaustive |u|,|w| in 0..4, |v| in 1..4: law (1) failures {bad1}; law (2) failures {bad2} of {cnt} triples")
random.seed(20261106)
bad = 0
for t in range(3000):
    u = [random.randint(0, 1) for _ in range(random.randint(0, 40))]
    v = [random.randint(0, 1) for _ in range(random.randint(1, 40))]
    w = [random.randint(0, 1) for _ in range(random.randint(0, 40))]
    L = len(v)
    rhs = H(u) + H(v) + H(w) - 2 + A(u) * C(v) + ((1 << L) * A(u) + A(v)) * C(w) + p(u) * eps(v) + p(v) * eps(w)
    bad += H(u + v + w) != rhs
out.append(f"(ii) random long triples (lengths up to 40): law (2) failures {bad} of 3000")

# (iii)
okA = okC = okbij = okcount = True
for a in range(2, 9):
    for pp in (0, 1):
        vals = sorted(A([1] + list(mid) + [pp]) for mid in itertools.product((0, 1), repeat=a - 2))
        okA &= vals == [(1 << (a - 1)) + pp + 2 * j for j in range(1 << (a - 2))]
for l in range(1, 9):
    for e in (0, 1):
        vals = sorted(C([1 - e] + list(rest)) for rest in itertools.product((0, 1), repeat=l - 1))
        okC &= vals == [e + 2 * i for i in range(1 << (l - 1))]
for L in range(1, 9):
    V = [list(v) for v in itertools.product((0, 1), repeat=L)]
    Avals = sorted(A(v) for v in V); Cvals = sorted(C(v) for v in V)
    okbij &= Avals == list(range(1 << L)) and Cvals == list(range(1 << L))
    from collections import Counter
    cC = Counter(C(v) - C(w) for v in V for w in V); cA = Counter(A(v) - A(w) for v in V for w in V)
    okcount &= all(cC[d] == (1 << L) - abs(d) and cA[d] == (1 << L) - abs(d) for d in range(-(1 << L) + 1, 1 << L))
    # relation C(v) = 2^L - 1 - rev(A(v))  (reverse complement)
    okbij &= all(C(v) == (1 << L) - 1 - A(v[::-1]) for v in V)
out.append(f"(iii) A(u) progression (diff 2, length 2^(a-2)): {okA}; C(w) progression (diff 2, length 2^(l-1)): {okC}; "
           f"A, C bijections onto [0,2^L) with C(v) = 2^L-1-A(reverse v): {okbij}; pair counts 2^L-|d|: {okcount}")

# (iv)
ok_transfer = True; ok_odd = True; ok_le = True
for n in range(1, 15):
    for x in itertools.product((0, 1), repeat=n):
        x = list(x)
        if x[0] != 1:
            continue
        hx = H(x)
        ok_odd &= hx % 2 == 1
        if 0 in x:
            m = max(i for i in range(n) if x[i] == 0) + 1
            xp = x[:m]
            ok_transfer &= (xp[0] == 1 and xp[-1] == 0 and m >= 2 and H(xp) == hx)
        else:
            ok_transfer &= hx == 1
for v in range(0, 1 << 18):
    if g(v) > v + 1:
        ok_le = False; break
out.append(f"(iv) x = x'1^t with x' strong and H(x) = H(x') for all x in X_n, n <= 14: {ok_transfer}; H odd: {ok_odd}; g(v) <= v+1 for v < 2^18: {ok_le}")
print("\n".join(out))
