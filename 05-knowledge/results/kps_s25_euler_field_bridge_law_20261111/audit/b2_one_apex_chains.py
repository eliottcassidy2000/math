"""b2_one_apex_chains.py - audit of Corollary B2 (chain law, i >= 1) and of the law D1 of addendum D (direct joins,
i >= 0) by brute force, with the auditor's own toolkit tt.py (no import from scripts/).

  K1  one-apex data, all words of length <= 8:  H = g(P), W+ = P + 1, W- = Q + 1, tau(Y - apex) = 2^n,
      o = 1 + P//2, e = (P+1)//2, (odd/even split of W-) O = 1 + Q//2, E = (Q+1)//2, g(Q) = g(P),
      Q + 1 = 2^k (R(p) + 1) and n = k + L(p) for k leading letters 0 (R = reversal of the complement on L(p) bits)
  K2  Corollary B2: random chains Y_0 => C_1 => Y_1 => ... => C_k => Y_k, depth 2..4, words of length 0..4, chain
      lengths 1..3, at most 17 vertices: H(T) and W+_T(last apex) against the product of the matrices M(p, i)
      (the product is formed with the NORMAL FORM: leading letters 0 moved into the exponent)
  K3  Law D1: random chains, depth 2..3 (and 4), chain lengths 0..2: (A, o, e) against brute force at EVERY stage
      (H, and the parity split of W+ at the new apex)
  K4  Redundancy of the parametrisation (used to judge the 'tuple' count of C5): trailing letters 1 of a block equal
      extra chain vertices after it:  M(p_2, i_2) M(2 p_1 + 1, i_1) = M(p_2, i_2 + 1) M(p_1, i_1)
"""
import random
from itertools import product

from tt import tables, invariants, ordinal_sum, reverse_arc, transitive, one_apex, g_fun, word_PQ


def L_of(p):
    return p.bit_length()


def R_of(p):
    L = p.bit_length()
    comp = ((1 << L) - 1) ^ p
    return int(format(comp, "0%db" % L)[::-1], 2) if L else 0


def M(p, i):
    return ((g_fun(p), (1 << i) * (R_of(p) + 1)), (p + 1, 1 << (i + L_of(p))))


def matmul(A, B):
    return tuple(tuple(sum(A[r][k] * B[k][c] for k in range(2)) for c in range(2)) for r in range(2))


def build(words, chains):
    blocks = []
    apex_block = []
    for j, w in enumerate(words):
        if j > 0 and chains[j - 1] >= 1:
            blocks.append(transitive(chains[j - 1]))
        apex_block.append(len(blocks))
        blocks.append(one_apex(w))
    n, adj, offs = ordinal_sum(blocks)
    apexes = [offs[b] for b in apex_block]
    for j in range(1, len(apexes)):
        reverse_arc(adj, apexes[j - 1], apexes[j])
    return n, adj, apexes


def k1():
    tot = bad = 0
    for n in range(0, 9):
        for w in product((0, 1), repeat=n):
            m, adj = one_apex(w)
            inv = invariants(m, adj, 0)
            P, Q = word_PQ(w)
            k = 0
            while k < n and w[k] == 0:
                k += 1
            ok = (inv["H"] == g_fun(P) and inv["Wp"] == P + 1 and inv["Wm"] == Q + 1 and inv["tau_minus_root"] == 1 << n
                  and inv["o"] == 1 + P // 2 and inv["e"] == (P + 1) // 2 and inv["Om"] == 1 + Q // 2 and inv["Em"] == (Q + 1) // 2
                  and g_fun(Q) == g_fun(P) and Q + 1 == (1 << k) * (R_of(P) + 1) and n == k + L_of(P))
            tot += 1
            bad += (not ok)
    print("K1  one-apex data on all words of length <= 8: %d words, %d failures" % (tot, bad))


def k2():
    random.seed(5148)
    tot = bad = badW = 0
    while tot < 400:
        depth = random.randint(2, 4)
        words = [tuple(random.randint(0, 1) for _ in range(random.randint(0, 4))) for _ in range(depth)]
        chains = [random.randint(1, 3) for _ in range(depth - 1)]
        n = sum(len(w) + 1 for w in words) + sum(chains)
        if n > 17:
            continue
        n, adj, apexes = build(words, chains)
        inv = invariants(n, adj, apexes[-1])
        P0, _ = word_PQ(words[0])
        prod = M(P0, 1)
        for j in range(1, depth):
            w = words[j]
            k = 0
            while k < len(w) and w[k] == 0:
                k += 1
            Pj, _ = word_PQ(w)
            prod = matmul(M(Pj, chains[j - 1] + k), prod)
        tot += 1
        bad += (prod[0][0] != inv["H"])
        badW += (prod[1][0] != inv["Wp"])
    print("K2  Corollary B2 on %d random chains (depth 2..4, words <= 4 letters, chains 1..3, <= 17 vertices): H mismatches %d, W+ mismatches %d" % (tot, bad, badW))


def d1_step(state, w, i):
    A, o, e = state
    if i >= 1:
        s = (o + e) << (i - 1)
        o, e = s, s
    n = len(w)
    P, Q = word_PQ(w)
    O, E = 1 + Q // 2, (Q + 1) // 2
    te, to = (1 << (n - 1), 1 << (n - 1)) if n >= 1 else (1, 0)
    return (A * g_fun(P) + 2 * (o * E + e * O),
            A * (1 + P // 2) + o * to + e * te,
            A * ((P + 1) // 2) + o * te + e * to)


def k3():
    random.seed(9675)
    tot = bad = stages = 0
    while tot < 600:
        depth = random.randint(2, 4)
        words = [tuple(random.randint(0, 1) for _ in range(random.randint(0, 4))) for _ in range(depth)]
        chains = [random.randint(0, 2) for _ in range(depth - 1)]
        n = sum(len(w) + 1 for w in words) + sum(chains)
        if n > 16:
            continue
        P0, _ = word_PQ(words[0])
        state = (g_fun(P0), 1 + P0 // 2, (P0 + 1) // 2)
        ok = True
        for j in range(1, depth):
            state = d1_step(state, words[j], chains[j - 1])
            m, adj, apexes = build(words[:j + 1], chains[:j])
            inv = invariants(m, adj, apexes[-1])
            stages += 1
            if state != (inv["H"], inv["o"], inv["e"]):
                ok = False
        tot += 1
        bad += (not ok)
    print("K3  law D1 on %d random chains (depth 2..4, words <= 4 letters, chain lengths 0..2, <= 16 vertices), %d stages, state (A, o, e) compared at every stage: %d failures" % (tot, stages, bad))
    # exhaustive small: depth 2 and 3, words of length <= 2, chain lengths 0, 1, 2
    tot = bad = 0
    for depth in (2, 3):
        for lens in product(range(0, 3), repeat=depth):
            for words in product(*[list(product((0, 1), repeat=l)) for l in lens]):
                for chains in product((0, 1, 2), repeat=depth - 1):
                    P0, _ = word_PQ(words[0])
                    state = (g_fun(P0), 1 + P0 // 2, (P0 + 1) // 2)
                    for j in range(1, depth):
                        state = d1_step(state, words[j], chains[j - 1])
                    m, adj, apexes = build(words, chains)
                    inv = invariants(m, adj, apexes[-1])
                    tot += 1
                    bad += (state != (inv["H"], inv["o"], inv["e"]))
    print("    exhaustive: depth 2 and 3, words of length <= 2, chain lengths 0, 1, 2: %d cases, %d failures" % (tot, bad))


def k4():
    bad = tot = 0
    for p1 in range(0, 40):
        for p2 in range(0, 40):
            for i1 in (1, 2):
                for i2 in (1, 2, 3):
                    tot += 1
                    bad += (matmul(M(p2, i2), M(2 * p1 + 1, i1)) != matmul(M(p2, i2 + 1), M(p1, i1)))
    print("K4  M(p2, i2) M(2 p1 + 1, i1) = M(p2, i2 + 1) M(p1, i1): %d checked, %d failures (the same tournament has several tuples)" % (tot, bad))


if __name__ == "__main__":
    k1()
    k2()
    k3()
    k4()
