"""Initial (leading-monomial) ideal of the arc ideal of R_q (wt x=2, y=3) under a lex order on jet variables,
weight by weight; prints minimal monomial generators of the initial ideal (as jet-variable multisets)."""
import sys, functools
import flint
from arc_hs import monos, mul, deriv, poly, RELS, WT, P

def order_key(order, W=64):
    """True lex monomial order. Variables ranked (largest first) by weight ascending; ties: x first ('wtx') or y
    first ('wty'); key(m) = exponent vector in that variable order (lexicographic comparison = lex order)."""
    allv = [(n, j) for n in ('x', 'y') for j in range(W)]
    if order == 'wtx':
        allv.sort(key=lambda v: (WT[v[0]] + v[1], 0 if v[0] == 'x' else 1))
    elif order == 'wty':
        allv.sort(key=lambda v: (WT[v[0]] + v[1], 0 if v[0] == 'y' else 1))
    elif order == 'rwtx':   # largest variable = highest weight
        allv.sort(key=lambda v: (-(WT[v[0]] + v[1]), 0 if v[0] == 'x' else 1))
    pos = {v: i for i, v in enumerate(allv)}
    def key(m):
        vec = [0] * len(allv)
        for v in m: vec[pos[v]] += 1
        return tuple(vec)
    return key

def initial_ideal_gens(rels, W, order):
    gens = [poly(r) for r in rels]
    gw = [sum(WT[v] + j for (v, j) in next(iter(g))) for g in gens]
    ders = []
    for g, w in zip(gens, gw):
        lst = [g]
        while w + len(lst) - 1 < W: lst.append(deriv(lst[-1]))
        ders.append((w, lst))
    key = order_key(order)
    lead_all = []
    for w in range(0, W + 1):
        M = monos(w)
        # order columns so that "larger" monomials come first: lex descending
        M = sorted(M, key=key, reverse=True)
        idx = {m: i for i, m in enumerate(M)}
        rows = []
        for (gw0, lst) in ders:
            for n, dg in enumerate(lst):
                rest = w - gw0 - n
                if rest < 0: break
                for m in monos(rest):
                    row = [0] * len(M)
                    for mo, co in dg.items():
                        k = idx[mul(mo, m)]; row[k] = (row[k] + co) % P
                    rows.append(row)
        if not rows: continue
        R, rk = flint.nmod_mat(len(rows), len(M), [x for r in rows for x in r], P).rref()
        for r in R.tolist()[:rk]:
            piv = next(j for j, x in enumerate(r) if int(x))
            lead_all.append(M[piv])
    # minimal generators: leading monomials not divisible by an earlier (lower-weight) leading monomial
    def divides(a, b):
        from collections import Counter
        ca, cb = Counter(a), Counter(b)
        return all(cb[v] >= c for v, c in ca.items())
    mins = []
    for m in lead_all:
        if not any(divides(g, m) for g in mins):
            mins.append(m)
    return mins

if __name__ == '__main__':
    q, W, order = int(sys.argv[1]), int(sys.argv[2]), sys.argv[3]
    mins = initial_ideal_gens(RELS[q], W, order)
    def fmt(m):
        from collections import Counter
        c = Counter(m)
        return "*".join(f"{v}{j}" + (f"^{e}" if e > 1 else "") for (v, j), e in sorted(c.items(), key=lambda t: (t[0][1], t[0][0])))
    print(f"q={q} order={order}: {len(mins)} minimal generators of the initial ideal up to weight {W}")
    for m in mins:
        print("  wt", sum(WT[v] + j for v, j in m), ":", fmt(m))
