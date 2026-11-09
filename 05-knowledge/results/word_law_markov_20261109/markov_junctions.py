"""THM-5021 (5): tight colourings at n = 3k (L = 3k - 2) that contain a junction, counted exactly by the greedy
automaton of markov_tight.py extended with the last four letters and a junction flag (factor 11011 = isolated-0
junction, factor 00100 = isolated-1 junction). Exact generating functions by the block transfer matrix (sympy), and a
comparison with brute force for small k.
usage: python markov_junctions.py KMAX
"""
import sys
import itertools
import sympy as sp


def step(state, col, phase):
    sA, sB, E, last, j0, j1 = state
    chosen = 0
    if col == 1:
        if sA >= 3:
            chosen, sA = 1, 1
        else:
            sA = min(3, sA + 1)
        sB = min(3, sB + 1)
    else:
        if sB >= 3:
            chosen, sB = 1, 1
        else:
            sB = min(3, sB + 1)
        sA = min(3, sA + 1)
    E = E + chosen - (1 if phase == 0 else 0)
    s = last + str(col)
    j0 = j0 or s[-5:] == '11011'
    j1 = j1 or s[-5:] == '00100'
    return (sA, sB, E, s[-4:], j0, j1)


def main():
    KMAX = int(sys.argv[1])
    start = (3, 3, 0, '', False, False)
    # reachable states
    states, frontier = {start}, [start]
    while frontier:
        nxt = []
        for s in frontier:
            for ph in (0, 1, 2):
                for col in (0, 1):
                    t = step(s, col, ph)
                    if t[2] in (0, 1) and t not in states:
                        states.add(t)
                        nxt.append(t)
        frontier = nxt
    states = sorted(states, key=str)
    idx = {s: i for i, s in enumerate(states)}
    S = len(states)

    def one(ph):
        M = sp.zeros(S, S)
        for s in states:
            for col in (0, 1):
                t = step(s, col, ph)
                if t[2] in (0, 1):
                    M[idx[s], idx[t]] += 1
        return M
    M0, M1, M2 = one(0), one(1), one(2)
    T = M0 * M1 * M2
    v0 = sp.zeros(1, S)
    v0[idx[start]] = 1

    def final(pred):
        f = sp.zeros(S, 1)
        for s in states:
            if s[2] == 0 and pred(s):
                f[idx[s]] = 1
        return f
    x = sp.symbols('x')
    R = (sp.eye(S) - x * T).inv()
    preds = {'0-junction': lambda s: s[4], '1-junction': lambda s: s[5], 'any junction': lambda s: s[4] or s[5],
             'both kinds': lambda s: s[4] and s[5]}
    print(f"automaton states (with last four letters and junction flags): {S}")
    for name, pr in preds.items():
        g = sp.factor(sp.simplify((v0 * R * M0 * final(pr))[0] * x))    # L = 3k - 2: v0 T^(k-1) M0
        ser = sp.series(g, x, 0, KMAX + 1).removeO()
        coeffs = [int(ser.coeff(x, k)) for k in range(1, KMAX + 1)]
        print(f"tight colourings at n = 3k with {name}: generating function {g}; k = 1..{KMAX}: {coeffs}")

    # brute-force comparison for small k
    def P(xs):
        c, last = 0, -10
        for xx in xs:
            if xx - last >= 3:
                c += 1
                last = xx
        return c
    for k in range(1, 7):
        L = 3 * k - 2
        c0 = c1 = 0
        for bits in itertools.product('01', repeat=L):
            w = ''.join(bits)
            A = [i for i in range(L) if w[i] == '1']
            B = [i for i in range(L) if w[i] == '0']
            if P(A) + P(B) == k:
                c0 += '11011' in w
                c1 += '00100' in w
        print(f"brute force k = {k}: with 0-junction {c0}, with 1-junction {c1}")


if __name__ == "__main__":
    main()
