"""Independent audit (auditor) of THM-5021: the tight-colouring automaton, its block transfer matrix, the
generating functions and the junction counts.

Packing DP per colour: f(i+1) = max(f(i), [i in X](1 + f(i-2))).  Its state is the difference vector
(d1, d2) = (f(i) - f(i-1), f(i-1) - f(i-2)) in {(0,0), (1,0), (0,1)}; a letter of that colour raises f iff
(d1, d2) = (0, 0).  (This is the DP, derived without the greedy rule; (0,0), (1,0), (0,1) correspond to the
session's capped distances 3, 1, 2.)  E = (P_A + P_B on the prefix) - ceil(m/3).

Checks:
  (a) the automaton with E in {0, 1} counts tight colourings: compared with brute force (all words, L <= 20);
  (b) the 18-state block matrix exactly as THM-5021 defines it (all 3 x 3 x 2 states): characteristic
      polynomial and Jordan structure, and the same for the trim automaton (states reachable from the start at
      a block boundary and co-reachable to acceptance);
  (c) the three generating functions, from the automaton (sympy, exact), compared with the claimed ones, and the
      closed forms checked for k <= 40;
  (d) the junction automaton (last four letters + flags): generating functions for 11011, 00100, both;
      number of reachable states with the phase tracked correctly versus "any phase at any step".
usage: python aud_markov.py
"""
import itertools
import sympy as sp

x = sp.symbols('x')
DSTATES = [(0, 0), (1, 0), (0, 1)]


def dp_step(d, mine):
    """difference state of one colour after reading a letter; returns (new_d, increased)"""
    d1, d2 = d
    if mine:
        inc = 1 if (d1, d2) == (0, 0) else 0
        return (inc, d1), inc
    return (0, d1), 0


def letter(state, col, m):
    """state = (dA, dB, E); m = number of letters read before this one"""
    dA, dB, E = state
    nA, iA = dp_step(dA, col == 1)
    nB, iB = dp_step(dB, col == 0)
    E2 = E + iA + iB - (1 if m % 3 == 0 else 0)
    return (nA, nB, E2)


def P(xs):
    f = [0, 0, 0]
    for b in xs:
        nf = f[2]
        if b and 1 + f[0] > nf:
            nf = 1 + f[0]
        f = [f[1], f[2], nf]
    return f[2]


def brute(L):
    k = -(-L // 3)
    c = 0
    for bits in itertools.product((0, 1), repeat=L):
        if P(bits) + P([1 - b for b in bits]) == k:
            c += 1
    return c


STATES = [(a, b, e) for a in DSTATES for b in DSTATES for e in (0, 1)]
IDX = {s: i for i, s in enumerate(STATES)}
START = ((0, 0), (0, 0), 0)


def one_step(phase, states=STATES, idx=IDX):
    M = sp.zeros(len(states), len(states))
    for s in states:
        for col in (0, 1):
            t = letter(s, col, phase)
            if t[2] in (0, 1) and t in idx:
                M[idx[s], idx[t]] += 1
    return M


def count_automaton(L):
    v = {START: 1}
    for m in range(L):
        nv = {}
        for s, c in v.items():
            for col in (0, 1):
                t = letter(s, col, m)
                if t[2] in (0, 1):
                    nv[t] = nv.get(t, 0) + c
        v = nv
    return sum(c for s, c in v.items() if s[2] == 0)


def jordan_sizes(T, lam):
    S = T.shape[0]
    A = T - lam * sp.eye(S)
    ranks = [S]
    P_ = sp.eye(S)
    for j in range(1, S + 1):
        P_ = P_ * A
        ranks.append(P_.rank())
        if ranks[-1] == ranks[-2]:
            break
    # number of blocks of size >= j = ranks[j-1] - ranks[j]
    ge = [ranks[j - 1] - ranks[j] for j in range(1, len(ranks))]
    sizes = {}
    for j in range(1, len(ge) + 1):
        exact = ge[j - 1] - (ge[j] if j < len(ge) else 0)
        if exact:
            sizes[j] = exact
    return sizes


def main():
    # (a) automaton versus brute force
    ok = True
    for L in range(1, 21):
        a, b = count_automaton(L), brute(L)
        if a != b:
            ok = False
            print("MISMATCH", L, a, b)
    print(f"(a) automaton count = brute-force count of tight colourings for L = 1..20: {ok}")

    # (b) the 18-state matrix as defined in THM-5021
    M0, M1, M2 = one_step(0), one_step(1), one_step(2)
    T = M0 * M1 * M2
    print("(b) 18-state T: characteristic polynomial", sp.factor(T.charpoly(x).as_expr()))
    print("    Jordan block sizes at eigenvalue 2:", jordan_sizes(T, 2), "; at eigenvalue 0:", jordan_sizes(T, 0))
    # trim automaton at block boundaries: reachable from START via whole blocks, co-reachable to E = 0 at the end of
    # some word of length 3k, 3k - 1 or 3k - 2
    reach = {START}
    frontier = [START]
    while frontier:
        nf = []
        for s in frontier:
            for t_ in STATES:
                if T[IDX[s], IDX[t_]] and t_ not in reach:
                    reach.add(t_)
                    nf.append(t_)
        frontier = nf
    # co-reachable: from s some continuation (whole blocks then 0, 1 or 2 letters) ends with E = 0
    fin = sp.Matrix([1 if s[2] == 0 else 0 for s in STATES])
    tails = [fin, M0 * fin, M0 * M1 * fin]
    co = set()
    vecs = list(tails)
    for _ in range(len(STATES) + 1):
        for vv in vecs:
            for s in STATES:
                if vv[IDX[s]] != 0:
                    co.add(s)
        vecs = [T * vv for vv in vecs]
    trim = [s for s in STATES if s in reach and s in co]
    print(f"    states reachable at block boundaries: {len(reach)}; co-reachable: {len(co)}; trim: {len(trim)}")
    Tt = sp.Matrix([[T[IDX[s], IDX[t_]] for t_ in trim] for s in trim])
    print("    trim T: characteristic polynomial", sp.factor(Tt.charpoly(x).as_expr()),
          "; Jordan sizes at 2:", jordan_sizes(Tt, 2))
    # states reachable at ANY position (phase tracked)
    seen = set()
    layer = {START}
    for m in range(60):
        seen |= layer
        nl = set()
        for s in layer:
            for col in (0, 1):
                t_ = letter(s, col, m)
                if t_[2] in (0, 1):
                    nl.add(t_)
        layer = nl
    print(f"    states reachable at some position (phase tracked): {len(seen)} of 18")

    # (c) generating functions
    S = len(STATES)
    v0 = sp.zeros(1, S)
    v0[IDX[START]] = 1
    R = (sp.eye(S) - x * T).inv()
    g0 = sp.factor(sp.simplify(x * (v0 * R * T * fin)[0]))
    g1 = sp.factor(sp.simplify(x * (v0 * R * M0 * M1 * fin)[0]))
    g2 = sp.factor(sp.simplify(x * (v0 * R * M0 * fin)[0]))
    claims = {'L = 3k': 2 * x / (1 - 2 * x), 'L = 3k - 1': 2 * x * (1 - x) / (1 - 2 * x) ** 2,
              'L = 3k - 2': 2 * x * (1 - x - x ** 2) / (1 - 2 * x) ** 3}
    for name, g in (('L = 3k', g0), ('L = 3k - 1', g1), ('L = 3k - 2', g2)):
        print(f"(c) sum_k N({name}) x^k = {g}; equals the claim: {sp.simplify(g - claims[name]) == 0}")
    print("    note: without the factor x, v0 (I - xT)^(-1) T fin =", sp.factor(sp.simplify((v0 * R * T * fin)[0])))
    # closed forms against the automaton for k <= 40
    good = True
    for k in range(1, 41):
        if count_automaton(3 * k) != 2 ** k:
            good = False
        if count_automaton(3 * k - 1) != (k + 1) * 2 ** (k - 1):
            good = False
        if 8 * count_automaton(3 * k - 2) != (k * k + 9 * k - 2) * 2 ** k:
            good = False
    print(f"    closed forms 2^k, (k+1) 2^(k-1), (k^2+9k-2) 2^(k-3) hold for k = 1..40: {good}")
    # the partial-fraction coefficient identity in the proof
    kk = sp.symbols('k', integer=True, positive=True)
    lhs = 2 * (sp.binomial(kk + 1, 2) * 2 ** (kk - 1) - sp.binomial(kk, 2) * 2 ** (kk - 2)
               - sp.binomial(kk - 1, 2) * 2 ** (kk - 3))
    print("    coefficient identity 2[C(k+1,2)2^(k-1) - C(k,2)2^(k-2) - C(k-1,2)2^(k-3)] = (k^2+9k-2)2^(k-3):",
          sp.simplify(sp.expand_func(lhs) - (kk ** 2 + 9 * kk - 2) * 2 ** (kk - 3)) == 0)
    # series check of the claimed GF of the n = 3k counts
    ser = sp.series(claims['L = 3k - 2'], x, 0, 12).removeO()
    print("    series of 2x(1-x-x^2)/(1-2x)^3:", [ser.coeff(x, j) for j in range(1, 12)])

    # (d) junction automaton
    def jletter(state, col, m):
        base, last, f0, f1 = state
        nb = letter(base, col, m)
        s = last + str(col)
        g0 = f0 or (len(s) >= 5 and s[-5:] == '11011')
        g1 = f1 or (len(s) >= 5 and s[-5:] == '00100')
        return (nb, s[-4:], g0, g1)

    jstart = (START, '', False, False)
    # reachable with the phase tracked (union over positions), and with any phase at any step
    seen, layer = set(), {jstart}
    for m in range(80):
        seen |= layer
        nl = set()
        for s in layer:
            for col in (0, 1):
                t_ = jletter(s, col, m)
                if t_[0][2] in (0, 1):
                    nl.add(t_)
        layer = nl
    anyph, fr = {jstart}, [jstart]
    while fr:
        nf = []
        for s in fr:
            for ph in (0, 1, 2):
                for col in (0, 1):
                    t_ = jletter(s, col, ph)
                    if t_[0][2] in (0, 1) and t_ not in anyph:
                        anyph.add(t_)
                        nf.append(t_)
        fr = nf
    print(f"(d) junction automaton: states reachable with the phase tracked {len(seen)}; "
          f"with any phase at any step {len(anyph)}")
    # counts by exact transfer over positions for L = 3k - 2, k <= 30, then fit-free comparison
    def jcount(L):
        v = {jstart: 1}
        for m in range(L):
            nv = {}
            for s, c in v.items():
                for col in (0, 1):
                    t_ = jletter(s, col, m)
                    if t_[0][2] in (0, 1):
                        nv[t_] = nv.get(t_, 0) + c
            v = nv
        out = {'11011': 0, '00100': 0, 'both': 0}
        for s, c in v.items():
            if s[0][2] == 0:
                if s[2]:
                    out['11011'] += c
                if s[3]:
                    out['00100'] += c
                if s[2] and s[3]:
                    out['both'] += c
        return out
    good = True
    rows = []
    for k in range(1, 31):
        o = jcount(3 * k - 2)
        want = (k - 2) * 2 ** (k - 3) if k >= 3 else 0
        if o['11011'] != want or o['00100'] != want or o['both'] != 0:
            good = False
        if k <= 8:
            rows.append((k, o['11011'], o['00100'], o['both']))
    print(f"    tight colourings at n = 3k with 11011 / 00100 / both, k = 1..8: {rows}")
    print(f"    (k - 2) 2^(k-3) each and 0 for both, for k = 3..30: {good}")
    # and at n = 3k + 1, 3k + 2 there are none (J = 0 there)
    none = all(jcount(L)['11011'] == 0 and jcount(L)['00100'] == 0 for L in range(1, 61) if L % 3 != 1)
    print(f"    no tight colouring with a junction at n != 0 (mod 3), L <= 60: {none}")


if __name__ == "__main__":
    main()
