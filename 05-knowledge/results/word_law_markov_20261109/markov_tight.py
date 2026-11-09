"""Tight colourings as a finite Markov chain (THM-5021), opus-2026-11-09-S88.

A word w of length L is read left to right. The greedy packing of each colour is a finite automaton:
  s_A, s_B in {1, 2, 3} = distance from the next position to the last chosen element of that colour, capped at 3
  (3 = "a new element of this colour may be chosen");
  E = c - ceil(m/3), where c = number of chosen elements so far and m = number of positions read.
The greedy is optimal for each colour, so after L steps c = P(A) + P(B), and w is tight iff E = 0.
Prefix bound: P(A) + P(B) >= c(m) + ceil((L - m - 2)/3), so a tight word has E <= 1 on every prefix. Hence the
states (s_A, s_B, E <= 1) and the phase m mod 3 give a finite automaton for the tight colourings.
The script builds the one-step transitions, the block transfer matrix T (three steps from phase 0 to phase 0),
exact counts, the rational generating functions F_r(x) = sum_k count(L = 3k - r) x^k (r = 0, 1, 2), and checks the
closed forms 2^k, (k + 1) 2^(k-1), (k^2 + 9k - 2) 2^(k-3) (both values of the first letter counted).
It also records the spectrum of T: the eigenvalue 2 (two independent colour choices per block) and its Jordan
structure, which produces the polynomial factors.
usage: python markov_tight.py KMAX
"""
import sys
from fractions import Fraction
import sympy as sp


def step(state, col, phase):
    sA, sB, E = state
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
    return (sA, sB, E)


def counts_bruteforce(L):
    """number of tight colourings of length L by direct greedy (both first letters)"""
    import itertools

    def P(xs):
        c, last = 0, -10
        for x in xs:
            if x - last >= 3:
                c += 1
                last = x
        return c
    k = -(-L // 3)
    return sum(1 for bits in itertools.product((0, 1), repeat=L)
               if P([i for i in range(L) if bits[i]]) + P([i for i in range(L) if not bits[i]]) == k)


def main():
    KMAX = int(sys.argv[1])
    start = (3, 3, 0)
    states = sorted({(a, b, e) for a in (1, 2, 3) for b in (1, 2, 3) for e in (0, 1)})
    idx = {s: i for i, s in enumerate(states)}
    S = len(states)

    def one_step(phase):
        M = sp.zeros(S, S)
        for s in states:
            for col in (0, 1):
                t = step(s, col, phase)
                if t[2] in (0, 1):            # prune E >= 2 (never tight) and E < 0 (impossible)
                    M[idx[s], idx[t]] += 1
        return M
    M0, M1, M2 = one_step(0), one_step(1), one_step(2)
    T = M0 * M1 * M2
    v0 = sp.zeros(1, S)
    v0[idx[start]] = 1
    fin = sp.zeros(S, 1)
    for s in states:
        if s[2] == 0:
            fin[idx[s]] = 1
    # exact counts by the automaton
    def count(L):
        v = v0
        for m in range(L):
            v = v * (M0 if m % 3 == 0 else M1 if m % 3 == 1 else M2)
        return int((v * fin)[0])
    ok_bf = all(count(L) == counts_bruteforce(L) for L in range(1, 17))
    print(f"automaton states {S}; automaton counts agree with brute force for L = 1..16: {ok_bf}")
    closed = {0: lambda k: 2 ** k, 2: lambda k: Fraction((k + 1) * 2 ** k, 2), 1: lambda k: Fraction((k * k + 9 * k - 2) * 2 ** k, 8)}
    names = {0: "L = 3k (n = 3k + 2)", 2: "L = 3k - 1 (n = 3k + 1)", 1: "L = 3k - 2 (n = 3k)"}
    for r in (0, 2, 1):
        good = True
        vals = []
        for k in range(1, KMAX + 1):
            L = 3 * k - (0 if r == 0 else (1 if r == 2 else 2))
            c = count(L)
            vals.append(c)
            if Fraction(c) != closed[r](k):
                good = False
        print(f"{names[r]}: counts k = 1..{KMAX}: {vals[:10]}...; closed form holds for all k <= {KMAX}: {good}")
    # generating functions F(x) = sum_k count x^k via the block transfer matrix
    x = sp.symbols('x')
    I = sp.eye(S)
    R = (I - x * T).inv()      # exact rational matrix
    gf = {}
    gf[0] = sp.simplify((v0 * R * T * fin)[0])                                  # L = 3k, k >= 1: v0 T^k fin
    gf[2] = sp.simplify((v0 * R * M0 * M1 * fin)[0] * x)                        # L = 3k - 1: v0 T^(k-1) M0 M1 fin
    gf[1] = sp.simplify((v0 * R * M0 * fin)[0] * x)                             # L = 3k - 2: v0 T^(k-1) M0 fin
    for r in (0, 2, 1):
        print(f"generating function, {names[r]}: {sp.factor(gf[r])}")
    print("characteristic polynomial of T:", sp.factor(T.charpoly(x).as_expr()))
    Jform = T.jordan_form(calc_transform=False)
    blocks = []
    J = Jform
    i = 0
    while i < S:
        lam = J[i, i]
        size = 1
        while i + size < S and J[i + size - 1, i + size] == 1:
            size += 1
        blocks.append((lam, size))
        i += size
    print("Jordan blocks of T (eigenvalue, size):", sorted(blocks, key=lambda b: (-abs(complex(b[0])), -b[1])))


if __name__ == "__main__":
    main()
