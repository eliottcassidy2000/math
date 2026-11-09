"""Which of THM-5021's 18 states (s_A, s_B, E) are ever used, and the spectrum of the part that is.
Reuses the auditor's automaton from aud_markov.py; prints states in the session's coordinates
(capped distance s = 3, 1, 2 for the DP difference states (0,0), (1,0), (0,1)).
usage: python aud_markov_states.py
"""
import sympy as sp
from aud_markov import STATES, IDX, START, letter, one_step, jordan_sizes, x

S_OF = {(0, 0): 3, (1, 0): 1, (0, 1): 2}


def show(s):
    return (S_OF[s[0]], S_OF[s[1]], s[2])


def main():
    by_phase = {0: set(), 1: set(), 2: set()}
    layer = {START}
    for m in range(90):
        by_phase[m % 3] |= layer
        nl = set()
        for s in layer:
            for col in (0, 1):
                t = letter(s, col, m)
                if t[2] in (0, 1):
                    nl.add(t)
        layer = nl
    allr = by_phase[0] | by_phase[1] | by_phase[2]
    print("states used at positions m = 0 (mod 3) (block boundaries):", sorted(show(s) for s in by_phase[0]))
    print("states used at some position:", len(allr), "of", len(STATES))
    print("never used:", sorted(show(s) for s in STATES if s not in allr))
    M0, M1, M2 = one_step(0), one_step(1), one_step(2)
    T = M0 * M1 * M2
    used = [s for s in STATES if s in by_phase[0]]
    Tu = sp.Matrix([[T[IDX[a], IDX[b]] for b in used] for a in used])
    print("T restricted to the block-boundary states:", sp.factor(Tu.charpoly(x).as_expr()),
          "; Jordan sizes at 2:", jordan_sizes(Tu, 2), "; at 0:", jordan_sizes(Tu, 0))
    # the generalised eigenvector chain of length 4 in the 18-state T: which states carry it?
    A = T - 2 * sp.eye(len(STATES))
    ker4 = (A ** 4).nullspace()
    ker3 = (A ** 3).nullspace()
    print("dim ker (T - 2)^3 =", len(ker3), ", dim ker (T - 2)^4 =", len(ker4))
    # a vector in ker^4 not in ker^3, and its support
    for v in ker4:
        if (A ** 3 * v) != sp.zeros(len(STATES), 1):
            supp = [show(STATES[i]) for i in range(len(STATES)) if v[i] != 0]
            print("a top vector of the size-4 chain is supported on:", supp)
            break
    # left (row) side: the generating functions use v0 (I - xT)^-1, i.e. the row vector side.
    B = (T - 2 * sp.eye(len(STATES))).T
    rk4 = (B ** 4).nullspace()
    for v in rk4:
        if (B ** 3 * v) != sp.zeros(len(STATES), 1):
            supp = [show(STATES[i]) for i in range(len(STATES)) if v[i] != 0]
            print("a top LEFT vector of the size-4 chain is supported on:", supp)
            break


if __name__ == "__main__":
    main()


def trim_check():
    """the trim block automaton: block-boundary states that are reachable AND can still reach acceptance"""
    M0, M1, M2 = one_step(0), one_step(1), one_step(2)
    T = M0 * M1 * M2
    fin = sp.Matrix([1 if s[2] == 0 else 0 for s in STATES])
    tails = [fin, M0 * M1 * fin, M0 * fin]
    reach, fr = {START}, [START]
    while fr:
        nf = []
        for s in fr:
            for t in STATES:
                if T[IDX[s], IDX[t]] and t not in reach:
                    reach.add(t)
                    nf.append(t)
        fr = nf
    co = set()
    vecs = list(tails)
    for _ in range(20):
        for v in vecs:
            co |= {s for s in STATES if v[IDX[s]] != 0}
        vecs = [T * v for v in vecs]
    dead = sorted(show(s) for s in reach - co)
    trim = [s for s in STATES if s in reach and s in co]
    print("dead reachable block-boundary states (never accepting):", dead,
          "; self-loop weight of (3,3,1):", T[IDX[((0, 0), (0, 0), 1)], IDX[((0, 0), (0, 0), 1)]])
    print("trim states:", [show(s) for s in trim])
    sub = lambda Mx: sp.Matrix([[Mx[IDX[a], IDX[b]] for b in trim] for a in trim])
    Tt = sub(T)
    print("trim T =", Tt.tolist())
    print("trim T: characteristic polynomial", sp.factor(Tt.charpoly(x).as_expr()),
          "; Jordan sizes at 2:", jordan_sizes(Tt, 2), "; at 0:", jordan_sizes(Tt, 0))
    v0 = sp.zeros(1, len(trim))
    v0[trim.index(START)] = 1
    R = (sp.eye(len(trim)) - x * Tt).inv()
    # tails restricted: compute full tail vectors then restrict (the tail letters may pass through untrimmed states)
    for name, tv in (("L = 3k", tails[0]), ("L = 3k - 1", tails[1]), ("L = 3k - 2", tails[2])):
        tr = sp.Matrix([tv[IDX[s]] for s in trim])
        g = (v0 * R * (Tt * tr if name == "L = 3k" else tr))[0] * x
        print(f"  trim generating function {name}: {sp.factor(sp.simplify(g))}")


trim_check()
