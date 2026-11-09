"""S19 audit: brute-force check of THM-5085 (1) (transfer), (2) (closed form and final port state) in bases 2, 3, 4.

Builds the actual tournaments (core + word over {S, B}, port u fixed), counts everything with tbrute (validated
against a naive counter by naive_check.py) and compares with
  (a) the step recursion S:(h,R,D)->(h,2R,0), B:(h,R,D)->(h + x(R+D)/2, 2R - hu, -hu) started from the core's
      brute-force state;
  (b) the closed form h = h0 + (x/2)[R0 N(I) + D0 [1 in I] - hu Q(I)] and the stated final port state
      R = 2^L R0 - hu sum_{i in I} 2^(L-i), D = -hu [L in I] (0 if last letter S, D0 if L = 0);
  (c) hu invariance (I(Omega(T - u), x) of the final tournament equals that of the core).
Q(I) is computed as the pair sum of THM-5085 and, independently, by Session 205's recursion
Q(0)=0, Q(2m)=Q(m), Q(2m+1)=Q(m)+m+(m mod 2); the two are compared.
Usage: python check_transfer.py TBRUTE_EXE SEED NRANDOM
"""
import itertools, random, subprocess, sys
from fractions import Fraction

XS = (2, 3, 4)

def all_tournaments(n):
    pairs = [(i, j) for i in range(n) for j in range(i + 1, n)]
    for bits in range(1 << len(pairs)):
        out = [0] * n
        for k, (i, j) in enumerate(pairs):
            if bits >> k & 1: out[i] |= 1 << j
            else: out[j] |= 1 << i
        yield out

def rand_tour(n, rng):
    out = [0] * n
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < 0.5: out[i] |= 1 << j
            else: out[j] |= 1 << i
    return out

def extend(core, u, word):
    out = list(core)
    for letter in word:
        z = len(out)
        for v in range(z):
            if letter == 'S' or v != u:
                out[v] |= 1 << z
        out.append(0 if letter == 'S' else (1 << u))
    return out

def Q_pairs(I):
    I = sorted(I); q = 0
    for a in range(len(I)):
        for b in range(a + 1, len(I)):
            d = I[b] - I[a]
            q += 2 if d == 1 else 2 ** (d - 1)
    return q

def Q_rec(m):
    if m == 0: return 0
    if m % 2 == 0: return Q_rec(m // 2)
    k = m // 2
    return Q_rec(k) + k + (k % 2)

def predict_steps(state, hu, word, x):
    h, R, D = state
    for letter in word:
        if letter == 'S': R, D = 2 * R, 0
        else: h, R, D = h + Fraction(x, 2) * (R + D), 2 * R - hu, -hu
    return h, R, D

def predict_closed(state, hu, word, x):
    h0, R0, D0 = state
    L = len(word)
    I = [i + 1 for i, c in enumerate(word) if c == 'B']
    N = sum(2 ** (i - 1) for i in I)
    Q = Q_pairs(I)
    assert Q == Q_rec(N), (I, Q, Q_rec(N))
    h = h0 + Fraction(x, 2) * (R0 * N + D0 * (1 if 1 in I else 0) - hu * Q)
    R = 2 ** L * R0 - hu * sum(2 ** (L - i) for i in I)
    if L == 0: D = D0
    elif word[-1] == 'S': D = 0
    else: D = -hu
    return h, R, D

def run(exe, items):
    inp = "".join(f"{len(out)} {u} " + " ".join(map(str, out)) + "\n" for (out, u) in items)
    res = subprocess.run([exe], input=inp, capture_output=True, text=True)
    if res.returncode != 0: raise SystemExit("tbrute failed: " + res.stderr)
    return [list(map(int, l.split())) for l in res.stdout.strip().split("\n")]

def stats_of(line, x):
    k = 2 + 4 * (x - 2)
    I, E, O, hu = line[k:k + 4]
    return I, E + O, E - O, hu

def main():
    exe = sys.argv[1]; seed = int(sys.argv[2]); nrand = int(sys.argv[3])
    rng = random.Random(seed)
    cases = []  # (core, u, word)
    # exhaustive: all labelled cores of order 1..4, every port, every word of length <= 5 (order <= 3: <= 6)
    for n in range(1, 5):
        maxL = 6 if n <= 3 else 5
        for core in all_tournaments(n):
            for u in range(n):
                for L in range(0, maxL + 1):
                    for word in itertools.product('SB', repeat=L):
                        cases.append((core, u, ''.join(word)))
    nex = len(cases)
    # the example of THM-5085 (3): core {u}, S^s B B
    for s in range(0, 10):
        cases.append(([0], 0, 'S' * s + 'BB'))
    # random: cores of order 5..7, random port, random words, total order <= 14
    for _ in range(nrand):
        n = rng.randint(5, 7)
        core = rand_tour(n, rng)
        u = rng.randrange(n)
        L = rng.randint(1, 14 - n)
        word = ''.join(rng.choice('SB') for _ in range(L))
        cases.append((core, u, word))
    cores = {}
    for core, u, word in cases: cores[(tuple(core), u)] = None
    core_list = list(cores)
    core_res = run(exe, [(list(c), u) for (c, u) in core_list])
    for key, line in zip(core_list, core_res): cores[key] = line
    tails = run(exe, [(extend(core, u, word), u) for (core, u, word) in cases])
    bad_step = bad_closed = bad_hu = bad_final = bad_ocf = 0
    example_lines = []
    for idx, ((core, u, word), line) in enumerate(zip(cases, tails)):
        cl = cores[(tuple(core), u)]
        if line[0] != line[2]: bad_ocf += 1
        for x in XS:
            h0, R0, D0, hu = stats_of(cl, x)
            h, R, D, hu2 = stats_of(line, x)
            if hu2 != hu: bad_hu += 1
            ps = predict_steps((h0, R0, D0), hu, word, x)
            pc = predict_closed((h0, R0, D0), hu, word, x)
            if ps != (h, R, D): bad_step += 1
            if pc[0] != h: bad_closed += 1
            if (pc[1], pc[2]) != (R, D): bad_final += 1
        if nex <= idx < nex + 10:
            s = len(word) - 2
            example_lines.append(f"  core {{u}}, S^{s} B B: H = {line[0]};  3*2^s - 1 = {3 * 2 ** s - 1};  3*2^(s+1) - 1 = {3 * 2 ** (s + 1) - 1}")
    print(f"instances: {len(cases)} ({nex} exhaustive: all labelled cores of order 1..4, all ports, all words of length <= 6/5;"
          f" 10 examples; {nrand} random cores of order 5..7 with words to total order 14)")
    print(f"max order {max(len(core) + len(word) for core, u, word in cases)}; bases x = {XS}")
    print(f"step-recursion mismatches (h,R,Delta): {bad_step}")
    print(f"closed-form h mismatches: {bad_closed}")
    print(f"closed-form final (R,Delta) mismatches: {bad_final}")
    print(f"hu invariance failures: {bad_hu}")
    print(f"OCF failures (H != I(Omega,2)): {bad_ocf}")
    print("THM-5085 (3) example (x = 2):")
    print("\n".join(example_lines))

if __name__ == "__main__":
    main()
