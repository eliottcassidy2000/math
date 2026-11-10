# a08: THM-5109 (F).
#  (a) physical check of the fixed-port letters (THM-5085, base 2) on random cores: S = append a sink,
#      B = append w beaten by all except the port u, with w -> u. State (h, R, Delta) = (H, E+O, E-O) at the port u,
#      Delta = -e (e = Hamiltonian paths ending at u). Check S: (h, 2R, 0), B: (h + R + Delta, 2R - c, -c), c = H(T0 - u),
#      and the post-bridge maps T_0 = B, T_j = S^j B as stated in PROOF.md section 6.
#  (b) lattice tiling: for several odd c, every (h', R') with h', R' odd and R' + c != 0 has exactly one preimage under
#      {T_j : j >= 0}; the point R' = -c has none (outside the physical range R' > 0).
import random
from a_lib import *
rng = random.Random(77)

def add_sink(A):
    n = len(A); B = [row + [1] for row in A]; B.append([0] * (n + 1)); return B
def add_bridge(A, u):
    n = len(A); B = [row + [0 if i == u else 1] for i, row in enumerate(A)]
    last = [0] * (n + 1); last[u] = 1; B.append(last); return B   # w -> u, every other vertex -> w

bad_letters = 0; checked = 0; bad_post = 0; post_checked = 0
for trial in range(40):
    n0 = rng.randrange(1, 6); T0 = rand_tour(n0, rng); u = rng.randrange(n0)
    Tu = [[T0[i][j] for j in range(n0) if j != u] for i in range(n0) if i != u]
    c = run([Tu])[0]['H'] if n0 > 1 else 1
    word = [rng.choice('SB') for _ in range(rng.randrange(1, 7))]
    A = T0; r = run([A])[0]; st = (r['H'], r['R'][u], -r['e'][u])
    for L in word:
        A2 = add_sink(A) if L == 'S' else add_bridge(A, u)
        r2 = run([A2])[0]; st2 = (r2['H'], r2['R'][u], -r2['e'][u])
        h, R, D = st
        pred = (h, 2 * R, 0) if L == 'S' else (h + R + D, 2 * R - c, -c)
        checked += 1
        if pred != st2: bad_letters += 1
        A, st = A2, st2
    # post-bridge maps: from a state with Delta = -c, apply S^j B physically and compare with T_j
    if st[2] == -c:
        for j in range(0, 4):
            A2 = A
            for _ in range(j): A2 = add_sink(A2)
            A2 = add_bridge(A2, u)
            r2 = run([A2])[0]; got = (r2['H'], r2['R'][u])
            h, R, _ = st
            exp = (h + R - c, 2 * R - c) if j == 0 else (h + 2 ** j * R, 2 ** (j + 1) * R - c)
            post_checked += 1
            if got != exp: bad_post += 1
print('(a) letters S/B checked on %d physical steps: mismatches %d; post-bridge T_j checked on %d: mismatches %d'
      % (checked, bad_letters, post_checked, bad_post))

def preimages(hp, Rp, c):
    out = []
    j = 0
    while 2 ** (j + 1) <= abs(Rp + c) and Rp + c != 0:
        if (Rp + c) % 2 ** (j + 1) == 0:
            R = (Rp + c) // 2 ** (j + 1)
            h = hp - R + c if j == 0 else hp - 2 ** j * R
            if R % 2 and h % 2 and ((h + R - c, 2 * R - c) if j == 0 else (h + 2 ** j * R, 2 ** (j + 1) * R - c)) == (hp, Rp):
                out.append((j, h, R))
        j += 1
    return out
for c in (1, 3, 5, 9, 15, 45):
    cnt = collections_bad = 0; zero = []
    for hp in range(-101, 102, 2):
        for Rp in range(-301, 302, 2):
            pre = preimages(hp, Rp, c)
            if Rp + c == 0:
                zero.append(len(pre)); continue
            cnt += 1
            if len(pre) != 1: collections_bad += 1
    print('(b) c=%d: %d lattice points (h, R odd, |h| <= 101, |R| <= 301, R != -c): points without exactly one preimage: %d; '
          'points with R = -c have %s preimages' % (c, cnt, collections_bad, sorted(set(zero))))
