"""
anti_redei_flip_parity_opus_S140.py

THE GF(2)-DERIVATIVE ROUTE to the Main Lemma (h_beta odd), following the structure:
  * beta-symmetric tournaments (fixed beta = rho) form a GF(2) cube with free bits:
      - INTERNAL: orientation of each pair {a, beta a}            (h bits, n=2h)
      - S-TIE:    {a->b  with  beta b -> beta a}                  (tied double-arc flips)
      - T-TIE:    {a->beta b  with  b -> beta a}
  * PROVED (mirror-position argument): a beta-symmetric Ham path uses an internal pair-arc
    only as its CENTER arc (exactly once), and uses each cross arc in a mirrored pair with
    its tie-partner.
  * Hence flipping an INTERNAL bit of pair A changes h_beta by W_{beta a} - W_a, where
    W_x = #(half-paths: one representative per pair != A, T-path, ending -> x appended).
    So internal-flip parity  <=>  MIRROR PARITY LEMMA:  W_a == W_{beta a}  (mod 2).
  * Cross flips (s/t ties) need their own pairing; tested here.
  * If all three flip types preserve h_beta mod 2: cube connectivity + base point
    (transitive T, h_rho = 1) prove the Main Lemma, hence Anti-Redei (with Lemma A).

TESTS: random beta-symmetric tournaments (n = 4, 6, 8), all flip types, exact Delta parity;
the mirror-parity lemma W_a vs W_{beta a} directly; n odd (fixed-point) variant probes.
"""
import sys, random, time
from itertools import permutations

sys.path.insert(0, ".")
from metagraph_anti_redei_test_opus_S139 import ham_paths, is_antiauto

def rand_sym_tournament(n, rng):
    beta = [n - 1 - i for i in range(n)]
    a = [[-1] * n for _ in range(n)]
    for i in range(n):
        a[i][i] = 0
    if n % 2 == 1:
        pass  # middle vertex fixed by beta
    for i in range(n):
        if i < beta[i]:
            b = rng.getrandbits(1)
            a[i][beta[i]] = b; a[beta[i]][i] = 1 - b
    for i in range(n):
        for j in range(n):
            if i == j or a[i][j] != -1: continue
            b = rng.getrandbits(1)
            a[i][j] = b; a[j][i] = 1 - b
            bi, bj = beta[i], beta[j]
            if a[bj][bi] == -1:
                a[bj][bi] = b; a[bi][bj] = 1 - b
    assert is_antiauto(a, n, beta), "construction"
    return a, beta

def hbeta(a, n, beta):
    cnt = 0
    for P in ham_paths(a, n):
        if all(beta[P[j]] == P[n - 1 - j] for j in range(n)):
            cnt += 1
    return cnt

def W_end(a, n, beta, A_pair, x):
    """#(sequences v_1..v_{h-1}, one rep per pair != A, T-path, with v_{h-1} -> x);
       for h=1 (no other pairs) define = 1 (empty path)."""
    pairs = []
    seen = set()
    for i in range(n):
        p = frozenset((i, beta[i]))
        if p not in seen and p != frozenset(A_pair):
            seen.add(p); pairs.append(tuple(sorted(p)))
    hh = len(pairs)
    if hh == 0:
        return 1
    cnt = 0
    def ext(seq, used):
        nonlocal cnt
        if len(seq) == hh:
            if a[seq[-1]][x]: cnt += 1
            return
        for pi, pr in enumerate(pairs):
            if used[pi]: continue
            for v in pr:
                if not seq or a[seq[-1]][v]:
                    used[pi] = True; seq.append(v)
                    ext(seq, used)
                    seq.pop(); used[pi] = False
    ext([], [False] * hh)
    return cnt

def main():
    rng = random.Random(140)
    print("=" * 96)
    print("flip-parity tests on beta-symmetric tournaments")
    print("=" * 96)
    stats = {"internal": [0, 0], "s-tie": [0, 0], "t-tie": [0, 0]}
    mirror_bad = mirror_tested = 0
    for trial in range(600):
        n = rng.choice([4, 6, 8])
        a, beta = rand_sym_tournament(n, rng)
        h0 = hbeta(a, n, beta)
        # INTERNAL flips
        for i in range(n // 2):
            u, bu = i, beta[i]
            a[u][bu] ^= 1; a[bu][u] ^= 1
            h1 = hbeta(a, n, beta)
            a[u][bu] ^= 1; a[bu][u] ^= 1
            stats["internal"][(h1 - h0) % 2] += 1
            # mirror parity lemma
            Wa = W_end(a, n, beta, (u, bu), u)
            Wb = W_end(a, n, beta, (u, bu), bu)
            mirror_tested += 1
            if (Wa - Wb) % 2: mirror_bad += 1
        # one S-TIE and one T-TIE flip per trial
        i, j = rng.sample(range(n // 2), 2) if n >= 4 else (0, 1)
        u, v = i, j
        bu, bv = beta[u], beta[v]
        # s-tie: arcs (u,v) and (bv,bu)
        for (p, q, r, s) in (((u, v), None, None, None),):
            pass
        def flip_pair(x, y):
            a[x][y] ^= 1; a[y][x] ^= 1
        # s-tie
        flip_pair(u, v); flip_pair(bv, bu)
        if is_antiauto(a, n, beta):
            h1 = hbeta(a, n, beta)
            stats["s-tie"][(h1 - h0) % 2] += 1
        flip_pair(u, v); flip_pair(bv, bu)
        # t-tie
        flip_pair(u, bv); flip_pair(v, bu)
        if is_antiauto(a, n, beta):
            h1 = hbeta(a, n, beta)
            stats["t-tie"][(h1 - h0) % 2] += 1
        flip_pair(u, bv); flip_pair(v, bu)
    for k, (ev, od) in stats.items():
        print(f"  {k:9s}: Delta even {ev}, Delta ODD {od}   "
              f"{'PARITY-SAFE' if od == 0 else '*** PARITY-VIOLATING FLIPS EXIST ***'}")
    print(f"  MIRROR PARITY LEMMA (W_a == W_beta_a mod 2): tested {mirror_tested}, "
          f"violations {mirror_bad}  {'HOLDS' if mirror_bad == 0 else '*** FAILS ***'}")
    # odd n probe: fixed-point variant, internal flips only exist for pairs
    print("\n  odd-n probe (n=5,7): h_beta odd on random symmetric tournaments:")
    bad = tested = 0
    for trial in range(400):
        n = rng.choice([5, 7])
        a, beta = rand_sym_tournament(n, rng)
        if hbeta(a, n, beta) % 2 == 0: bad += 1
        tested += 1
    print(f"    tested {tested}, even h_beta: {bad}  {'(Main Lemma consistent)' if bad == 0 else '*** COUNTEREXAMPLE ***'}")

if __name__ == "__main__":
    t0 = time.time()
    main()
    print(f"[{time.time()-t0:.0f}s]")
