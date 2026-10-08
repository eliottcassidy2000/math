"""The full lower block of the t = 23 fan bottom and its clusters (session opus-2026-10-07-S21, audit corrections).
Note: 05-knowledge/results/mersenne_line_barriers_20261007.md (section 2).  Uses coins_of/par/xres/stepf from
mersenne_fan_escape_20261007.py.  Source M_99708993677 (the fan bottom), source-reference chains D = 1..4000, run with
collapse of identical states (identical states have identical futures, since every chain reads the same coins).
  (J) all 4000 chains to step 2,741,940: the group of D = 1 is exactly {1..1910} u {1929..1932, 1935, 1936}
      (1916 members; its last merge is at step 2,741,933); the other 2084 chains form two clusters with interleaved
      D ranges; no absorption yet.
  (K) the surviving clusters continued to 19,000,766: the two lower clusters merge at 4,474,989 into one cluster of
      2084 members (level 2); the D = 1 group is absorbed at Terras time 19,000,765 with all 1916 members, so
      M_99708993677 ~> M_(99708993677 - D) for each of them (deepest D = 1936: M_99708991741); the level-2 cluster is
      not absorbed by 19,000,766.  Merge events of the marked chains D = 1, 1910, 1911, 2895, 2902, 4000; levels at
      checkpoints; the maximum level of the D = 1911 chain over [0, 19,000,765].
  (N) odd-step counts over 19,000,765 Terras steps of x_E0, y_1, y_256 from the fast parity vectors.
Prints ALL CHECKS PASSED.  Runtime about 2.5 minutes (single process); memory about 1 GB."""
import sys, os, time
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import mersenne_fan_escape_20261007 as F          # noqa: E402

E0 = F.E0
NABS = F.NABS                                       # 19000765
NJ = 2741940
DMAX = 4000
EXTRA = {1929, 1930, 1931, 1932, 1935, 1936}
MARK = [1, 1910, 1911, 2895, 2902, 4000]
FAIL = []


def check(c, m):
    print(('  ok   ' if c else '  FAIL ') + m, flush=True)
    if not c:
        FAIL.append(m)


if __name__ == '__main__':
    t0 = time.time()
    cb = F.coins_of(E0, NABS + 16).tolist()
    print(f'coins of the fan bottom to {len(cb)} steps [{time.time() - t0:.0f}s]', flush=True)
    alive = {D: (-D, 1 - F.P3[D]) for D in range(1, DMAX + 1)}
    groups = {D: {D} for D in alive}
    rep = {m: m for m in MARK}                      # current representative of each marked chain
    events = []                                     # (step, marked chains of the surviving group, of the joining group)
    absorbed = {}                                   # representative -> absorption time
    last_merge1 = None
    max1911 = -1911; cps = {NJ, 4474989, 8388608, 16777216, NABS}; levels = {}
    merge2 = None
    n = 0
    print('(J)/(K) collapsing search', flush=True)
    while n <= NABS + 1:
        b = cb[n]
        new = {}
        for D, (k, av) in alive.items():
            if k == 0 and av == 0:
                absorbed.setdefault(D, n)
            new[D] = F.stepf(k, av, b)              # (0, 0) is a fixed point of the update
        seen = {}; alive = {}
        for D, st in new.items():
            if st in seen:
                r = seen[st]
                mr = sorted(m for m in MARK if rep[m] == r); md = sorted(m for m in MARK if rep[m] == D)
                if mr and md:
                    events.append((n + 1, mr, md))
                for m in md:
                    rep[m] = r
                if 1 in groups[r] or 1 in groups[D]:
                    last_merge1 = n + 1
                groups[r] |= groups.pop(D)
            else:
                seen[st] = D; alive[D] = st
        n += 1
        if n <= NABS:
            k1911 = alive[rep[1911]][0]
            if k1911 > max1911:
                max1911 = k1911
        if merge2 is None and n > NJ and len(alive) == 2:
            merge2 = n
        if n in cps:
            levels[n] = {m: alive[rep[m]][0] for m in (1, 1911, 2895)}
            print(f'     step {n}: levels of the chains of D = 1, 1911, 2895: {levels[n]}  (live groups {len(alive)}) [{time.time() - t0:.0f}s]', flush=True)
        if n == NJ:
            g1 = groups[rep[1]]
            others = sorted((min(g), max(g), len(g)) for r, g in groups.items() if r != rep[1])
            print(f'(J) at step {NJ}: group of D = 1 has {len(g1)} members, those above 1910: {sorted(m for m in g1 if m > 1910)}; '
                  f'its last merge at step {last_merge1}; other clusters (min D, max D, size): {others}')
            check(not absorbed and g1 == set(range(1, 1911)) | EXTRA and last_merge1 == 2741933 and len(others) == 2
                  and sum(o[2] for o in others) == DMAX - 1916,
                  'the group of D = 1 is exactly {1..1910} u {1929..1932, 1935, 1936} (1916 members, last merge at 2,741,933); '
                  'the other 2084 chains form two clusters; no absorption yet')
    print('(K) to step 19,000,766')
    print(f'     merge events of marked chains (step, surviving, joining): {events}')
    r1 = rep[1]; r2 = rep[1911]
    g2 = groups[r2]
    print(f'     the two lower clusters merge at step {merge2}: level-2 cluster of {len(g2)} members, D in [{min(g2)}, {max(g2)}]')
    print(f'     D = 1 group absorbed at Terras time {absorbed.get(r1)}; level-2 cluster absorbed: {absorbed.get(r2)}; '
          f'maximum level of the D = 1911 chain over [0, 19000765]: {max1911}; its level at 19000765: {levels[NABS][1911]}')
    check(merge2 == 4474989 and rep[2895] == r2 and rep[4000] == r2 and rep[2902] == r2 and len(g2) == DMAX - 1916
          and absorbed.get(r1) == NABS and len(groups[r1]) == 1916 and r2 not in absorbed and len(alive) == 2,
          'the two lower clusters merge at 4,474,989 (2084 members: level 2); the D = 1 group (1916 members) is absorbed at '
          'Terras time 19,000,765, certifying M_99708993677 ~> M_(99708993677 - D) for D in {1..1910, 1929..1932, 1935, 1936} '
          '(deepest: M_99708991741); the level-2 cluster is not absorbed by 19,000,766')

    print('(N) odd-step counts at Terras time 19,000,765')
    cnt = {}
    for D in (0, 1, 256):
        _, a, _ = F.par(F.xres(E0 - D, NABS + 8), NABS)
        cnt[D] = a
    print(f'     odd steps of x_E0, y_1, y_256 over 19,000,765 Terras steps: {cnt[0]}, {cnt[1]}, {cnt[256]}')
    check(cnt[1] - cnt[0] == 1 and cnt[256] - cnt[0] == 256, 'the odd-step counts differ by exactly D (relation level 0 at absorption)')
    print(f'\n{"ALL CHECKS PASSED" if not FAIL else "FAILURES: " + str(FAIL)}  ({time.time() - t0:.0f}s)')
