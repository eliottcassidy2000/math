"""audit D2 (THM-5056), items 3-4: symbolic sweep of the induction for every odd k < 2^LMAX.
Gamma_odd(c) is computed as the closure of the first-PH orders of the seeds under +2 / +3 (Lemma O2); each k is
assigned to exactly one case and the case hypotheses are checked; the recursion in case 2 must bottom out in the base
(odd k < 512, minus {1, 5, 37}) or in an offset / special case."""
import sys
LMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 24
seed_first = {5: [4 + 3], 1: [6 + 3], 37: [7 + 3, 8 + 3]}          # |seed| + 3 after PH(w, 2)
gam = {}
for c, starts in seed_first.items():
    S = set()
    stack = list(starts)
    while stack:
        n = stack.pop()
        if n in S or n > LMAX + 5: continue
        S.add(n); stack += [n + 2, n + 3]
    gam[c] = S
for c in (5, 1, 37):
    print(f'Gamma_odd({c}) up to {LMAX}:', sorted(x for x in gam[c] if x <= LMAX))
gaps = [(c, L) for c in (1, 5, 37) for L in range(10, LMAX + 1) if L not in gam[c]]
print('(x, L) gaps with L >= 10:', gaps)
assert 513 == 3 * 171 and (171).bit_length() + 2 == 10
ok = {}
cases = {'base': 0, 'c1': 0, 'c2': 0, 'c3': 0, 'c3_513': 0}
def good(k):
    if k in ok: return ok[k]
    L = k.bit_length()
    if k < 512:
        r = k not in (1, 5, 37); cases['base'] += r
    elif k >= 3 * 2 ** (L - 2):
        t = k - 3 * 2 ** (L - 2); r = (t % 2 == 1 and 1 <= t < 2 ** (L - 2) and t.bit_length() <= L - 2); cases['c1'] += 1
    else:
        x = k - 2 ** (L - 1)
        assert x % 2 == 1 and 1 <= x < 2 ** (L - 2)
        if x not in (1, 5, 37):
            r = good(x) and x.bit_length() <= L - 2; cases['c2'] += 1
        elif L in gam[x]:
            r = True; cases['c3'] += 1
        elif (x, L) == (1, 10):
            r = (k == 513); cases['c3_513'] += 1
        else:
            r = False
    ok[k] = r; return r
bad = [k for k in range(1, 2 ** LMAX, 2) if not good(k)]
print('odd k < 2^%d failing the case analysis:' % LMAX, bad[:10], '(count %d)' % len(bad))
print('case counts:', cases)
