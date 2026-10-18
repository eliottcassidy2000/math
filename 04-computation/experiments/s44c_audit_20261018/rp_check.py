#!/usr/bin/env python3
"""Audit check (e)(vii): the abstracted final step of the RP argument and the copyright constraint.

1. Lemma: integers M_n >= 0 with (i) M_0 = M_2 = 0, (ii) 2 M_(a+b) <= M_(2a) + M_(2b), (iii) M_(2a) in mu Z,
   (iv) M_n <= sigma n.  If mu > 2 sigma then M = 0.
   - the sharpness family M_n = m max(0, n-2) (mu = 2m) satisfies (i)-(iv) whenever m <= sigma, i.e. mu <= 2 sigma.
2. Copyright: longest run of consecutive identical words between section 3 of the merge note and the extracted
   text of the source PDF (untrusted data, read only, path given on the command line).
"""
import sys, re, os, json, itertools

HERE = os.path.dirname(os.path.abspath(__file__))


def satisfies(M, mu, sigma, K, quarter=None):
    if M[0] != 0 or M[2] != 0:
        return False
    for n in range(K + 1):
        if M[n] < 0 or M[n] > sigma * n:
            return False
        if quarter and M[n] % quarter:
            return False
    for a in range(K // 2 + 1):
        if M[2 * a] % mu:
            return False
    for a in range(K // 2 + 1):
        for b in range(K // 2 + 1):
            if a + b <= K and 2 * M[a + b] > M[2 * a] + M[2 * b]:
                return False
    return True


def brute(mu, sigma, K, quarter=None):
    """Search all sequences with M_n in [0, sigma n] (step 'quarter' if given). Returns count of nonzero solutions."""
    step = quarter or 1
    ranges = [range(0, sigma * n + 1, step) for n in range(K + 1)]
    ranges[0] = range(0, 1)
    ranges[2] = range(0, 1)
    cnt = 0
    example = None
    # DFS with pruning on (ii) and (iii) for prefixes
    M = [0] * (K + 1)

    def ok_prefix(n):
        if n % 2 == 0 and M[n] % mu:
            return False
        # constraints involving indices <= n: 2 M_(a+b) <= M_2a + M_2b with 2a, 2b, a+b <= n
        for a in range(n // 2 + 1):
            for b in range(a, n // 2 + 1):
                if 2 * b <= n and 2 * M[a + b] > M[2 * a] + M[2 * b]:
                    return False
        return True

    def rec(n):
        nonlocal cnt, example
        if n > K:
            if any(M):
                cnt += 1
                if example is None:
                    example = list(M)
            return
        for x in ranges[n]:
            M[n] = x
            if ok_prefix(n):
                rec(n + 1)
        M[n] = 0
    rec(0)
    return cnt, example


def longest_common_word_run(a, b):
    wa = re.findall(r"[A-Za-z0-9]+", a.lower())
    wb = re.findall(r"[A-Za-z0-9]+", b.lower())
    # dynamic programming on words (rolling)
    best, pos = 0, None
    prev = [0] * (len(wb) + 1)
    for i in range(1, len(wa) + 1):
        cur = [0] * (len(wb) + 1)
        for j in range(1, len(wb) + 1):
            if wa[i - 1] == wb[j - 1]:
                cur[j] = prev[j - 1] + 1
                if cur[j] > best:
                    best, pos = cur[j], (i, j)
        prev = cur
    return best, ' '.join(wa[pos[0] - best:pos[0]]) if pos else ''


if __name__ == '__main__':
    out = {}
    # (A brute-force search over truncated sequences is not informative: for odd n > K/2 the truncated system has no
    #  constraint 2 M_n <= M_0 + M_(2n), so nonzero tails always exist.  The lemma is checked by hand in the report.)
    fam = []
    for m in range(1, 9):
        for sigma in range(1, 10):
            M = [m * max(0, n - 2) for n in range(0, 61)]
            fam.append({'m': m, 'sigma': sigma, 'mu': 2 * m, 'satisfies': satisfies(M, 2 * m, sigma, 60),
                        'mu_le_2sigma': 2 * m <= 2 * sigma})
    agree = all(f['satisfies'] == f['mu_le_2sigma'] for f in fam)
    out['sharpness_family_satisfies_iff_mu_le_2sigma'] = agree
    print('family M_n = m max(0, n-2): satisfies (i)-(iv) on n <= 60 exactly when mu = 2m <= 2 sigma:', agree)
    M = [max(0, 4 * n - 8) for n in range(0, 201)]
    out['note_example_mu8_sigma4'] = satisfies(M, 8, 4, 200, quarter=4)
    print('note example M_n = max(0, 4n-8), mu = 8, sigma = 4 (also M_n in 4Z, M_1 = 0):', out['note_example_mu8_sigma4'])
    if len(sys.argv) > 1:
        src = open(sys.argv[1], encoding='latin-1').read()
        note = open(os.path.join(HERE, '..', '..', '..', '05-knowledge', 'results',
                                 'pythagorean_form_conic_graphs_rp_margin_20261018.md'), encoding='utf-8').read()
        sec3 = note[note.index('## 3.'):note.index('## 4.')]
        best, run = longest_common_word_run(sec3, src)
        out['copyright_longest_common_word_run'] = best
        out['copyright_run_text'] = run
        print(f'longest run of identical consecutive words, note section 3 vs source text: {best} words: "{run}"')
    json.dump(out, open(os.path.join(HERE, 'rp_check.json'), 'w'), indent=1)
