"""Auditor: deficit Delta = h - b (HYP-9369 section 5) with independent code.
h = H/|Aut|; b = #blue/|Aut|; a Hamiltonian word v_1..v_n (v_i -> v_(i+1)) is blue iff the positional reversal
psi: v_i -> v_(n+1-i) is an anti-automorphism: A[x][y] == A[psi y][psi x] for all x != y.
|Aut| and self-converseness by backtracking (tourn.all_automorphisms).  Also checks the counting bound
#blue <= #(anti-involutions) * 2^floor(n/2) floor(n/2)!  and  b <= 2^floor(n/2) floor(n/2)!.
Usage: python deficit_check.py > deficit_check.out"""
import os, sys, math
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import tourn

HERE = os.path.dirname(os.path.abspath(__file__))
UP = os.path.join(HERE, '..')


def ham_words(A):
    n = A.shape[0]
    out = [list(np.nonzero(A[v])[0]) for v in range(n)]
    words = []
    w = []

    def dfs(v, used):
        w.append(v)
        if len(w) == n:
            words.append(tuple(w))
        else:
            for x in out[v]:
                if not (used >> x) & 1:
                    dfs(x, used | (1 << x))
        w.pop()
    for v in range(n):
        dfs(v, 1 << v)
    return words


def is_blue(A, word):
    n = len(word)
    psi = [0] * n
    for i in range(n):
        psi[word[i]] = word[n - 1 - i]
    # psi anti-automorphism: A[x][y] == A[psi y][psi x] for all x != y
    P = np.array(psi)
    return bool(np.all(A == A[np.ix_(P, P)].T)), psi


def anti_involutions(A):
    """count anti-automorphisms p with p o p = id, by backtracking over anti-automorphisms."""
    n = A.shape[0]
    s = A.sum(1); ts = A.T.sum(1)
    cnt_all = 0; cnt_inv = 0
    p = [-1] * n; used = [False] * n

    def bt(i):
        nonlocal cnt_all, cnt_inv
        if i == n:
            cnt_all += 1
            if all(p[p[x]] == x for x in range(n)):
                cnt_inv += 1
            return
        for y in range(n):
            if used[y] or ts[y] != s[i]:
                continue
            if all(A[i, j] == A[p[j], y] for j in range(i)):      # A[i][j] == A[p j][p i]
                p[i] = y; used[y] = True
                bt(i + 1)
                used[y] = False; p[i] = -1
    bt(0)
    return cnt_all, cnt_inv


def deficit(A):
    n = A.shape[0]
    words = ham_words(A)
    blue = [wd for wd in words if is_blue(A, wd)[0]]
    aut = tourn.all_automorphisms(A)
    anti, ainv = anti_involutions(A)
    H = len(words); nb = len(blue)
    assert H % aut == 0 and nb % aut == 0
    psis = set(tuple(is_blue(A, wd)[1]) for wd in blue)
    bound = 2 ** (n // 2) * math.factorial(n // 2)
    return dict(H=H, aut=aut, sc=anti > 0, anti=anti, antiinv=ainv, blue=nb, h=H // aut, b=nb // aut,
                D=H // aut - nb // aut, distinct_psi=len(psis), bound_ok=(nb <= ainv * bound and nb // aut <= bound))


def fmt(name, d):
    return (f"{name}: H={d['H']} |Aut|={d['aut']} self-converse={d['sc']} anti-auts={d['anti']} anti-involutions={d['antiinv']} "
            f"#blue={d['blue']} (distinct psi {d['distinct_psi']}) h={d['h']} b={d['b']} Delta={d['D']} bound_ok={d['bound_ok']}")


def main():
    K2 = np.array([[0, 1], [0, 0]])
    C3 = tourn.circulant(3, [1])
    print(fmt('C_5', deficit(tourn.circulant(5, [1, 2]))))
    print(fmt('C_3[K_2]', deficit(tourn.lex(C3, K2))))
    print(fmt('P_7', deficit(tourn.circulant(7, [1, 2, 4]))))
    print(fmt('C_7(1,2,3)', deficit(tourn.circulant(7, [1, 2, 3]))))
    print(fmt('C_3[C_3]', deficit(tourn.lex(C3, C3))))
    # all classes at n = 6 (file), regular 7 (3 classes via my canonical list in zeros_small), n = 8 (file vs deficit_ar8.txt)
    for n in (6,):
        print(f'-- almost-regular {n} (all classes, ar{n}_classes.txt)')
        for l in open(os.path.join(UP, f'ar{n}_classes.txt')):
            A = tourn.from_code(l.split()[0]); d = deficit(A)
            print('  ' + fmt(l.split()[0], d))
    print('-- regular 7: the three classes')
    for code in ('000111010010010100010', '000111000110001000000', '000111000110100000010'):
        print('  ' + fmt(code, deficit(tourn.from_code(code))))
    print('-- almost-regular 8: every row of deficit_ar8.txt (current version) recomputed (file columns: canon H aut sc blue h b Delta)')
    rows = [l.split() for l in open(os.path.join(UP, 'deficit_ar8.txt')) if l.strip() and not l.startswith('#')]
    mism = 0; sc = 0; minD = None; maxb = 0
    for r in rows:
        A = tourn.from_code(r[0]); d = deficit(A)
        mine = [d['H'], d['aut'], int(d['sc']), d['blue'], d['h'], d['b'], d['D']]
        theirs = list(map(int, r[1:8]))
        if mine != theirs:
            mism += 1; print('  MISMATCH', r, mine)
        if d['sc']:
            sc += 1
            minD = d['D'] if minD is None else min(minD, d['D'])
            maxb = max(maxb, d['b'])
        assert d['bound_ok']
    print(f'  rows {len(rows)}, mismatches {mism}, self-converse {sc}, min Delta (self-converse) {minD}, max b {maxb}')
    print('-- regular 9: all 15 rows of deficit_reg9.txt')
    rows = [l.split() for l in open(os.path.join(UP, 'deficit_reg9.txt')) if l.strip() and not l.startswith('#')]
    mism = 0; minD = None; maxb = 0
    for r in rows:
        A = tourn.from_code(r[0]); d = deficit(A)
        mine = [d['H'], d['aut'], int(d['sc']), d['blue'], d['h'], d['b'], d['D']]
        if mine != list(map(int, r[1:8])):
            mism += 1; print('  MISMATCH', r, mine)
        if d['sc']:
            minD = d['D'] if minD is None else min(minD, d['D']); maxb = max(maxb, d['b'])
    print(f'  rows {len(rows)}, mismatches {mism}, min Delta {minD}, max b {maxb}')
    # file-level summaries for n = 10, 11 and spot checks of the extreme rows
    for fn, n in (('deficit_ar10.txt', 10), ('deficit_reg11.txt', 11)):
        rows = [l.split() for l in open(os.path.join(UP, fn)) if l.strip() and not l.startswith('#')]
        scr = [r for r in rows if r[3] == '1']
        Ds = [int(r[7]) for r in scr]; bs = [int(r[6]) for r in scr]
        print(f'-- {fn}: rows {len(rows)}, self-converse {len(scr)}, min Delta {min(Ds)}, max b {max(bs)}, '
              f'nonfree rows {sum(1 for r in rows if "NONFREE" in r)}, Delta = 6 rows {Ds.count(6)}')
        picks = [scr[int(np.argmin(Ds))], scr[int(np.argmax(bs))]]
        if n == 10:
            picks += [r for r in scr if int(r[2]) > 1][:2]
        for r in picks:
            d = deficit(tourn.from_code(r[0]))
            mine = [d['H'], d['aut'], int(d['sc']), d['blue'], d['h'], d['b'], d['D']]
            print(f'   spot check {r[0]}: file {r[1:8]} mine {mine} agree {mine == list(map(int, r[1:8]))}')


if __name__ == '__main__':
    main()
