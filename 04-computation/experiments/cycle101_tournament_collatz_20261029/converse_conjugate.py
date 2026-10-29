"""THM-4848 (v): the converse labelling conjugates the block determinant.  M_k(T^op) = M_k(T)^T gives
B_j(T^op) = B_(-j)(T)^T, so det B_1(T(NQR)) = conj(det B_1(T(QR))) in Q(sqrt(-q)).  Checked for P7 (k = 2, 3) and
P11 (k = 2): with beta = a + b*eta, conj(beta) = (a - b) - b*eta, since eta + conj(eta) = -1."""
import sys, itertools
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from fast_runs import m2_fast
from runs import circulant, run_matrix, gp_batch


def det_b1(q, S, k):
    A = circulant(q, set(S))
    if k == 2:
        T, _ = m2_fast(A)
        tuples = [(a, b) for a in range(q) for b in range(q) if a != b]
        M = {(s, t): int(T[s[0], s[1], t[0], t[1]]) for s in tuples for t in tuples}
    else:
        tuples, Ml = run_matrix(A, k, list(itertools.permutations(range(q), k)))
        M = {(s, t): Ml[i][j] for i, s in enumerate(tuples) for j, t in enumerate(tuples)}
    shift = lambda t, s: tuple((x + s) % q for x in t)
    reps, seen = [], set()
    for t in tuples:
        if t not in seen:
            seen.update(shift(t, s) for s in range(q))
            reps.append(t)
    entry = lambda r, rp: '(' + ('+'.join('%d*z^%d' % (M[(r, shift(rp, s))], s) for s in range(q) if M[(r, shift(rp, s))]) or '0') + ')'
    B1 = 'Mat([' + ';'.join(','.join(entry(r, rp) for rp in reps) for r in reps) + '])'
    QR = sorted({(x * x) % q for x in range(1, q)})
    out = gp_batch(['my(d=lift(Mod(matdet(%s), polcyclo(%d, z)))); [polcoeff(d, 0), polcoeff(d, 1), d - polcoeff(d,0) - polcoeff(d,1)*(%s)]'
                    % (B1, q, '+'.join('z^%d' % j for j in QR))], timeout=7200)[0]
    a, b, rest = out.strip('[]').split(',', 2)
    assert rest.strip() == '0'
    return int(a), int(b)


for q, k in ((7, 2), (7, 3), (11, 2)):
    QR = {(x * x) % q for x in range(1, q)}
    NQR = set(range(1, q)) - QR
    a, b = det_b1(q, QR, k)
    c, d = det_b1(q, NQR, k)
    print('P%d k=%d: beta(QR) = %d + %d eta ; beta(NQR) = %d + %d eta ; conj(beta(QR)) = %d + %d eta ; equal: %s'
          % (q, k, a, b, c, d, a - b, -b, (c, d) == (a - b, -b)))
