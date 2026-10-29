"""P19 and P23: verify that the nontrivial-character determinant of the two-run observer lies in Q(sqrt(-q)),
so that its absolute norm is an exact ((q-1)/2)-th power of a quadratic norm; partially factor that quadratic norm
and the trivial block.  (Companion to paley_fourier.py, which fully factored P7 and P11.)"""
import sys, json
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from fast_runs import m2_fast
from runs import circulant, gp_batch
res = {}
for q in [int(x) for x in sys.argv[1:]] or [19]:
    QR = sorted({(x * x) % q for x in range(1, q)})
    g = 2 if q == 19 else 5          # a generator of Z_q^*
    h2 = pow(g, 2, q)                # a quadratic residue of order (q-1)/2
    A = circulant(q, set(QR))
    T, _ = m2_fast(A)
    pairs = [(s1, s2) for s1 in range(q) for s2 in range(q) if s1 != s2]
    reps, seen = [], set()
    for t in pairs:
        if t in seen:
            continue
        seen.update(((t[0] + s) % q, (t[1] + s) % q) for s in range(q))
        reps.append(t)
    def entry(r, rp, j):
        terms = []
        for s in range(q):
            val = int(T[r[0], r[1], (rp[0] + s) % q, (rp[1] + s) % q])
            if val:
                terms.append('%d*z^%d' % (val, (j * s) % q))
        return '(' + ('+'.join(terms) if terms else '0') + ')'
    B1 = 'Mat([' + ';'.join(','.join(entry(r, rp, 1) for rp in reps) for r in reps) + '])'
    B0 = 'Mat([' + ';'.join(','.join(entry(r, rp, 0) for rp in reps) for r in reps) + '])'
    P = 'polcyclo(%d, z)' % q
    cmd = ('my(d=lift(Mod(matdet(%s), %s)), s=lift(Mod(subst(d, z, z^%d), %s)), c=lift(Mod(subst(d, z, z^%d), %s)),'
           ' qn=lift(Mod(d*c, %s)), n=norm(Mod(d, %s)), d0=matdet(subst(%s, z, 1)));'
           ' [s==d, type(qn), qn, n==qn^%d, factor(qn, 10^6), factor(d0, 10^6)]'
           % (B1, P, h2, P, q - 1, P, P, P, B0, (q - 1) // 2))
    out = gp_batch([cmd], timeout=7200)[0]
    print('P_%d (h = %d of order %d):' % (q, h2, (q - 1) // 2), out[:1500])
    res[q] = out
json.dump(res, open('out/paley_fourier2.json', 'w'))
