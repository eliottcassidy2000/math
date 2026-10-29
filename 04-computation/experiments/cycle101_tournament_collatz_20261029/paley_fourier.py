"""Fourier-Galois structure of the two-run observer of Paley tournaments P_q (q = 3 mod 4 prime).

Translations are automorphisms, so M_2(P_q) splits over the characters of Z_q into blocks B_j (q-1 x q-1 over
Z[zeta_q]); the quadratic residues act by automorphisms (multipliers), so B_{hj} ~ B_j for h a QR, and
det B_1 is fixed by Gal(Q(zeta_q)/Q(sqrt(-q))).  Prediction: det M_2 = det B_0 * N_{Q(sqrt-q)/Q}(det B_1)^((q-1)/2),
so every prime outside det B_0 has exponent divisible by (q-1)/2, and primes with odd exponent in the
quadratic norm split or ramify in Q(sqrt(-q)).  Run for q = 7, 11 (default; q from argv).  P19 is checked by
paley_fourier2.py, which avoids sympy's factorint on the 300-digit norm."""
import sys, itertools, json
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from fast_runs import m2_fast
from runs import circulant, gp_batch
from sympy import factorint, legendre_symbol

res = {}
for q in [int(x) for x in sys.argv[1:]] or [7, 11]:
    QR = {(x * x) % q for x in range(1, q)}
    A = circulant(q, QR)
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
    cmds = []
    for j in (0, 1):
        mat = 'Mat([' + ';'.join(','.join(entry(r, rp, j) for rp in reps) for r in reps) + '])'
        if j == 0:
            cmds.append('matdet(subst(%s, z, 1))' % mat)
        else:
            cmds.append('my(d=lift(Mod(matdet(%s), polcyclo(%d, z)))); [norm(Mod(d, polcyclo(%d,z)))]' % (mat, q, q))
    full = 'Mat([' + ';'.join(','.join(str(int(T[a[0], a[1], b[0], b[1]])) for b in pairs) for a in pairs) + '])'
    cmds.append('matdet(%s)' % full)
    out = gp_batch(cmds, timeout=7200)
    d0 = int(out[0]); n1 = int(out[1].strip('[]')); dM = int(out[2])
    h = (q - 1) // 2
    f1 = factorint(abs(n1)) if n1 else {}
    ok_exp = all(e % h == 0 for e in f1.values())
    print('P_%d: det B_0 = %s; det M_2 = det B_0 * N(det B_1): %s' % (q, factorint(abs(d0)), d0 * n1 == dM))
    print('   N(det B_1) = %s ; all exponents divisible by (q-1)/2 = %d: %s' % (f1, h, ok_exp))
    for p, e in f1.items():
        if p != q:
            quad = e // h
            print('      p=%d exponent %d (quadratic-norm exponent %d), (-q/p) = %d' % (p, e, quad, legendre_symbol(-q % p, p) if p != 2 else 0))
    res[q] = {'det0': d0, 'norm1': n1, 'det': dM}
json.dump(res, open('out/paley_fourier.json', 'w'), default=str)
