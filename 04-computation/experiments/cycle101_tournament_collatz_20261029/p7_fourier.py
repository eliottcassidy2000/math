"""Paley P7 three-run observer under its automorphism group Z_7 x| C_3: where the primes of its Smith form come from.

THREE_RUN_GLUE (R9, R12) gives the Smith form of M_3(P7) (210 x 210): order 2^80 3^18 5^12 11^6 17^4 79^3 5861^6
392809^3.  Translation t: i -> i+1 is an automorphism, so M_3 commutes with it and splits over the characters of Z_7
into blocks B_j (30 x 30, over Z[zeta_7]); det M_3 = det B_0 * N(det B_1), and the multiplier i -> 2i (also an
automorphism) makes B_1, B_2, B_4 conjugate, so det B_1 lies in Q(sqrt(-7)) up to a unit and every prime of
N(det B_1) = N_{Q(sqrt-7)}(det B_1)^3 appears to a power divisible by 3 and splits or ramifies in Q(sqrt(-7))
(p = 7 or p = 1, 2, 4 mod 7).  This script computes the determinants exactly and checks the claims.
"""
import sys, itertools, json
sys.path.insert(0, '../growth_gf_runs_fermat_20261029')
from runs import run_matrix, circulant, gp_batch, gp_matrix, parse_vec
from sympy import factorint

A = circulant(7, {1, 2, 4})
tuples, M = run_matrix(A, 3, list(itertools.permutations(range(7), 3)))
idx = {t: i for i, t in enumerate(tuples)}
# orbits of ordered triples under translation (free: 7 prime, triples of distinct points)
reps = []
seen = set()
for t in tuples:
    if t in seen:
        continue
    orb = [tuple((x + s) % 7 for x in t) for s in range(7)]
    seen.update(orb)
    reps.append(t)
assert len(reps) == 30
# block B_j[r, r'] = sum_s M[r, r'+s] zeta^(j s); represent as a gp expression in Mod(z, polcyclo(7))
def block_entry(r, rp, j):
    terms = []
    for s in range(7):
        val = M[idx[r]][idx[tuple((x + s) % 7 for x in rp)]]
        if val:
            terms.append('%d*z^%d' % (val, (j * s) % 7))
    return '(' + ('+'.join(terms) if terms else '0') + ')'
cmds = []
for j in (0, 1):
    rows = [','.join(block_entry(r, rp, j) for rp in reps) for r in reps]
    mat = 'Mat([' + ';'.join(rows) + '])'
    if j == 0:
        cmds.append('matdet(subst(%s, z, 1))' % mat)
    else:
        cmds.append('my(B=%s, d=lift(Mod(matdet(B), polcyclo(7, z)))); [d, norm(Mod(d, polcyclo(7,z)))]' % mat)
cmds.append('matdet(%s)' % gp_matrix(M))
out = gp_batch(cmds)
det0 = int(out[0])
d1, n1 = out[1].strip('[]').rsplit(',', 1)
n1 = int(n1)
detM = int(out[2])
print('det B_0 (trivial character, 30 x 30) =', det0, '=', factorint(abs(det0)))
print('det B_1 =', d1.strip())
print('N(det B_1) =', n1, '=', factorint(abs(n1)))
print('det M_3(P7) =', detM, '=', factorint(abs(detM)))
print('det B_0 * N(det B_1) == det M_3 ?', det0 * n1 == detM)
for p, e in factorint(abs(n1)).items():
    print('   prime %d: exponent %d in N(det B_1), divisible by 3: %s, residue mod 7: %d (QR mod 7: %s)'
          % (p, e, e % 3 == 0, p % 7, p % 7 in (1, 2, 4)))
json.dump({'det0': det0, 'norm1': n1, 'det': detM}, open('out/p7_fourier.json', 'w'))
