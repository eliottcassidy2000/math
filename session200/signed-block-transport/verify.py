"""Exact short-block slope transport and profile-fiber normalization."""
from pathlib import Path
from math import gcd
import hashlib
import json
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE / 'dependencies'))
from source_helpers import apply, graph, counts, g, T

def need(ok, message):
    if not ok:
        raise RuntimeError(message)

def parameters(a, b):
    need(b % 2 == 1 and gcd(a, b) == 1 and b < 2 * a < 2 * b, 'admissible input')
    j = 1
    while 2 * (b - (1 << j) * (b - a)) > b:
        j += 1
    t = b - (1 << j) * (b - a)
    ell = 1
    while 2 * (1 << ell) * t < b:
        ell += 1
    aa = (1 << ell) * t
    c = (1 << j) * (b - a) + a - 2 * b
    return ell, j, aa, c

def word_from_U(n):
    need(n >= 2, 'nontrivial canonical source')
    out = []
    while n != 2:
        if n % 2:
            out.append('Z')
            n = (n + 1) // 2
        else:
            out.append('O')
            n //= 2
    return ''.join(reversed(out))

def main():
    basis = [tuple(int(i == j) for i in range(4)) for j in range(4)]
    slopes, raw_checks = [], 0
    for b in range(3, 102, 2):
        for a in range(b // 2 + 1, b):
            if gcd(a, b) > 1:
                continue
            ell, j, aa, c = parameters(a, b)
            W = 'O' * ell + 'Z' * j
            need(aa % 2 == 0 and gcd(aa, b) == 1 and b < 2 * aa < 2 * b, 'effective numerator')
            need(ell >= 1 and j >= 1 and ell + j <= 2 * (b - 1).bit_length(), 'vertex-cost bound')
            need(c % 2 == a % 2, 'response parity shift')
            for s in basis + [(17, 4, 29, 3), (-9, 5, -3, -2)]:
                D, E, U, d = s
                q = apply(s, W)
                need(q == (D + ((1 << j) - 1) * (1 << ell) * U - ((1 << j) - 2) * d,
                           d, (1 << (ell + j)) * U - ((1 << j) - 1) * d, d), 'full block profile')
                need(a * q[2] - b * q[0] == aa * U - b * D + c * d, 'full signed identity')
                raw_checks += 1
            slopes.append(dict(a=a, b=b, ell=ell, j=j, effective_a=aa, shift=c, word=W))
    literal = []
    for a, b in ((3, 5), (5, 7), (5, 9), (7, 9)):
        ell, j, aa, c = parameters(a, b)
        W = 'Z' + 'O' * ell + 'Z' * j
        rows, root = graph(W)
        prof, start, cyc = counts(rows, root)
        need(prof == apply(T, W) and prof[1] == prof[3] == start == 1 and cyc == 1, 'literal block source')
        need(a * prof[2] - b * prof[0] == aa * 3 - b * 3 + c, 'literal signed transfer')
        literal.append(dict(a=a, b=b, word=W, profile=prof, order=len(rows)))
    # Exhaustive declared finite universe, separate from the all-order proof.
    cap, checked, positive = 100001, 0, 0
    for n in range(3, cap + 1, 2):
        left = 5 * n - 9 * g(n) == 24
        q = (n + 1) // 16
        right = (n + 1) % 16 == 0 and q >= 3 and q % 2 == 1 and 8 * q - 9 * g(q) == 29
        need(left == right, 'bounded exact fiber iff')
        positive += left
        checked += 1
    for q in range(3, 10002, 2):
        need(g(16 * q - 1) == g(q) + 8 * q, 'independent g block recurrence')
    q, n = 675385, 10806159
    qword, nword = word_from_U(q), word_from_U(n)
    need(g(q) == 600339 and g(n) == 6003419, 'large g witness')
    need(8 * q - 9 * g(q) == 29 and 5 * n - 9 * g(n) == 24, 'both exact fibers')
    need(n == 16 * q - 1 and nword == qword + 'OOOZ', 'physical canonical suffix')
    need(apply(T, qword) == (600339, 1, 675385, 1), 'physical predecessor profile')
    need(apply(T, nword) == (6003419, 1, 10806159, 1), 'physical successor profile')
    jj = (q - 7) // 18
    need(jj == 37521 and q == 18 * jj + 7 and g(q) == 16 * jj + 3, 'target index')
    need(79232 * jj + 30133 == 2972894005, 'common positive target')
    need(g(7) == 5 and g(7) != 3, 'small index nonphysical control')
    witness = dict(q_word=qword, q_profile=(600339, 1, 675385, 1),
                   successor_word=nword, successor_profile=(6003419, 1, 10806159, 1),
                   index=jj, common_C2_target=2972894005)
    out = dict(status='PASS', max_odd_denominator=101, slope_count=len(slopes),
               slopes=slopes, exact_basis_and_signed_controls=raw_checks,
               literal_sources=literal, bounded_fiber_domain=dict(odd_n_from=3, odd_n_through=cap, checked=checked, positive=positive),
               large_exact_witness=witness,
               dependency_hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                                  for p in sorted((HERE / 'dependencies').glob('*')) if p.is_file()},
               scope='Fixed short physical signed-slope transport and exact profile-level fiber normal form; no universal positive-target coverage.')
    content = json.dumps(out, sort_keys=True, indent=2) + '\n'
    (HERE / 'VALIDATION.json').write_text(content)
    print(content, end='')

if __name__ == '__main__':
    main()
