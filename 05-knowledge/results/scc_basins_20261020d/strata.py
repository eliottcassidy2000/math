"""Independent finite controls for fixed-mass SCC word rationality.

All atoms are genuine tournaments. The alphabet is explicitly restricted;
G3=z5 etc below are not claims about unrestricted strong-mass classes.
"""
from collections import defaultdict
import importlib.util
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def convolve(a, b, limit):
    c = [0] * (limit + 1)
    for i, x in enumerate(a):
        for j, y in enumerate(b[:limit + 1 - i]):
            c[i + j] += x * y
    return c


def factorizations(q):
    if q == 1:
        yield ()
    for m in (3, 5, 9):
        if q % m == 0:
            for suffix in factorizations(q // m):
                yield (m,) + suffix


def omega(q):
    count, p = 0, 2
    while p * p <= q:
        while q % p == 0:
            count += 1
            q //= p
        p += 1
    return count + (q > 1)


def main():
    spec = importlib.util.spec_from_file_location('centre_controls', HERE / 'centers.py')
    centers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(centers)
    atoms = [('1', centers.ONE, 1), ('C3', centers.C3, 1),
             ('R4', centers.R4, 5), ('C5', centers.C5, 3),
             ('Paley7', centers.decode(238754, 7), 9)]
    checked_atoms = []
    for name, adj, mass in atoms:
        data = centers.literal_masses(adj)
        need(data['h'] == mass, f'actual atom mass {name}')
        checked_atoms.append(dict(name=name, n=len(adj), **data))

    limit, maxmass = 30, 81
    # Independent first method: append any actual alphabet atom.
    dp = [defaultdict(int) for _ in range(limit + 1)]
    dp[0][1] = 1
    for n in range(limit + 1):
        for mass, count in list(dp[n].items()):
            for _, adj, atom_mass in atoms:
                if n + len(adj) <= limit and mass * atom_mass <= maxmass:
                    dp[n + len(adj)][mass * atom_mass] += count

    # Second method: ordered nonrigid factorization and independent rigid gaps.
    rigid = [0] * (limit + 1)
    rigid[0] = 1
    for n in range(1, limit + 1):
        rigid[n] = rigid[n - 1] + (rigid[n - 3] if n >= 3 else 0)
    powers = [[1] + [0] * limit]
    for _ in range(6):
        powers.append(convolve(powers[-1], rigid, limit))
    nonrigid_size = {3: 5, 5: 4, 9: 7}
    tables = {}
    for mass in range(1, maxmass + 1, 2):
        coefficients = [0] * (limit + 1)
        factors = list(factorizations(mass))
        for word in factors:
            degree = sum(nonrigid_size[m] for m in word)
            for n in range(degree, limit + 1):
                coefficients[n] += powers[len(word) + 1][n - degree]
        need(coefficients == [dp[n][mass] for n in range(limit + 1)],
             f'rational word formula mass {mass}')
        if any(coefficients):
            tables[mass] = dict(nonrigid_mass_words=factors, coefficients=coefficients)

    threshold_checks = []
    for j in range(1, 9):
        q0 = 3 ** j
        all_smaller = [(q, omega(q)) for q in range(1, q0, 2)]
        need(all(k < j for _, k in all_smaller), 'minimum odd product for j factors')
        beta = q0 * (q0 - 1)
        need(omega(q0) == j, 'attained prime-factor threshold')
        # Equality in the centre equation forces c=d=1.
        possibilities = [(c, q0 * c - beta // q0)
                         for c in range(1, 2 * q0 + 1, 2)
                         if 1 <= q0 * c - beta // q0 <= c]
        need(possibilities == [(1, 1)], 'unique leading centre mass')
        threshold_checks.append(dict(j=j, q=q0, beta=beta,
                                     exact_positive_pole_order=j + 1,
                                     allowed_leading_centres=['empty', '1', 'C3']))
    need([d['beta'] for d in threshold_checks[:4]] == [6, 72, 702, 6480], 'first thresholds')
    receipt = dict(status='PASS', scope='restricted genuine strong alphabet; formal word checks, not unrestricted atom census',
                   checked_atoms=checked_atoms, n_max=limit, mass_max=maxmass,
                   fixed_mass_tables=tables, exact_threshold_checks=threshold_checks)
    (HERE / 'strata.json').write_text(json.dumps(receipt, indent=2, sort_keys=True) + '\n', encoding='utf-8')
    print('PASS actual atom masses: 1,C3,R4,C5,Paley7')
    print('PASS independent word-DP versus rational rigid-gap expansion: sizes0..30, odd masses1..81')
    print('PASS exact pole thresholds j1..8; first beta values6,72,702,6480')
    print('Scope: restricted finite controls; unrestricted rationality and cutoff follow from strata.md proofs.')


if __name__ == '__main__':
    main()
