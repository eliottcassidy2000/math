"""Independent finite-difference audit of the sibling full-context loop proof."""
from math import comb, factorial, gcd, lcm
import json
from check import count_paths, require


def factors(q):
    out, p = [], 2
    while p*p <= q:
        a = 0
        while q % p == 0:
            q //= p
            a += 1
        if a:
            out.append((p, a))
        p += 1
    if q > 1:
        out.append((q, 1))
    return out


def claimed(q):
    fs = factors(q)
    return max(a for p, a in fs), lcm(*(p**(a-1)*(p-1) if p != 2 else (1 if a <= 2 else 2**max(1, a-3)) for p, a in fs))


def actual_orbit(q):
    K = next(k for k in range(1, q+2) if factorial(k-1)*factorial(k) % q == 0)
    coeffs = [[factorial(k-1)*(-1)**(k-j)*comb(k, j) % q for j in range(1, k+1)] for k in range(1, K)]
    powers = [j % q for j in range(1, K)]
    seen, n = {}, 1
    while True:
        state = tuple(sum(a*powers[j] for j, a in enumerate(row)) % q for row in coeffs)
        if state in seen:
            return seen[state], n-seen[state], len(seen), K
        seen[state] = n
        powers = [v*(j+1) % q for j, v in enumerate(powers)]
        n += 1
        require(n < 30000, 'unexpected audit horizon')


def separator(n, m):
    size = n+2*m
    rows = [0]*size
    for i in range(size):
        for j in range(i+1, size):
            if j < n or i >= n or (j-n) % 2 == 0:
                rows[i] |= 1 << j
            else:
                rows[j] |= 1 << i
    return rows


def main():
    universe = sorted(set(range(2, 97)) | {128, 256, 512, 1024, 121, 343, 625, 729, 196, 5929})
    rows = []
    for q in universe:
        actual = actual_orbit(q)
        tau, period = claimed(q)
        require(actual[:3] == (tau, period, tau-1+period), f'loop theorem failed modulo {q}: {actual}')
        rows.append([q, *actual])
    controls = []
    for n in range(1, 6):
        d1, d2, d3 = 1, 2**n-2, 2*3**n-6*2**n+6
        values = [count_paths(separator(n, m))[0][-1] for m in range(3)]
        require(values == [d1, 3*d1+d2, 9*d1+10*d2+d3], 'literal triangular separator probe failed')
        controls.append([n, *values])
    print(json.dumps({'method': 'Direct finite-difference power sums; first repeated full profile found without an injected tail or period.', 'moduli': len(universe), 'orbits_q_tau_period_states_cutoff': rows, 'literal_separator_controls': controls, 'verdict': 'PASS'}, sort_keys=True, indent=2))


if __name__ == '__main__':
    main()
