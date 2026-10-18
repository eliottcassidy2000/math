#!/usr/bin/env python3
"""From the Paley tournament P_23 to the Leech lattice and the moonshine vertex algebra V^natural (S43).

The Golay code C is generated from the quadratic residues mod 23 (m24_paley_golay_eta_20261018.py). The Leech lattice,
scaled by sqrt(8), is the set of x in Z^24 with
    x == m 1 (mod 2) for some m in {0, 1},  {i : x_i == m + 2 (mod 4)} pattern: (x - m 1)/2 mod 2 in C,  sum x_i == 4m (mod 8)
(Conway-Sloane, SPLAG ch. 4, sec. 11). Minimal vectors (norm 32 in this scaling, norm 4 after dividing by 8):
    shape (+-2^8, 0^16) on octads with an even number of minus signs:    759 * 2^7   =  97152
    shape (-+3, +-1^23) with signs from a Golay codeword (and the rule):   24 * 2^12  =  98304
    shape (+-4^2, 0^22):                                                 C(24,2) * 4  =   1104
We enumerate them directly from the code, check the norm and lattice conditions, and count 196560.
Then the weight-2 space of the Leech lattice vertex algebra V_Lambda has dimension 196560 + 24 + 300 = 196884, and the
Frenkel-Lepowsky-Meurman orbifold V^natural has (V_Lambda^+)_2 + (twisted)_2 = (98280 + 300) + 98304 = 196884,
the first nontrivial coefficient of J(q) = j(q) - 744 (moonshine); checked against q-expansions of j.
Usage: python paley_golay_leech_moonshine_20261018.py
Session opus-2026-10-17-S43.
"""
import os, sys, json, itertools
from math import comb
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from m24_paley_golay_eta_20261018 import golay


def in_leech(x, code):
    m = x[0] % 2
    if any(v % 2 != m for v in x):
        return False
    y = [((v - m) // 2) % 2 for v in x]
    word = sum(b << i for i, b in enumerate(y))
    if word not in code:
        return False
    return sum(x) % 8 == (4 * m) % 8


def main():
    basis, code = golay()
    octads = [c for c in code if bin(c).count('1') == 8]
    count = {}
    # shape 2^8: octad supports, even number of minus signs
    n1 = 0
    for o in octads:
        supp = [i for i in range(24) if (o >> i) & 1]
        for signs in range(256):
            if bin(signs).count('1') % 2:
                continue
            x = [0] * 24
            for j, i in enumerate(supp):
                x[i] = -2 if (signs >> j) & 1 else 2
            assert in_leech(x, code)
            n1 += 1
    count['2^8'] = n1
    # shape (-+3, +-1^23): x = 1-vector with signs flipped on a Golay codeword, then one coordinate -3 (adjusted)
    n2 = 0
    seen = set()
    for c in code:
        for pos in range(24):
            x = [(-1 if (c >> i) & 1 else 1) for i in range(24)]
            x[pos] = 3 * x[pos] * -1      # flip and triple at one place
            if in_leech(x, code):
                t = tuple(x)
                if t not in seen:
                    seen.add(t)
                    n2 += 1
            x2 = [-v for v in x]
            if in_leech(x2, code):
                t = tuple(x2)
                if t not in seen:
                    seen.add(t)
                    n2 += 1
    count['3 1^23'] = n2
    n3 = 0
    for i, j in itertools.combinations(range(24), 2):
        for si in (4, -4):
            for sj in (4, -4):
                x = [0] * 24
                x[i] = si
                x[j] = sj
                assert in_leech(x, code)
                n3 += 1
    count['4^2'] = n3
    total = n1 + n2 + n3
    print('Leech minimal vectors from the Paley-generated Golay code:', count, 'total', total, '(expected 196560)')
    # j-function check: j = E4^3/Delta; J = j - 744 = q^-1 + 196884 q + ...
    N = 4
    def sigma3(n):
        return sum(d ** 3 for d in range(1, n + 1) if n % d == 0)
    E4 = [1] + [240 * sigma3(n) for n in range(1, N + 2)]
    E43 = [0] * (N + 2)
    for a in range(N + 2):
        for b in range(N + 2 - a):
            for c in range(N + 2 - a - b):
                E43[a + b + c] += E4[a] * E4[b] * E4[c]
    # Delta / q = prod (1 - q^n)^24
    D = [1] + [0] * (N + 1)
    for n in range(1, N + 2):
        for _ in range(24):
            for m in range(N + 1, n - 1, -1):
                D[m] -= D[m - n]
    # j * q = E4^3 / (Delta/q): series division
    jq = [0] * (N + 2)
    for n in range(N + 2):
        s = E43[n] - sum(D[k] * jq[n - k] for k in range(1, n + 1))
        jq[n] = s // D[0]
    J = [jq[0], jq[1] - 744] + jq[2:]
    print('j(q) - 744 = q^-1 (%d) + %d + %d q + %d q^2 + ...' % (J[0], J[1], J[2], J[3]))
    vlam2 = total + 24 + 300
    vnat2 = total // 2 + 300 + 2 ** 12 * 24
    print('dim (V_Lambda)_2 = 196560 + 24 + 300 = %d ; dim (V^natural)_2 = 196560/2 + 300 + 2^12 * 24 = %d ; J coefficient %d'
          % (vlam2, vnat2, J[2]))
    here = os.path.dirname(os.path.abspath(__file__))
    json.dump(dict(counts=count, total=total, J=J, dim_Vlambda_2=vlam2, dim_Vnatural_2=vnat2),
              open(os.path.join(here, '..', '..', '05-knowledge', 'results', 'paley_golay_leech_moonshine_20261018.json'), 'w'), indent=1)


if __name__ == '__main__':
    main()
