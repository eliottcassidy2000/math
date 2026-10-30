"""Exact finite controls and symbolic actual-schedule audit. No factoring benchmark."""
from pathlib import Path
from itertools import product
from math import comb, gcd, lcm
import argparse
import hashlib
import importlib.util
import json

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[2]
CHECKS = 0


def need(ok, label):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise RuntimeError(label)


def rectangles(g):
    # Colexicographic order has the closed pair rank b(b-1)/2+a.
    pairs = [(a, b) for b in range(g) for a in range(b)]
    return [(a, b, c, d) for a, b in pairs for c, d in pairs]


def pair_rank(a, b):
    return b * (b - 1) // 2 + a


def positions(g):
    return [((a*g+b)*g+c)*g+d for a, b, c, d in rectangles(g)]


def digits(code, q, length):
    result = []
    for _ in range(length):
        code, digit = divmod(code, q)
        result.append(digit)
    need(code == 0, 'dense code fits its declared coordinate universe')
    return result


def embed(ds, g, q):
    return sum(d*q**p for d, p in zip(ds, positions(g)))


def restricted(code, g, q):
    return tuple((code//q**p) % q for p in positions(g))


def encodings():
    tested = []
    for g, q in [(0, 2), (1, 3), (2, 2), (2, 3), (3, 2), (3, 3), (4, 2), (4, 5)]:
        r = comb(g, 2)**2
        need(len(rectangles(g)) == r, 'rectangle cardinality')
        for j, (a, b, c, d) in enumerate(rectangles(g)):
            need(pair_rank(a, b)*comb(g, 2)+pair_rank(c, d) == j, 'closed rank')
        count = q**r
        # All small codes; a declared deterministic sample for larger spaces.
        codes = list(range(count)) if count <= 20000 else sorted({0, count-1} | {
            int.from_bytes(hashlib.sha256(str(i).encode()).digest(), 'big') % count
            for i in range(512)})
        images = set()
        for code in codes:
            ds = digits(code, q, r)
            full = embed(ds, g, q)
            need(restricted(full, g, q) == tuple(ds), 'restriction is a left inverse')
            images.add(full)
        need(len(images) == len(codes), 'canonical representatives are distinct')
        tested.append({'g': g, 'q': q, 'universe': str(count), 'codes_checked': len(codes),
                       'exhaustive': len(codes) == count})

    # Exhaust every ambient binary code on a 2x2 grid: 16 coordinates, 1 rectangle.
    fibres = [0, 0]
    totals = [0, 0, 0]
    # Signed, three-port response vectors; no absolute-value replacement.
    f = [(2, -3, 5), (-7, 11, -13)]
    for code in range(2**16):
        key = restricted(code, 2, 2)[0]
        fibres[key] += 1
        for j in range(3):
            totals[j] += f[key][j]
    need(fibres == [2**15, 2**15], 'full exhaustive fibres')
    need(totals == [2**15*(a+b) for a, b in zip(*f)], 'all named signed ports preserved')
    need(restricted(0, 2, 2) == restricted(1, 2, 2), 'same active data')
    need((0 % 2) != (1 % 2), 'hostile: an off-rectangle observer blocks erasure')
    need(sum((1, 0)) == sum((0, 1)) and (1, 0) != (0, 1),
         'equal scalar totals do not identify boundary tables')
    need(2*3+1*5 != 3+5, 'nonuniform actual fibres require retained multiplicities')
    return {'coding': tested, 'ambient_binary_codes': 2**16, 'fibres': fibres,
            'signed_response_sum': totals}


def observers():
    # Generic algebra control only: G=2, lambda=4 is NOT the uniform schedule.
    old_scale, new_scale = 3**15*4, 4
    need(old_scale % 3 == 0 and new_scale % 3 == 1, 'producer removes a bad odd factor')
    need((old_scale*0) % 3 == (old_scale*1) % 3, 'old readings cannot be repaired by division')
    need((new_scale*0) % 3 != (new_scale*1) % 3, 'recomputed reduced readings distinguish')
    cases = 0
    for padding, exponent, bank in product([1, 3, 6, 15, 27, 100], [0, 1, 5, 50],
                                          [(3, 5), (4, 9), (9, 25, 49)]):
        modulus = lcm(*bank)
        full = padding*2**exponent
        qo = modulus//gcd(modulus, full)
        qn = modulus//gcd(modulus, 2**exponent)
        need(qn % qo == 0, 'reduced observer refines full observer')
        need(qn//qo == gcd(modulus, full)//gcd(modulus, 2**exponent), 'exact gain')
        for c in range(2*modulus+1):
            need((full*c % modulus == 0) == (c % qo == 0), 'full fibre')
            need(((2**exponent)*c % modulus == 0) == (c % qn == 0), 'reduced fibre')
        cases += 1
    for b in range(25):
        bank = tuple(range(3, 2*b+9, 2))
        need(lcm(*bank) > 2**b, 'inherited dyadic bounded-count cutoff')
    return {'generic_g2_lambda4': {'old_scale': old_scale, 'new_scale': new_scale,
                                 'old_mod3': old_scale % 3, 'new_mod3': new_scale % 3},
            'effective_modulus_cases': cases, 'dyadic_bounds_b': [0, 24]}


def integer_receipt(x):
    raw = x.to_bytes(max(1, (x.bit_length()+7)//8), 'big')
    result = {'bit_length': x.bit_length(), 'sha256_unsigned_big_endian': hashlib.sha256(raw).hexdigest(),
              'low64': x % 2**64}
    if x.bit_length() <= 256:
        result['decimal'] = str(x)
    return result


def umax(g, width):
    return (12*g**3*width+11*g**4*width+24*g**4)//12


def lnext(k, g, width):
    return (k*umax(g, width)+g**4*width)//k**3


def scheduled(depth):
    """Actual uniform depth thresholds, evaluated as integers, without tables."""
    kb = [48]
    for _ in range(depth+1):
        kb.append(9*kb[-1]**16)
    n = kb[depth]**40
    need(n < kb[depth+1]**40, 'chosen n has precisely the declared ellOf depth')
    # sOf is the least s with n <= 2^(4s+2).
    ceil_log = (n-1).bit_length()
    s = max(0, (ceil_log-2+3)//4)
    need(n <= 2**(4*s+2) and (s == 0 or 2**(4*(s-1)+2) < n), 'sOf leastness')
    pa, pl = 2**(s+1), 2**(2*s+1)
    fw = 4*pl+1
    ktop = kb[depth]
    gtop = 3*ktop**8
    clp = 7*(pa//gtop+2)*pa*fw+3*(pa*pa//gtop+2)*fw
    width = (ktop*(10*pa*pa*fw)+gtop*clp)//ktop**3
    rows = []
    moduli = list(range(3, 40, 2))
    full_residues = {m: 1 for m in moduli}
    reduced_residues = {m: 1 for m in moduli}
    padding_residues = {m: 1 for m in moduli}
    for d in reversed(range(depth)):
        k = kb[d]
        g = 3*k**8
        t = (width+2)//3
        q = t+1
        r = comb(g, 2)**2
        a = g**4-r
        need(4*r == g*g*(g-1)**2 and 0 < r < g**4, 'exact branch-exponent ratio')
        v2 = (q & -q).bit_length()-1
        oddpart = q >> v2
        for m in moduli:
            p = pow(q, a, m)
            binary = pow(2, t*r, m)
            full_residues[m] = full_residues[m]*p*binary % m
            reduced_residues[m] = reduced_residues[m]*binary % m
            padding_residues[m] = padding_residues[m]*p % m
            need(gcd(binary, m) == 1, 'all odd moduli remain units for reduced scale')
        rows.append({'d': d, 'K': integer_receipt(k), 'G': integer_receipt(g),
                     'lambda': integer_receipt(width), 'q': integer_receipt(q), 'v2_q': v2,
                     'odd_part_q': integer_receipt(oddpart),
                     'full_branch_exponent_G4': integer_receipt(g**4),
                     'reduced_branch_exponent_R': integer_receipt(r),
                     'padding_exponent': integer_receipt(a),
                     'q_divisors_in_odd_bank': [m for m in moduli if q % m == 0]})
        width = lnext(k, g, width)
    bank = lcm(*moduli)
    need(all(full_residues[m] == padding_residues[m]*reduced_residues[m] % m
             for m in moduli), 'full schedule factorization')
    qo = lcm(*(m//gcd(m, full_residues[m]) for m in moduli))
    qn = lcm(*(m//gcd(m, reduced_residues[m]) for m in moduli))
    need(qn == bank, 'full effective modulus retained by reduced schedule')
    return {'depth': depth, 'input_bit_length_n': integer_receipt(n), 'sOf': s,
            'levels': rows, 'odd_bank': moduli, 'full_scale_residues': full_residues,
            'reduced_scale_residues': reduced_residues, 'old_effective_modulus': qo,
            'new_effective_modulus': qn, 'information_gain_factor': qn//qo,
            'scope': 'Exact schedule arithmetic only; no tables, leaf counts, or input integer N materialized.'}


def hostile():
    path = ROOT/'05-knowledge/results/regular_parent_symmetry_20261029/verify.py'
    spec = importlib.util.spec_from_file_location('inherited_hostile', path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod.hostile_checks()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--output', type=Path)
    args = ap.parse_args()
    result = {'encodings': encodings(), 'observers': observers(),
              'actual_uniform_schedules': [scheduled(1), scheduled(2)],
              'inherited_81_27': hostile(), 'local_checks': CHECKS,
              'scope': 'FINITE-EXACT. Native controls complement, not replace, ReducedEngine.lean.'}
    out = json.dumps(result, indent=2, sort_keys=True)+'\n'
    if args.output:
        args.output.write_text(out, encoding='utf-8', newline='\n')
    else:
        print(out, end='')


if __name__ == '__main__':
    main()
