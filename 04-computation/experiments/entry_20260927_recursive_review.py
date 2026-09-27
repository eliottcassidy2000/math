"""Independent arithmetic audit of the all-k return cylinders and nested rise.
Imports no implementation from the certificate compiler being audited.
"""
from collections import Counter
from fractions import Fraction
from pathlib import Path
import json

checks = Counter()


def need(ok, label):
    checks[label] += 1
    if not ok:
        raise RuntimeError(label)


def strip(n):
    if n <= 0:
        raise ValueError('positive integer required')
    a = 0
    while n % 2 == 0:
        n //= 2
        a += 1
    return n, a


def U(n):
    return strip(3*n+1)[0]


def threshold(k):
    A, B = 8**(k+1), 9**(k+1)
    t = 1
    while 2**t*(A-5) <= B-5:
        t += 1
    return t


def main():
    cylinder_rows = []
    controls = 0
    for k in range(257):
        A, B = 8**(k+1), 9**(k+1)
        t = threshold(k)
        need(2**t*A-B > 5*(2**t-1), 'strict comparison at coefficient1')
        need(t == 1 or 2**(t-1)*(A-5) <= B-5, 'least stated sufficient threshold')
        b0 = 5*pow(B, -1, 2**t) % 2**t
        need(b0 % 2 == 1 and b0 > 0, 'coefficient residue positive odd')
        residue = b0*A-5
        need(0 < residue < A*2**t, 'entire positive cylinder starts in-domain')
        need(Fraction(1, A*2**t) < Fraction(1, B), 'density tail majorant')
        if k < 20:
            cylinder_rows.append(dict(k=k, t=t, b=b0, residue=residue, modulus_power=3*(k+1)+t))
        for j in (*range(8), 10**25+37):
            b = b0+j*2**t
            start = b*A-5
            n = start
            for i in range(k):
                middle, a1 = strip(3*n+1)
                n, a2 = strip(3*middle+1)
                need((a1, a2) == (1, 2), 'literal repeated pair guard')
                need(middle > start and n > start, 'repeated pair above original source')
                need(n == b*8**(k-i)*9**(i+1)-5, 'independent repeated-pair endpoint')
            middle, a1 = strip(3*n+1)
            need(a1 == 1 and middle == 12*b*9**k-7 and middle > start, 'penultimate step and strict rise')
            n, a2 = strip(3*middle+1)
            raw = b*B-5
            expected, extra = strip(raw)
            need(a2 == 2+extra and extra >= t, 'final exact division exponent')
            need(n == expected and n < start, 'exact first descent at2kplus2')
            controls += 1

    need([threshold(k) for k in range(11)] == [1]*5+[2]*6, 'initial threshold blocks')
    hostiles = []
    for k, b, s in [(5, 3, 1), (11, 1, 2)]:
        start = b*8**(k+1)-5
        raw = b*9**(k+1)-5
        end, actual_s = strip(raw)
        need(actual_s == s and s == threshold(k)-1 and end > start, 'lower threshold hostile endpoint')
        n = start
        for _ in range(2*k+2):
            n = U(n)
            need(n > start, 'hostile all preceding values above source')
        need(n == end, 'hostile actual endpoint')
        hostiles.append(dict(k=k, b=b, source=start, endpoint=end, steps=2*k+2))

    # Audit the older bank through its frozen explicit rows, without regenerating it.
    table = Path('05-knowledge/results/reset_20260926_swaplift.out').read_text().splitlines()
    rows = []
    for line in table:
        fields = line.split()
        if len(fields) == 9 and all(f.isdigit() for f in fields):
            rows.append(tuple(map(int, fields)))
    need(len(rows) == 171, 'old bank explicit table universe')
    for q,J,A,P,b,clip,R,residue,power in rows:
        common = 2**min(power, 12)
        need((residue+5) % common != 0, 'old bank disjoint from minus5mod4096')
    for k, q in enumerate((7, 51, 83)):
        prior = next(row for row in rows if row[0] == q)
        fresh = cylinder_rows[k]
        need((prior[-2], prior[-1]) == (fresh['residue'], fresh['modulus_power']), 'first three cylinders exactly old rows')
    old_kept = []
    for row in sorted(rows, key=lambda row: row[-1]):
        residue, power = row[-2:]
        if all(residue % 2**old_power != old_residue for old_residue, old_power in old_kept):
            old_kept.append((residue, power))
    old_mass = sum((Fraction(1, 2**p) for r,p in old_kept), Fraction())
    new_mass = sum((Fraction(1, 2**r['modulus_power']) for r in cylinder_rows), Fraction())
    overlap = sum((Fraction(1, 2**r['modulus_power']) for r in cylinder_rows[:3]), Fraction())
    need(len(old_kept) == 65 and overlap == Fraction(73,1024), 'bank overlap density')
    need(new_mass-overlap == Fraction(2553380107527241,2**64), 'exact first20 added density')
    need(old_mass+new_mass-overlap == Fraction(6990313556829423891,2**65), 'exact augmented density')

    # Fixed-coefficient27 family: literal rise, exact gate, and index spacing.
    gate_data = []
    successful = []
    for k in range(1, 513):
        start = 4*8**k-5
        n = start
        for _ in range(2*k):
            n = U(n)
            need(n > start, 'fixed27 initial rise')
        need(n == 4*9**k-5, 'fixed27 pair endpoint')
        vk = strip(k)[1]
        need(strip(n+1)[1] == 5+vk, 'shifted endpoint valuation')
        length = 4+vk
        for _ in range(length):
            n, a = strip(3*n+1)
            need(a == 1 and n > start, 'fixed27 nested a1 rise')
        if k % 2:
            need(n == (81*9**k-85)//4, 'odd-k nested endpoint')
            end, a = strip(3*n+1)
            need(a == strip(243*9**k-251)[1]-2, 'next gate exact valuation')
            need((243*9**k-251) > 32*start, 'low gate exponent always insufficient')
            if k % 4 == 1:
                need(a == 2 and end > start, 'quarter-class failed gate')
            gate_data.append((k, a))
            if end < start:
                need(a >= 4 and 2**a*16*8**k > 243*9**k, 'successful gate necessary exponential budget')
                successful.append((k, a, end))
    for i, (k, a) in enumerate(gate_data):
        for ell, c in gate_data[i+1:]:
            need(strip(ell-k)[1] >= min(a,c)-1, 'all pairwise valuation spacing')
            if a != c:
                need(strip(ell-k)[1] == min(a,c)-1, 'unequal valuation exact spacing')
    for (k,a,end), (ell,c,other_end) in zip(successful, successful[1:]):
        need(32*(ell-k)*8**k > 243*9**k, 'successful gate exponential spacing')

    print(json.dumps(dict(
        status='INDEPENDENT AUDIT PASS; no universal-entry or policy-totality claim',
        universes=dict(cylinder_k_inclusive=[0,256], coefficients_per_k=9,
                       old_bank_rows=171, fixed27_k_inclusive=[1,512],
                       pairwise_odd_gate_indices=256),
        cylinder_controls=controls,
        first20_added_density=str(new_mass-overlap),
        first20_augmented_density=str(old_mass+new_mass-overlap),
        cylinder_first20=cylinder_rows,
        strict_boundary_hostiles=hostiles,
        successful_next_gates_in_audit=successful,
        checks=dict(sorted(checks.items())), total=sum(checks.values()), result='PASS'
    ), indent=2))


if __name__ == '__main__':
    main()
