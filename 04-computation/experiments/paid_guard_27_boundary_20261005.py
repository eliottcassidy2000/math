"""Exact source-relative boundary exits; no universal or ROOT claim.

Run: python -B this_file.py; repeat with -O. No import-time computation.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from pathlib import Path
import hashlib
import json
from math import gcd

ROOT = Path(__file__).resolve().parents[2]


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(x, minimum=1):
    need(type(x) is int and x >= minimum, 'exact integer in range')


def coefficients(word):
    p = q = 1
    b = 0
    for a in word:
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def step(n):
    integer(n)
    need(n % 2 == 1 and n != 1, 'positive nonroot odd source')
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def replay(n, word):
    for a in word:
        n, seen = step(n)
        need(a == seen, 'literal valuation guard')
    return n


@dataclass(frozen=True)
class Row:
    q: int
    r: int
    a: int
    P: int
    D: int
    B: int
    residue: int
    modulus: int
    minimum: int

    @property
    def word(self):
        return (1, 2)*self.q+(1,)*self.r+(self.a,)


def row(q, r):
    integer(q)
    integer(r)
    P = 3**(2*q+r+1)
    a = max(2, P.bit_length()-3*q-r)
    D = 2**(3*q+r+a)
    B = 5*P-12*3**r*8**q-2**(3*q+r+1)
    need(D > P and a >= 3 and D//2 < P, 'least contracting exit')
    modulus = 2*D
    residue = ((D-B)*pow(P, -1, modulus)) % modulus
    gap = D-P
    first_lift = max(0, (B-gap*residue)//(gap*modulus)+1)
    return Row(q, r, a, P, D, B, residue, modulus,
               residue+first_lift*modulus)


def canonical(item):
    need(type(item) is Row, 'row type')
    need(all(type(v) is int for v in item.__dict__.values()), 'exact row fields')
    need(item == row(item.q, item.r), 'canonical row fields')


def apply(n, item):
    canonical(item)
    integer(n)
    need(n % item.modulus == item.residue, 'exact source cylinder')
    need((item.D-item.P)*n > item.B, 'immutable-source payment')
    y = (item.P*n+item.B)//item.D
    need(0 < y < n and y % 2 == 1, 'positive odd strict descent')
    need(replay(n, item.word) == y, 'independent primitive replay')
    return y


def power_phase(item):
    canonical(item)
    bits = item.modulus.bit_length()-1
    need(item.residue % 8 == 3, 'power-three image')
    e = 1
    for b in range(3, bits):
        candidates = (e, e+2**(b-2))
        valid = [v for v in candidates if pow(3, v, 2**(b+1)) == item.residue % 2**(b+1)]
        need(len(valid) == 1, 'unique exponent lift')
        e = valid[0]
    period = 2**(bits-2)
    # Find the smallest exponent whose actual source satisfies the height guard.
    least_e, power = 0, 1
    while (item.D-item.P)*power <= item.B:
        least_e += 1
        power *= 3
    lift = max(0, (least_e-e+period-1)//period)
    return e, period, e+lift*period


def quarter_child(n):
    integer(n)
    need(n % 256 == 219, 'quarter-child source guard')
    x = replay(n, (1, 2, 1, 1))
    h = (81*n+53)//128
    y, a = step(x)
    z, b = step(h)
    need(x == 4*h+1 and 0 < h < n, 'quarter-child payment')
    need(a >= 3 and b == a-2 and y == z, 'actual common future')
    return h, y, a


def contains(cell, n):
    return n % cell['modulus'] == cell['residue']


def disjoint(x, y):
    return (x.residue-y['residue']) % gcd(x.modulus, y['modulus']) != 0


def rejected(call):
    try:
        call()
    except (ValueError, TypeError):
        return 1
    raise ValueError('hostile was accepted')


def experiment():
    controls = 0
    for q in range(1, 13):
        for r in range(1, 17):
            item = row(q, r)
            need(coefficients(item.word) == (item.P, item.D, item.B), 'closed coefficients')
            need(item.residue % 32 == 27, 'source phase')
            for lift in (0, 1, 7):
                apply(item.minimum+lift*item.modulus, item)
                controls += 1
    new = row(1, 2)
    need((new.P, new.D, new.B, new.residue, new.modulus) == (243, 256, 287, 219, 512), 'new exact cell')
    need(power_phase(new) == (83, 128, 83), 'new exponent cell')
    for t in range(128):
        n = 219+256*t
        h, y, a = quarter_child(n)
        need(h == 139+162*t and y < n, 'coarse all-height source formulas')
        need((a == 3) == (t % 2 == 0), 'new versus inherited exit half')
    for k in (83, 211, 339):
        apply(3**k, new)
    need(pow(3, 27, 256) == 187, 'inherited exponent27 type')
    need(replay(27, (1, 2, 1, 1)) == 71 and 71 > 27, 'literal27 hostile')

    # Freeze comparison through the existing artifacts, without importing their mains.
    oldpath = ROOT/'05-knowledge/results/collatz_paid_portrait_controllers_20261004.json'
    fourpath = ROOT/'05-knowledge/results/collatz_four_slot_compression_20261004.json'
    old, four = json.loads(oldpath.read_text()), json.loads(fourpath.read_text())
    oldcells = old['rational_bank']['rows']+old['switched_bank']['rows']
    oldcells += [x for x in old['safe_cylinders'] if x['c'] == 17]
    oldcells += [old['switched_bank']['checkpoint']]
    oldcells += [{'residue': 187, 'modulus': 256}, {'residue': 7, 'modulus': 256}]
    comparisons = 0
    for q in range(1, 9):
        for r in range(1, 13):
            if (q, r) == (1, 1):
                continue
            item = row(q, r)
            for cell in oldcells:
                need(disjoint(item, cell), 'finite independently reconstructed old-cell exclusion')
                comparisons += 1

    Qmax, Rmax = 24, 64
    selected = [row(q, r) for q in range(1, Qmax+1) for r in range(1, Rmax+1) if (q, r) != (1, 1)]
    lo = sum((F(1, z.D) for z in selected), F(0))
    tail = F(4, 7*2**(3*(Qmax+1)))+F(1, 14*2**Rmax)
    hi = lo+tail
    oldsw = sum((F(c['odd_density']) for c in old['switched_bank']['rows']), F(0))
    need(lo == oldsw+F(1,1024)-F(1,128), 'independent density reconstruction')
    full_boundary = sum(z.minimum != z.residue for z in selected)
    phase_cuts = 0
    for z in selected:
        e, period, emin = power_phase(z)
        need(e % 8 == 3 and pow(3, e, z.modulus) == z.residue, 'entire exponent progression guard')
        phase_cuts += emin != e
    originlo, originhi = map(F, old['coverage']['strengthened_origin_density_interval'])
    prevlo, prevhi = map(F, four['coverage']['named_residual_interval'])
    exlo, exhi = map(F, four['coverage']['exponent_3_mod8_coverage_interval'])
    residual = (prevlo-hi*(1-originlo), prevhi-lo*(1-originhi))
    exponents = (exlo+16*lo, exhi+16*hi)
    need(F('0.53744803446924') < exponents[0] < exponents[1] < F('0.53744803446925'), 'outward exponent interval')
    need(F('0.10584246338775') < residual[0] < residual[1] < F('0.10584246338776'), 'outward residual interval')
    hostiles = sum((rejected(lambda: row(True, 2)), rejected(lambda: row(1, 2.0)),
                    rejected(lambda: apply(219, replace(new, B=285))),
                    rejected(lambda: apply(27, new)), rejected(lambda: quarter_child(27)),
                    rejected(lambda: apply(219.0, new))))
    print('Type distinction: exponent27mod64 -> source187mod256; literal source27mod64 starts1211.')
    print('New coarse family:219+256t -> quarter-child139+162t; actual1211a, a>=3 pays source.')
    print('New exact half:219+512t ->209+486t via12113; powers3^(83+128t).')
    print('All q,r>=1: first contracting exit a0 plus exact source guard(D-P)n>B gives infinite paid AP tail.')
    print('No assumption that every uncut least source pays; inherited128x128 probe remains finite evidence.')
    print('Independent literal controls q1..12,r1..16,lifts0,1,7:', controls)
    print('Quarter-child controls:', 128, '; explicit power sources:', 3)
    print('Prior-cell comparison q1..8,r1..12 except1,1:', comparisons)
    print('Density truncation q1..24,r1..64 except1,1:', len(selected), 'rows; binary/phase initial cuts:', full_boundary, phase_cuts)
    print('Added odd density:', float(lo), '; exact tail bound:', tail)
    print('Named exponent coverage:0.53744803446924 < delta <0.53744803446925.')
    print('Named fused residual:0.10584246338775 < delta <0.10584246338776.')
    print('Frozen prior SHA256:', hashlib.sha256(oldpath.read_bytes()).hexdigest())
    print('Frozen four-slot SHA256:', hashlib.sha256(fourpath.read_bytes()).hexdigest())
    print('Typed/source/corrupted-row hostiles rejected:', hostiles)
    print('Payment is a smaller dependency; universal guard coverage and completed ROOT coverage remain OPEN.')


if __name__ == '__main__':
    experiment()
