"""Earlier guarded sibling reroutes with their exact least inverse clock.

The algebraic proofs and finite universes are in the companion note. An output
dependency is not a root certificate unless its actual child proof is supplied.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from functools import lru_cache
import json

import collatz_complement_routing_20261004 as previous
import collatz_branch_toll_rank_20261004 as energy
import inverse_ray_ternary_addresses_20261004 as codec


CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def odd(n):
    need(type(n) is int and n > 0 and n % 2 == 1, 'positive exact odd integer')


def v2(n):
    need(type(n) is int and n > 0, 'positive valuation argument')
    return (n & -n).bit_length()-1


def step(n):
    odd(n)
    z = 3*n+1
    a = v2(z)
    return z >> a, a


def literal(n, word):
    """Independent division replay, retaining first-hit root semantics."""
    odd(n)
    states = [n]
    for a in word:
        need(type(a) is int and a >= 1 and n != 1, 'exact word before first root')
        z, actual = 3*n+1, 0
        while z % 2 == 0:
            z //= 2
            actual += 1
        need(a == actual, 'literal valuation guard')
        n = z
        states.append(n)
    return tuple(states)


@lru_cache(None, typed=True)
def clock(k):
    need(type(k) is int and k >= 1, 'positive sibling index')
    d = 1
    while 3**d <= 2**(d+2*k-1):
        d += 1
    return d


@dataclass(frozen=True)
class Row:
    position: int
    k: int
    d: int
    M: int
    residue: int
    P: int
    C: int
    intercept: int

    @property
    def ell(self):
        return self.position+self.d


@lru_cache(None, typed=True)
def row(position, k):
    need(type(position) is int and position >= 0, 'nonnegative exact checkpoint')
    need(type(k) is int and k >= 1 and k % 3**position == 0,
         'position-specific ternary feasibility')
    d = clock(k)
    M, P = 3**d, 2**(d+2*k-1)
    raw = 2**(position+1)*(4**k-1)
    need(raw % 3**(position+1) == 0, 'cancel only the proved ternary factor')
    C = raw//3**(position+1)
    residue = (C*pow(4**k, -1, M)-1) % M
    intercept = P-2**(d-1)*C-M
    need(P < M and intercept < 0, 'strict uniform source payment')
    return Row(position, k, d, M, residue, P, C, intercept)


def apply(n, spec):
    """Return the exact join, or None; no orbit or child-certificate search."""
    odd(n)
    need(type(spec) is Row and all(type(x) is int for x in spec.__dict__.values()),
         'typed reroute row')
    need(spec == row(spec.position, spec.k), 'canonical row coefficients and guards')
    if n == 1 or spec.position > v2(n+1)-1 or n % spec.M != spec.residue:
        return None
    numerator = spec.P*n+spec.intercept
    need(numerator % spec.M == 0, 'integer child')
    child = numerator//spec.M
    x = 3**spec.position*(n+1)//2**spec.position-1
    target, a = step(x)
    source_word = (1,)*spec.position+(a,)
    child_word = (1,)*(spec.ell-1)+(2, a+2*spec.k-2)
    need(0 < child < n and child % 2 and v2(child+1) == spec.ell,
         'positive odd child, original-source payment, exact initial run')
    return dict(source=n, child=child, endpoint=target, position=spec.position,
                k=spec.k, d=spec.d, source_word=source_word, child_word=child_word)


def candidates(n):
    """Complete finite applicability test for this whole early-reroute grammar."""
    odd(n)
    if n == 1:
        return ()
    cap = n.bit_length()-1
    result = []
    for r in range(v2(n+1)):
        if 2*3**r > cap-r:
            break
        for k in range(3**r, (cap-r)//2+1, 3**r):
            if r+clock(k) > cap:
                break
            hit = apply(n, row(r, k))
            if hit is not None:
                result.append(hit)
    return tuple(result)


def height_obstruction(n):
    """Proved sufficient failure predicate for the entire early-reroute grammar."""
    odd(n)
    unit, ternary = n, 1
    while unit % 3 == 0:
        unit //= 3
        ternary *= 3
    cap = n.bit_length()-1
    return ternary*2**cap >= 2*3**cap


def progression(spec, binary=123, period=128):
    need(type(binary) is int and type(period) is int and period > 0
         and period & (period-1) == 0 and binary % 2 == 1,
         'dyadic source cylinder')
    n = spec.residue+spec.M*((binary-spec.residue)*pow(spec.M, -1, period) % period)
    hit = apply(n, spec)
    need(hit is not None, 'binary cylinder retains the requested initial ones')
    return n, period*spec.M, hit['child'], period*spec.P


def audit_hit(hit):
    left = literal(hit['source'], hit['source_word'])
    right = literal(hit['child'], hit['child_word'])
    need(left[-1] == right[-1] == hit['endpoint'], 'independent actual common future')
    x = left[-2]
    sibling = 4**(hit['k']-1)*x+(4**(hit['k']-1)-1)//3
    need(right[-2] == sibling, 'retained intermediate sibling, not a substituted source')


def attach_supplied_child(hit, certificate):
    codec.audit_certificate(certificate)
    need(codec.expand(certificate, bit_cap=codec.bit_bounds(certificate)[1]) == hit['child'],
         'supplied first-hit certificate has the exact child identity')
    suffix = certificate
    for a in hit['child_word']:
        need(suffix != codec.ROOT and codec.exponent(suffix) == a, 'consume actual child prefix')
        suffix = suffix.parent
    need(codec.expand(suffix, bit_cap=codec.bit_bounds(suffix)[1]) == hit['endpoint'],
         'supplied common-future suffix')
    states = literal(hit['source'], hit['source_word'])
    for i in range(len(hit['source_word'])-1, -1, -1):
        source, target = states[i:i+2]
        a = hit['source_word'][i]
        least = codec.kappa(target % 9, source % 3)
        need(a >= least and (a-least) % 6 == 0, 'exact inverse certificate address')
        suffix = codec.extend(suffix, source % 3, (a-least)//6)
    need(codec.expand(suffix, bit_cap=codec.bit_bounds(suffix)[1]) == hit['source'],
         'retained original source certificate')
    return suffix


def phase_log2(target, depth):
    """Exact ternary digit lift, using ord_(3^d)(2)=2*3^(d-1)."""
    need(type(depth) is int and depth >= 1 and type(target) is int and target % 3,
         'unit target at positive ternary precision')
    exponent = 0 if target % 3 == 1 else 1
    for level in range(2, depth+1):
        modulus, old_period = 3**level, 2*3**(level-2)
        choices = [exponent+j*old_period for j in range(3)
                   if pow(2, exponent+j*old_period, modulus) == target % modulus]
        need(len(choices) == 1, 'unique next ternary exponent digit')
        exponent = choices[0]
    return exponent, 2*3**(depth-1)


def ternary_density(rows):
    primitive, total = [], Fraction(0)
    for key, depth, residue in sorted(rows, key=lambda item: item[1]):
        if not any((residue-old[2]) % 3**old[1] == 0 for old in primitive):
            primitive.append((key, depth, residue))
            total += Fraction(1, 3**depth)
    return total, tuple(item[0] for item in primitive)


def decimal_interval(lo, hi, places=16):
    scale = 10**places
    lower = lo.numerator*scale//lo.denominator
    upper = -((-hi.numerator*scale)//hi.denominator)
    def fmt(n):
        return f'{n//scale}.{n%scale:0{places}d}'
    return fmt(lower), fmt(upper)


def main():
    print('EARLY REROUTES: exact least branch clock and guarded larger coverage; universal home OPEN')
    bank = previous.load_bank()
    controls = 0
    for r in range(4):
        binary, B = previous.debt.cylinder((1,)*r+(2,))
        for u in (1, 2, 4):
            spec = row(r, 3**r*u)
            n0, period, h0, increment = progression(spec, binary, B)
            for t in (0, 1, 17):
                n, h = n0+period*t, h0+increment*t
                hit = apply(n, spec)
                need(hit is not None and hit['child'] == h, 'full arithmetic family')
                audit_hit(hit)
                need(any(z['position'] == r and z['k'] == spec.k for z in candidates(n)),
                     'finite applicability search retains this exact family')
                controls += 1
    print('All-height schema controls:', controls, '(r0..3,u1/2/4,three parameters; k=3^r*u)')

    # Independent finite check of the iff clock: every legal shallower construction
    # fails payment; all legal minimal-clock constructions pay it.
    shallow = 0
    for r in range(3):
        for k in range(1, 13):
            for ell in range(1, min(r+clock(k)+1, 13)):
                for n in range(3, 400, 2):
                    if r > v2(n+1)-1:
                        continue
                    x = 3**r*(n+1)//2**r-1
                    S = 4**k*x+(4**k-1)//3
                    if (S+1) % 3**ell:
                        continue
                    h = 2**(ell-1)*((S+1)//3**ell)-1
                    if h < n:
                        need(k % 3**r == 0 and ell-r >= clock(k), 'necessary feasibility and clock')
                    else:
                        need(ell-r < clock(k), 'all legal paid-clock rows really pay')
                    shallow += 1
    print('Independent legal clock comparisons:', shallow, '(n3..399,r0..2,k1..12,ell<=12)')

    origin, one_step = row(0, 2), row(1, 3)
    first = progression(origin)
    later = progression(one_step)
    need(first == (89211, 93312, 62655, 65536), 'origin family constants')
    need(later == (2424699, 2519424, 2018303, 2097152), 'one-step family constants')
    print('ORIGIN [source,period,child,increment]:', first, ';words(1) versus1^5,2,3')
    print('ONE_STEP [source,period,child,increment]:', later, ';words(1,2) versus1^9,2,6')
    refined = []
    for digit in (0, 1):
        n0, h0 = first[0]+digit*first[1], first[2]+digit*first[3]
        refined.append((n0, 3*first[1], h0, 3*first[3]))
        for j in range(3):
            _, M, residue = previous.old_row(j)
            need(n0 % M != residue, 'first three old ternary guards excluded')
        for t in range(64):
            n = n0+3*first[1]*t
            hit = apply(n, origin)
            audit_hit(hit)
            need(all(x > n for x in literal(n, (1, 2, 1, 2))[1:]), 'first four actual steps all grow')
            need(previous.fusion.select_debt(n, bank) is None, 'binary16 exclusion')
            need(n % 3 == 0 and energy.critical(n) and energy.rank(hit['child']) < energy.rank(n),
                 'rank cancellation at a critical source with no odd predecessor')
    need((first[0]+2*first[1]) % 3**7 == previous.old_row(2)[2], 'third digit is inherited')
    old_tail = Fraction(1, 78)
    prior_tail = Fraction(1, 3**55)/(1-Fraction(1, 3**30))
    fresh = 1-old_tail-prior_tail
    need(fresh > Fraction(98717, 100000), 'strict positive surviving fraction')
    prior10 = previous.family(10)
    need(first[0] % 9 == 3 and prior10.residue % 9 == 0, 'separate ternary source classes')
    need(one_step.residue != prior10.residue % one_step.M, 'one-step family differs from prior k10')
    need(previous.family(19).ell-2 == 62, 'first conservative prior-composition tail depth')
    print('REFINED_ORIGIN_ROWS', json.dumps(refined))
    print('Each refined row:oldbank occupancy<=1/78; previous composition extra<=3^-55/(1-3^-30).')
    print('Surviving relative density >0.98717; union odd-relative density lower bound:',
          str(2*fresh/Fraction(64*3**7)))

    # The position-one family supplies an additional cell outside the strengthened
    # origin bank; the first four origin rows disagree at available precision.
    for k in range(4):
        if k == 0:
            d, residue = 1, 2
        else:
            q = row(0, k)
            d, residue = q.d, q.residue
        need(one_step.residue % 3**d != residue, 'position-one excludes first four origin rows')
    need(clock(4) == 12, 'first remaining origin-bank depth')
    need(Fraction(3**9, 3**12)/(1-Fraction(1, 27)) == Fraction(1, 26), 'position-one tail')
    print('Position-one row:at least25/26 lies outside the full strengthened origin bank.')
    phases = []
    for spec in (origin, one_step):
        slope = spec.P*pow(spec.M, -1, 19) % 19
        intercept = spec.intercept*pow(spec.M, -1, 19) % 19
        n0, period, _, _ = progression(spec)
        observed = set()
        for t in range(19):
            n = n0+period*t
            hit = apply(n, spec)
            h = hit['child']
            need((h-slope*n-intercept) % 19 == 0, 'row-indexed affine phase transport')
            need((pow(slope, -1, 19)*(h-intercept)-n) % 19 == 0,
                 'retained row makes phase transport reversible')
            need(energy.critical(n) and energy.rank(h) < energy.rank(n),
                 'both early rows cancel original critical rank')
            observed.add(n % 19)
        need(len(observed) == 19, 'all source phases occur in each family')
        phases.append((spec.position, spec.k, slope, intercept))
    print('MOD19 [position,k,slope,intercept]:', phases, ';invertible transport, not identity assumed')

    old_rows, new_rows = [], [(0, 1, 2)]
    for k in range(32):
        d, _, residue = previous.old_row(k)
        old_rows.append((k, d, residue))
        if k:
            spec = row(0, k)
            new_rows.append((k, spec.d, spec.residue))
            need(previous.old_row(k)[0]-spec.d in (1, 2), 'one or two fewer ternary digits')
    prior_density, _ = ternary_density(old_rows+[('inverse12', 2, 4)])
    new_density, primitive = ternary_density(new_rows)
    prior_error = Fraction(27, 26*3**previous.old_row(32)[0])
    new_error = Fraction(27, 26*3**clock(32))
    print('Origin-bank density interval:', decimal_interval(new_density, new_density+new_error))
    print('Oldbank+inherited12 density interval:', decimal_interval(prior_density, prior_density+prior_error))
    print('Origin-bank primitive indices through31:', primitive)
    need(new_density > prior_density+prior_error, 'strict coverage increase after inherited12 credit')

    mersennes = []
    for k in (1, 4, 7):
        spec = row(0, k)
        exponent, period = phase_log2(spec.residue+1, spec.d)
        for t in (0, 1, 17):
            p = exponent+period*t
            need((pow(2, p, spec.M)-1) % spec.M == spec.residue, 'symbolic Mersenne guard')
            need(p % 2 == 1 and p > spec.d, 'hard-Mersenne reset-two type and sufficient height')
        mersennes.append((k, spec.d, exponent, period, exponent % 6))
    need(mersennes == [(1, 2, 5, 6, 5), (4, 12, 75379, 354294, 1),
                       (7, 23, 36969036897, 62762119218, 3)], 'exact lifted exponent addresses')
    p = mersennes[1][2]
    huge = (1 << p)-1
    literal_large = apply(huge, row(0, 4))
    audit_hit(literal_large)
    need(literal_large['child'] < huge and huge.bit_length() == 75379, 'actual large Mersenne reduction')
    print('MERSENNE [k,d,p0,p_period,p_mod6]:', mersennes)
    print('Literal large control:75379-bit source,13 child odd edges; largest exponent family stays symbolic.')

    power_controls = 0
    for a in range(1, 33):
        n = 3**a
        need(height_obstruction(n) and not candidates(n), 'all-height powers-of-three hostile')
        power_controls += 1
    wedge = []
    for unit in (1, 5, 7, 11, 19, 101):
        a = 1
        while 9**a < 32*unit**3:
            a += 1
        wedge.append((unit, a))
        for lift in range(4):
            n = 3**(a+lift)*unit
            need(height_obstruction(n) and not candidates(n), 'whole fixed-cofactor tower tail excluded')
    print('Infinite hostile:all powers3^a,a>=1; literal controls:', power_controls)
    print('Ternary wedge9^a>=32u^3; [u,first sufficient a]:', wedge, ';four heights each')
    divisions = dict(direct=0, half_child=0, three_steps=0, residual=0)
    for a in range(1, 65):
        n = 3**a
        if a % 2 == 0:
            need(step(n)[0] < n, 'even-exponent immediate source descent')
            divisions['direct'] += 1
        elif a % 4 == 1:
            x, one = step(n)
            _, reset = step(x)
            need(one == 1 and reset >= 3, 'inherited half-child reset applies')
            divisions['half_child'] += 1
        elif a % 8 == 7:
            x = literal(n, (1, 2))[-1]
            endpoint, last = step(x)
            need(last >= 2 and endpoint < n, 'three-step source-relative descent')
            divisions['three_steps'] += 1
        else:
            need(a % 8 == 3 and n % 96 == 27, 'remaining critical exponent class')
            divisions['residual'] += 1
    print('Other constructions on powers3^a,a1..64:', divisions,
          ';all-height survivor for these three rules is a=3mod8.')

    for spec, n in ((origin, first[0]), (origin, refined[1][0]), (one_step, later[0])):
        hit = apply(n, spec)
        child = codec.encode_source(hit['child'], step_cap=5000)
        result = attach_supplied_child(hit, child)
        need(result == codec.encode_source(n, step_cap=5000), 'independent complete source certificate')
    print('Supplied-child demonstrations:3; premise search cap5000, not universal child grounding.')
    for n in (1, 7, 27, 703):
        need(candidates(n) == (), 'retained root or uncovered-source controls')
    rejected = 0
    for thunk in (lambda: candidates(True), lambda: candidates(7.0), lambda: candidates(-7),
                  lambda: candidates(8), lambda: row(1, 2), lambda: row(False, 1),
                  lambda: clock(True),
                  lambda: apply(121, replace(row(0, 1), intercept=67)),
                  lambda: attach_supplied_child(apply(first[0], origin), codec.ROOT)):
        try:
            thunk()
        except ValueError:
            rejected += 1
        else:
            raise ValueError('invalid API input accepted')
    print('Hostiles:7/27/703 still uncovered; k1 inherited12; new type/source rejections:', rejected)
    print('PASS:', CHECKS, 'always-active checks; no coverage or home claim beyond the proved guards.')


if __name__ == '__main__':
    main()
