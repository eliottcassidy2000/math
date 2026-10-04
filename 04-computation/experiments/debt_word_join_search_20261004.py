"""Exact affine-word search for reset-two common futures.

Complete source-word universe: all positive compositions of total cost 1..16.
No numerical trajectory search supplies a proof of an infinite family.
"""
from fractions import Fraction
from itertools import combinations
import json
import inverse_ray_ternary_addresses_20261004 as codec


SELECTED = {
    (1, 1): ((1, 6, 1), (1, 2, 1, 3)),
    (1, 2): ((1, 6, 4, 1), (2, 1, 1, 1, 3, 3)),
    (1, 3): ((1, 2, 9, 1), (1, 1, 1, 1, 1, 2, 5)),
    (1, 4): ((1, 14, 1), (2, 2, 1, 3, 3, 1, 3)),
    (2, 1): ((6, 1), (1, 1, 3)),
    (2, 2): ((10, 1), (3, 2, 1, 3)),
    (2, 3): ((10, 1), (1, 1, 1, 3, 3)),
    (2, 4): ((3, 1, 11, 1), (1, 1, 2, 1, 1, 1, 2, 5)),
}


def need(ok, message):
    if not ok:
        raise ValueError(message)


def compositions(total):
    if total == 0:
        yield ()
    for first in range(1, total+1):
        for tail in compositions(total-first):
            yield (first,)+tail


def data(word):
    p = q = 1
    b = 0
    for a in word:
        need(type(a) is int and a >= 1, 'positive exact valuation')
        p, q, b = 3*p, q*2**a, 3*b+q
    return p, q, b


def decode(carry, cost, length):
    """Unique positive word, if one exists, with this length/cost/carry."""
    if length < 1 or cost < length:
        return None
    word = []
    for remaining in range(length, 1, -1):
        difference = carry-3**(remaining-1)
        if difference <= 0:
            return None
        exponent = (difference & -difference).bit_length()-1
        if exponent < 1:
            return None
        word.append(exponent)
        cost -= exponent
        carry = difference >> exponent
    if carry != 1 or cost < 1:
        return None
    return tuple(word)+(cost,)


def cylinder(word):
    p, q, b = data(word)
    return (q-b)*pow(p, -1, 2*q) % (2*q), 2*q


def replay(n, word):
    need(type(n) is int and n > 0 and n % 2, 'positive odd input')
    path = [n]
    for a in word:
        need(n != 1, 'first root stop')
        numerator = 3*n+1
        actual = (numerator & -numerator).bit_length()-1
        need(actual == a, 'literal exact valuation')
        n = numerator >> a
        path.append(n)
    return path


def words_at_reset(r, e, s, u, v):
    return (1,)*r+(2,)*e+u, (1,)*(r-1)+(2*e+s,)+v


def select_join(n):
    """Recognize a supplied source; return a dependency, never invent its proof."""
    need(type(n) is int and n > 0 and n % 2, 'positive exact odd source')
    run = ((n+1) & -(n+1)).bit_length()-2
    if n == 1 or run < 1:
        return None
    t = (n+1) >> (run+1)
    raw = 3**run*t-1
    child_reset = (raw & -raw).bit_length()  # 1+v2(raw).
    if child_reset < 3:
        return None
    e = (child_reset-1)//2
    s = child_reset-2*e
    selected = SELECTED.get((s, e))
    if selected is None:
        return None
    u, v = selected
    prefix = (1,)*run+(2,)*e+u[:-1]
    residue, period = cylinder(prefix)
    if n % period != residue:
        return None
    p, q, b = data(prefix)
    w = (p*n+b)//q
    numerator = 3*w+1
    a = (numerator & -numerator).bit_length()-1
    parent, child = words_at_reset(run, e, s, u[:-1]+(a,),
                                    v[:-1]+(a+v[-1]-u[-1],))
    return dict(source=n, r=run, s=s, e=e, a=a, child=(n-1)//2,
                endpoint=numerator >> a, parent_word=parent, child_word=child)


def transport_certificate(n, child_certificate):
    """Consume an existing child AST; this function performs no home search."""
    record = select_join(n)
    need(record is not None, 'source matches a proved join')
    codec.audit_certificate(child_certificate)
    need(codec.expand(child_certificate) == record['child'], 'actual child source identity')
    suffix = child_certificate
    for exponent in record['child_word']:
        need(suffix != codec.ROOT and codec.exponent(suffix) == exponent,
             'supplied child certificate contains the exact join word')
        suffix = suffix.parent
    need(codec.expand(suffix) == record['endpoint'], 'rooted common endpoint suffix')
    path = replay(n, record['parent_word'])
    need(path[-1] == record['endpoint'], 'original source reaches common endpoint')
    for i in range(len(record['parent_word'])-1, -1, -1):
        row = path[i] % 3
        least = codec.kappa(path[i+1] % 9, row)
        exponent = record['parent_word'][i]
        need(exponent >= least and (exponent-least) % 6 == 0, 'inverse edge guard')
        suffix = codec.extend(suffix, row, (exponent-least)//6)
    need(codec.expand(suffix) == n, 'transported original source identity')
    odd_rank, ordinary_rank = codec.ranks(child_certificate)
    need(codec.ranks(suffix) == (odd_rank, ordinary_rank+1), 'exact transported ranks')
    return suffix


def main():
    bank = {}
    counts = {(s, e): 0 for s in (1, 2) for e in range(1, 5)}
    decoded = enumerated = 0
    for cost in range(1, 17):
        for u in compositions(cost):
            enumerated += 1
            p, q, b = data(u)
            need(decode(b, cost, len(u)) == u, 'independent carry round trip')
            decoded += 1
            for s in (1, 2):
                if (s == 1 and u[0] != 1) or (s == 2 and u[0] < 3):
                    continue  # Exact parity of Y=2^s*3^e*M+1 forces these guards.
                for e in range(1, 5):
                    need((p+b) % 2**s == 0, 'candidate integer carry')
                    v = decode((p+b)//2**s, cost-s, len(u)+e)
                    if v is None:
                        continue
                    pv, qv, bv = data(v)
                    need(pv == p*3**e and q == 2**s*qv and 2**s*bv == p+b,
                         'independent affine common-future identity')
                    counts[s, e] += 1
                    if (s, e) not in bank:
                        bank[s, e] = (u, v)
    need(enumerated == decoded == 65535, 'complete positive-composition universe')
    need(bank == SELECTED, 'selected minimal-cost joins')

    # Independent decoder control: construct words by separator subsets.
    for cost in range(1, 11):
        for length in range(1, cost+1):
            for cuts in combinations(range(1, cost), length-1):
                ends = (0,)+cuts+(cost,)
                word = tuple(y-x for x, y in zip(ends, ends[1:]))
                p, q, b = data(word)
                direct = sum(3**(length-i-1)*2**sum(word[:i]) for i in range(length))
                need(direct == b and decode(b, cost, length) == word,
                     'separator generation and direct carry independently agree')

    rows = []
    controls = 0
    for (s, e), (u, v) in sorted(bank.items()):
        delta = v[-1]-u[-1]
        need(delta >= 2 and delta % 2 == 0, 'positive sibling scale')
        pre_p, pre_q, pre_b = data(u[:-1])
        root_m = Fraction(pre_q-pre_p-pre_b, 2**s*3**e*pre_p)
        need(root_m.denominator != 1, 'source preterminal root requires noninteger M')
        child_p, child_q, child_b = data(v[:-1])
        need(Fraction(child_p, child_q) == 2**delta*Fraction(2**s*3**e*pre_p, pre_q)
             and Fraction(child_b, child_q) ==
             2**delta*Fraction(pre_p+pre_b, pre_q)+Fraction(2**delta-1, 3),
             'child preterminal is larger sibling, with correct orientation')
        for a in (1, 2, 3, 7):
            ua, va = u[:-1]+(a,), v[:-1]+(a+delta,)
            residue, period = cylinder(va)
            for k in range(4):
                m = residue+period*k
                left = replay(2**s*3**e*m+1, ua)
                right = replay(m, va)
                need(left[-1] == right[-1], 'literal boundary common future')
                controls += 1
            for r in (1, 2, 8, 13, 21):
                parent, child = words_at_reset(r, e, s, ua, va)
                p, q, b = data(parent)
                pc, qc, bc = data(child)
                need(p == pc and q == 2*qc and b == 2*bc-p,
                     'whole source-to-smaller-child affine identity')
                residue, period = cylinder(parent)
                for k in range(4):
                    n = residue+period*k
                    m = (n-1)//2
                    left, right = replay(n, parent), replay(m, child)
                    need(0 < m < n and left[-1] == right[-1], 'guarded smaller dependency')
                    need(len(parent) == len(child) and sum(parent) == sum(child)+1,
                         'odd rank preserved, ordinary rank increases by one')
                    controls += 1
        constant = 2*e+sum(u[:-1])+1
        growth = next(r for r in range(1, 100)
                      if all(data(words_at_reset(r, e, s, u, v)[0][:j])[0] >
                             data(words_at_reset(r, e, s, u, v)[0][:j])[1]
                             for j in range(1, r+e+len(u)+1)))
        parent, child = words_at_reset(growth, e, s, u, v)
        n, period = cylinder(parent)
        actual = replay(n, parent)
        need(all(x > n for x in actual[1:]), 'growing positive source prefix')
        if growth > 1:
            hostile, _ = words_at_reset(growth-1, e, s, u, v)
            need(any(data(hostile[:j])[0] < data(hostile[:j])[1]
                     for j in range(1, len(hostile)+1)), 'sharp coefficient growth cutoff')
        rows.append(dict(s=s, e=e, u=u, v=v, minimum_cost=sum(u),
                         joins_through_cost16=counts[s, e], last_offset=delta,
                         preterminal_root_requires_M=str(root_m),
                         union_cylinder_exponent_constant=constant,
                         absolute_union_density=str(Fraction(1, 2**constant)),
                         first_all_prefix_growth_r=growth,
                         growing_example=n, child=(n-1)//2, common_future=actual[-1],
                         growing_word=parent))

    # Distinct (s,e) families have different child first-reset values 2e+s.
    # Their tails have summably vanishing upper density via initial runs.
    density = sum((Fraction(row['absolute_union_density']) for row in rows), Fraction())
    kernel = Fraction(1, 8)-density
    selected_count = 0
    for n in range(1, 20001, 2):
        observed, current, run = [], n, 0
        # Independent literal parser: initial ones, then seven more letters.
        while current != 1:
            numerator = 3*current+1
            a = (numerator & -numerator).bit_length()-1
            observed.append(a)
            current = numerator >> a
            if a != 1:
                break
            run += 1
        for _ in range(6):
            if current == 1:
                break
            numerator = 3*current+1
            a = (numerator & -numerator).bit_length()-1
            observed.append(a)
            current = numerator >> a
        literal = []
        if run >= 1:
            for (s, e), (u, v) in SELECTED.items():
                prefix = (1,)*run+(2,)*e+u[:-1]
                if tuple(observed[:len(prefix)]) == prefix:
                    literal.append((s, e))
        selected = select_join(n)
        need(len(literal) <= 1 and bool(literal) == (selected is not None),
             'independent literal selector applicability and disjointness')
        if selected is not None:
            need(literal == [(selected['s'], selected['e'])], 'same selected boundary')
            need(replay(n, selected['parent_word'])[-1] ==
                 replay(selected['child'], selected['child_word'])[-1] == selected['endpoint'],
                 'selected join independently replayed')
            selected_count += 1

    for row in rows:
        n = row['growing_example']
        # Only these bounded controls search for a premise; transport never does.
        supplied = codec.encode_source((n-1)//2, step_cap=5000)
        transported = transport_certificate(n, supplied)
        codec.literal_certificate_check(transported)
        need(transported == codec.encode_source(n, step_cap=5000), 'independent canonical route')
    for bad in (True, 0, -1, 2):
        try:
            select_join(bad)
        except ValueError:
            pass
        else:
            raise ValueError('malformed selector source accepted')
    try:
        transport_certificate(rows[0]['growing_example'], codec.ROOT)
    except ValueError:
        pass
    else:
        raise ValueError('unrelated supplied child certificate accepted')
    print(json.dumps(dict(status='PROVED affine joins and guarded families; FINITE-EXACT bounded minima',
                          source_words=enumerated, maximum_source_cost=16,
                          controls=controls, families=rows,
                          selector_universe='all10000 odd sources through20000',
                          selected_sources=selected_count, supplied_child_AST_transports=len(rows),
                          malformed_sources_rejected=4, unrelated_child_certificate_rejected=True,
                          absolute_removed_density=str(density),
                          absolute_remaining_necessary_kernel=str(kernel),
                          odd_relative_remaining_necessary_kernel=str(2*kernel),
                          scope='Smaller dependency requires its own home certificate; not universal coverage'),
                     indent=2))
    print('PASS: all checks remain active under -O')


if __name__ == '__main__':
    main()
