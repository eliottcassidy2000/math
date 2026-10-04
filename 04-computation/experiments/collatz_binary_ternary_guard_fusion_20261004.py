"""Compose the incoming dyadic debt bank with exact ternary sibling guards.

Independent literal replay, exhaustive minimal-cost search, and finite CRT
controls. A necessary least-counterexample domain is not a convergence proof.
"""
from fractions import Fraction
from pathlib import Path
import json

import debt_word_join_search_20261004 as inherited

CHECKS = 0


def check(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ArithmeticError(why)


def val(n, p):
    check(n != 0, 'finite valuation argument')
    a = 0
    while n % p == 0:
        n //= p
        a += 1
    return a


def literal(n, word):
    for a in word:
        z, actual = 3*n+1, 0
        while z % 2 == 0:
            z //= 2
            actual += 1
        check(actual == a, 'independent actual word')
        n = z
    return n


def sibling_row(k):
    ell = 1
    while 3**ell <= 2**(ell+2*k):
        ell += 1
    modulus = 3**ell
    residue = -((4**k+2)//3)*pow(4**k, -1, modulus) % modulus
    return dict(k=k, ell=ell, modulus=modulus, residue=residue)


def sibling_child(n):
    # The proof gives ell+1 < log2(n+1), hence this safe finite k bound.
    for k in range((n+1).bit_length()//2+1):
        row = sibling_row(k)
        if n % row['modulus'] != row['residue']:
            continue
        S = 4**k*n+(4**k-1)//3
        m = 2**row['ell']*((S+1)//row['modulus'])-1
        check(0 < m < n and literal(m, (1,)*row['ell']) == S,
              'sibling rank and exact route')
        check((3*S+1)//(3*n+1) == 4**k, 'same odd future')
        return dict(child=m, **row)
    return None


def select_debt(n, bank):
    """Recognize supplied n using the discovered words; no child proof is invented."""
    check(type(n) is int and n > 0 and n % 2 == 1, 'selector source type')
    r = val(n+1, 2)-1
    if n == 1 or r < 1:
        return None
    t = (n+1)//2**(r+1)
    reset = 1+val(3**r*t-1, 2)
    if reset < 3:
        return None
    e, s = (reset-1)//2, 1+(reset-1) % 2
    if (s, e) not in bank:
        return None
    u, v = bank[s, e]
    prefix = (1,)*r+(2,)*e+u[:-1]
    residue, period = inherited.cylinder(prefix)
    if n % period != residue:
        return None
    p, q, b = inherited.data(prefix)
    x = (p*n+b)//q
    a = val(3*x+1, 2)
    w, z = inherited.words_at_reset(r, e, s, u[:-1]+(a,),
                                   v[:-1]+(a+v[-1]-u[-1],))
    return dict(source=n, child=(n-1)//2, s=s, e=e, r=r,
                parent_word=w, child_word=z, endpoint=(3*x+1)//2**a)


def find_debt_bank():
    bank, counts = {}, []
    enumerated = 0
    for cost in range(1, 25):
        count = 0
        for u in inherited.compositions(cost):
            enumerated += 1
            count += 1
            s = 1 if u[0] == 1 else 2 if u[0] >= 3 else 0
            es = [e for e in range(1, 9)
                  if s and (s, e) not in bank and cost-s >= len(u)+e]
            if not es:
                continue
            p, q, b = inherited.data(u)
            check((p+b) % 2**s == 0, 'boundary carry integral')
            for e in es:
                v = inherited.decode((p+b)//2**s, cost-s, len(u)+e)
                if v is not None:
                    pv, qv, bv = inherited.data(v)
                    check(pv == p*3**e and q == 2**s*qv and p+b == 2**s*bv,
                          'independent affine identity')
                    bank[s, e] = (u, v)
            if len(bank) == 16:
                break
        counts.append(dict(cost=cost, words=count))
        if len(bank) == 16:
            break
    check(len(bank) == 16 and counts[-1]['cost'] == 24, 'full requested debt bank')
    check(all(counts[c-1]['words'] == 2**(c-1) for c in range(1, 24)),
          'complete lower-cost universe')
    check(all(bank[key] == value for key, value in inherited.SELECTED.items()),
          'recover incoming eight controls')
    rows = []
    for (s, e), (u, v) in sorted(bank.items()):
        c = 2*e+sum(u[:-1])+1
        example = None
        for r in (1, 2, 8, 31):
            for a in (1, 2, 7):
                uw, vw = u[:-1]+(a,), v[:-1]+(a+v[-1]-1,)
                w, z = inherited.words_at_reset(r, e, s, uw, vw)
                residue, period = inherited.cylinder(w)
                for t in (0, 1, 17):
                    n = residue+period*t
                    check(n > 1 and n % 4 == 3, 'positive source and odd child')
                    m = (n-1)//2
                    check(literal(n, w) == literal(m, z), 'all-height family control')
                    check(sum(w) == sum(z)+1 and len(w) == len(z), 'half-source diagram')
                    selected = select_debt(n, bank)
                    check(selected is not None and selected['parent_word'] == w and
                          selected['child_word'] == z and selected['child'] == m,
                          'recognize exact supplied family instance')
                    if example is None:
                        example = selected
        rows.append(dict(s=s, e=e, source_word=u, child_word=v, cost=sum(u),
                         binary_c=c, odd_relative_density=str(Fraction(1, 2**(c-1))),
                         first_tested_example=example))
    for n in range(1, 20000, 2):
        ours, old = select_debt(n, inherited.SELECTED), inherited.select_join(n)
        check((ours is None) == (old is None), 'independent selector compatibility')
        if ours:
            check(ours['parent_word'] == old['parent_word'] and
                  ours['child_word'] == old['child_word'], 'exact incoming word compatibility')
    check(all(select_debt(n, bank) is None for n in (1, 7, 27, 703)),
          'extended debt bank retains hostile controls')
    return rows, counts, enumerated


def ternary_atlas():
    rows = [sibling_row(k) for k in range(21)]
    primitive, density = [], Fraction(0)
    for row in rows:
        covered = [old['k'] for old in primitive
                   if (row['k']-old['k']) % old['modulus'] == 0]
        row['covered_by'] = covered
        if not covered:
            primitive.append(row)
            density += Fraction(1, row['modulus'])
        for old in rows[:row['k']]:
            g = old['modulus']
            check(((row['residue']-old['residue']) % g == 0) ==
                  ((row['k']-old['k']) % g == 0), 'isometric nested/disjoint law')
    for h in range(1, 7):
        modulus = 3**h
        images = [-((4**k+2)//3)*pow(4**k, -1, modulus) % modulus
                  for k in range(modulus)]
        check(len(set(images)) == modulus, 'finite ternary address bijection')
    for k in range(32):
        for j in range(k):
            check(val(4**(k-j)-1, 3)-1 == val(k-j, 3), 'integer isometry/LTE control')
    tail = Fraction(27, 26*3**sibling_row(21)['ell'])
    observed = [n for n in range(1, 20000, 2) if sibling_child(n)]
    check(all(sibling_child(n) is None for n in (1, 7, 27, 703)), 'retained hostile sources')
    return dict(rows=rows, lower=str(density), upper=str(density+tail),
                tail_bound=str(tail), odd_sources_through_19999=len(observed),
                hostiles=[1, 7, 27, 703])


def crt_audit():
    cases = 0
    for b in range(1, 5):
        B = 2**b
        for t in range(1, 4):
            T = 3**t
            for a in range(1, B, 2):
                for r in range(T):
                    matches = [n for n in range(1, B*T, 2) if n % B == a and n % T == r]
                    check(len(matches) == 1, 'independent exact CRT count')
                    cases += 1
    return cases


def main():
    root = Path(__file__).resolve().parents[2]
    rows, counts, enumerated = find_debt_bank()
    removed = sum((Fraction(row['odd_relative_density']) for row in rows), Fraction(0))
    binary = Fraction(1, 4)-removed
    ternary = ternary_atlas()
    lower, upper = Fraction(ternary['lower']), Fraction(ternary['upper'])
    fusion = [binary*(1-upper), binary*(1-lower)]
    seeds = json.loads((root/'05-knowledge/results/checked_switch_phase19_20261004.json').read_text())['compiler'][-1]['seed_sources']
    incoming_hits = [n for n in seeds if inherited.select_join(n)]
    check(incoming_hits == [6783], 'incoming eight-family frozen-universe control')
    report = dict(status='PROVED scoped family and address mechanisms; FINITE-EXACT controls; Collatz OPEN',
                  debt_rows=rows, search_counts=counts, words_examined=enumerated,
                  binary_necessary_density=str(binary), ternary=ternary,
                  fused_necessary_density_interval=list(map(str, fusion)),
                  fused_decimal_interval=list(map(float, fusion)), crt_cases=crt_audit(),
                  incoming_eight_family_seed_hits=incoming_hits, checks=CHECKS)
    stem = root/'05-knowledge/results/collatz_binary_ternary_guard_fusion_20261004'
    stem.with_suffix('.json').write_text(json.dumps(report, indent=2)+'\n')
    lines = ['BINARY/TERNARY GUARD FUSION: universal coverage OPEN']
    lines += [json.dumps(row, sort_keys=True) for row in rows]
    lines += ['TERNARY_DENSITY '+str(float(lower))+' tail '+str(ternary['tail_bound']),
              'BINARY_REMAINDER '+str(binary), 'FUSED_INTERVAL '+str(list(map(float, fusion))),
              f'PASS: {CHECKS} checks; {enumerated} words; all costs below24 exhaustive.']
    out = '\n'.join(lines)+'\n'
    stem.with_suffix('.out').write_text(out)
    print(out, end='')


if __name__ == '__main__':
    main()
