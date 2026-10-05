"""Source-paid sibling reroutes at arbitrary checked forward checkpoints.

The search is complete at each pre-descent checkpoint. The forward horizon is
explicit and finite. Exact guarded families remain smaller-child obligations.
"""
from fractions import Fraction
import json

import collatz_early_reroute_20261004 as early
import collatz_complement_routing_20261004 as prior

CHECKS = 0


def need(ok, why):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(why)


def stats(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'exact tuple of positive valuations')
    A, B = 0, 0
    for a in word:
        B = 3*B+2**A
        A += a
    return len(word), A, B


def phase(word):
    r, A, B = stats(word)
    modulus = 3**(r+1)
    target = -2**(A+1)*pow(3*B+2**A, -1, modulus) % modulus
    kappa = 0
    for level in range(1, r+1):
        m, period = 3**(level+1), 3**(level-1)
        choices = [kappa+j*period for j in range(3)
                   if pow(4, kappa+j*period, m) == target % m]
        need(len(choices) == 1, 'unique principal-unit branch phase')
        kappa = choices[0]
    return kappa, 3**r


def guard_data(word, ell, k):
    r, A, B = stats(word)
    need(type(ell) is int and ell >= 1 and type(k) is int and k >= 1,
         'positive exact inverse depth and sibling index')
    H = 4**k*(3*B+2**A)+2**(A+1)
    if ell <= r:
        if H % 3**(ell+1):
            return None
        ternary, residue = 1, 0
    else:
        if H % 3**(r+1):
            return None
        ternary = 3**(ell-r)
        residue = -(H//3**(r+1))*pow(4**k, -1, ternary) % ternary
    rho = Fraction(2**(ell+2*k-1)*3**r, 2**A*3**ell)
    beta = Fraction(2**(ell-1)*H, 2**A*3**(ell+1))
    return dict(r=r, A=A, B=B, H=H, modulus=ternary, residue=residue,
                rho=rho, beta=beta)


def compile_family(word, ell, k, last):
    """Compile fixed exact source/child words and their all-height paid ray."""
    need(type(last) is int and last >= 1, 'exact terminal source valuation')
    data = guard_data(word, ell, k)
    if data is None or data['rho'] >= 1:
        return None
    source_word = word+(last,)
    r, A, B = stats(source_word)
    binary = 2**(A+1)
    b = (2**A-B)*pow(3**r, -1, binary) % binary
    M, residue = data['modulus'], data['residue']
    n = residue+M*((b-residue)*pow(M, -1, binary) % binary)
    period = binary*M
    threshold = (data['beta']-1)/(1-data['rho'])
    t = max(0, (threshold-n)//period+1, (2-n)//period+1)
    n += period*t
    child = data['rho']*n+data['beta']-1
    increment = data['rho']*period
    need(child.denominator == increment.denominator == 1 and 0 < child < n,
         'guarded integral all-height source payment')
    child_word = (1,)*(ell-1)+(2, last+2*k-2)
    def arithmetic_chain(value, letters):
        states = [value]
        for letter in letters:
            value, actual = early.step(value)
            need(actual == letter, 'compiled final-oddness cylinder has exact valuations')
            states.append(value)
        return states
    left = arithmetic_chain(n, source_word)
    right = arithmetic_chain(int(child), child_word)
    trimmed_root_member = 1 in left[:-1] or 1 in right[:-1]
    if trimmed_root_member:
        # Every fixed-word intermediate increases strictly with the family
        # parameter. One increment moves every possible root prefix above1.
        n += period
        child += increment
    need(early.literal(n, source_word)[-1] == early.literal(int(child), child_word)[-1],
         'first-hit-safe least family member after finite root exception')
    all_growth = all(3**j > 2**sum(source_word[:j]) for j in range(1, len(source_word)+1))
    return dict(source=n, period=period, child=int(child), child_period=int(increment),
                source_word=source_word, child_word=child_word, ell=ell, k=k,
                rho=data['rho'], intercept=data['beta']-1,
                trimmed_root_member=trimmed_root_member,
                all_source_prefixes_grow=all_growth)


def cut_family_hit(n):
    """The widened later quarter-child switch, preserving the supplied source."""
    early.odd(n)
    if n % 2048 != 155:
        return None
    word = (1, 2, 1, 1, 1, 2)
    x = early.literal(n, word)[-1]
    h = (729*n+669)//1024
    target, a = early.step(x)
    need(x == 4*h+1 and h % 2 and 0 < h < n and a >= 3,
         'actual quarter-child identity with original-source payment')
    need(early.literal(h, (a-2,))[-1] == target, 'variable-terminal common future')
    return dict(source=n, child=h, source_word=word+(a,), child_word=(a-2,), endpoint=target)


def scan(n, horizon=128):
    """Retain every paid reroute strictly before the first ordinary descent."""
    early.odd(n)
    need(type(horizon) is int and horizon >= 0, 'nonnegative exact observation horizon')
    if n == 1:
        return dict(status='ROOT', steps=0, word=(), hits=[], guard_tests=0, frontier=1)
    cap = n.bit_length()-1
    x, word, hits, tests = n, (), [], 0
    for count in range(horizon):
        y, a = early.step(x)
        if y < n:
            return dict(status='ORDINARY_DESCENT', steps=count+1, word=word+(a,),
                        hits=hits, guard_tests=tests, frontier=y)
        for ell in range(2, cap+1):
            for k in range(1, ell//2+1):
                if ell < early.clock(k):
                    continue
                tests += 1
                raw = 4**k*x+(4**k+2)//3
                actual = raw % 3**ell == 0
                data = guard_data(word, ell, k)
                compiled = data is not None and n % data['modulus'] == data['residue']
                need(actual == compiled, 'direct checkpoint and prefix-carry guards agree')
                if not actual:
                    continue
                child = 2**(ell-1)*(raw//3**ell)-1
                if 0 < child < n:
                    need(data['rho']*n+data['beta']-1 == child, 'independent rational child map')
                    v = (1,)*(ell-1)+(2, a+2*k-2)
                    need(early.literal(child, v)[-1] == y, 'actual checked child join')
                    hits.append(dict(checkpoint=count, source_depth=count+1, ell=ell, k=k,
                                     child=child, source_word=word+(a,), child_word=v))
        word += (a,)
        x = y
    return dict(status='PENDING', steps=horizon, word=word, hits=hits,
                guard_tests=tests, frontier=x)


def exponent3(target, bits):
    need(type(target) is int and target % 8 in (1, 3) and type(bits) is int and bits >= 3,
         'target in the powers-of-three subgroup')
    exponent = 0 if target % 8 == 1 else 1
    for level in range(4, bits+1):
        period = 2**(level-3)
        choices = [exponent+j*period for j in range(2)
                   if pow(3, exponent+j*period, 2**level) == target % 2**level]
        need(len(choices) == 1, 'unique binary exponent digit')
        exponent = choices[0]
    return exponent, 2**(bits-2)


def main():
    print('ARBITRARY CHECKPOINT REROUTES: finite complete local search; universal coverage OPEN')
    w = (1, 2, 1, 1, 1, 2, 3, 1, 1, 2)
    family = compile_family(w, 4, 1, 1)
    need((family['source'], family['period'], family['child'], family['child_period']) ==
         (155, 131072, 111, 93312), 'new exact fixed-word family')
    need(family['rho'] == Fraction(729, 1024) and family['intercept'] == Fraction(669, 1024),
         'ordered source-relative affine child map')
    need(family['all_source_prefixes_grow'], 'all eleven source prefixes grow at every height')
    r, A, B = stats(w)
    data = guard_data(w, 4, 1)
    kappa, modulus = phase(w)
    need((r, A, B) == (10, 15, 120749) and data is not None
         and data['H'] % 3**5 == 0 and data['H'] % 3**6 != 0,
         'strictly shorter inverse-depth guard')
    need((1-kappa) % modulus != 0, 'full-depth phase would falsely reject this valid family')
    bank = prior.load_bank()
    for t in (0, 1, 2, 17, 10**6):
        n = family['source']+family['period']*t
        h = family['child']+family['child_period']*t
        left = early.literal(n, family['source_word'])
        right = early.literal(h, family['child_word'])
        need(left[-1] == right[-1] and all(x > n for x in left[1:]) and h < n,
             'independent all-height diagram and immutable-source comparisons')
        need(prior.fusion.select_debt(n, bank) is None, 'inherited binary16 excluded by fixed prefix')
    print('NEW_FAMILY', json.dumps({k: list(v) if isinstance(v, tuple) else str(v) if isinstance(v, Fraction)
                                   else v for k, v in family.items()}, sort_keys=True))
    print('SHORT_GUARD r10,ell4,H:', data['H'], ';full phase kappa:', kappa, 'mod', modulus)
    epsilon = Fraction(1, 3**9)/(1-Fraction(1, 3**10))
    older = Fraction(1, 3**31)/(1-Fraction(1, 3**30))
    residual = 1-Fraction(4458180914749392, 10**16)-epsilon-older
    need(residual > Fraction(55413, 100000), 'strict new relative density after all named banks')
    print('New relative density >0.55413 after originbank,allposition1rows,andpriorcomposedrows.')
    print('Exact conservative relative lower bound:', residual,
          ';odd-relative lower bound:', residual/65536)

    for t in range(128):
        n = 155+2048*t
        hit = cut_family_hit(n)
        need(hit is not None and hit['child'] == 111+1458*t, 'minimal child-odd dyadic family')
        values = early.literal(n, hit['source_word'])
        need(all(x > n for x in values[1:7]), 'all six retained source prefixes grow')
        if t % 2 == 0:
            need(hit['source_word'][-1] == 3 and hit['endpoint'] > n,
                 'even parameter is the genuine pre-descent bypass')
        else:
            need(hit['source_word'][-1] >= 4 and hit['endpoint'] < n,
                 'odd parameter already first-descends at the final edge')
        need(prior.fusion.select_debt(n, bank) is None, 'widened binary16 exclusion')
        need(early.energy.rank(hit['child']) < early.energy.rank(n), 'strict original graph rank')
    print('CUT_FAMILY source155+2048t,child111+1458t;variable last a>=3,child a-2.')
    print('Even t:source155+4096u,child111+2916u;all7source steps grow.')
    print('Odd t:already first-descends at7;minimal union guard is64x broader than discovery.')
    print('Nontrivial even-row odd-relative new density lower bound:', residual/2048)

    exponent, period = exponent3(155, 17)
    need((exponent, period) == (27107, 32768) and exponent % 8 == 3,
         'new subfamily inside the previous power-of-three survivor')
    for t in (0, 1, 17):
        need(pow(3, exponent+period*t, 131072) == 155, 'all-height exponent source guard')
    n = 3**exponent
    h = (729*n+669)//1024
    left = early.literal(n, family['source_word'])
    right = early.literal(h, family['child_word'])
    need(left[-1] == right[-1] and all(x > n for x in left[1:]) and 0 < h < n,
         'literal large power-of-three bypass before source descent')
    need(early.height_obstruction(n), 'new checkpoint leaves the whole early grammar behind')
    print('POWER3 exponent a=27107+32768t;least source bits:', n.bit_length(),
          ';11 source and5 child literal edges; child remains an obligation.')
    need(exponent3(155, 11) == (483, 512) and exponent3(155, 12) == (483, 1024),
         'first-common-future and terminal-guard widening of exponent family')
    n = 3**483
    hit = cut_family_hit(n)
    need(hit is not None and hit['source_word'][-1] == 3 and hit['endpoint'] > n
         and early.height_obstruction(n), 'literal widened power-of-three bypass')
    print('WIDENED_POWER3 a=483+512t;nontrivial growing-join half a=483+1024u;least bits:', n.bit_length())

    targets = [('3^'+str(a), 3**a) for a in range(3, 68, 8)]+[('7', 7), ('703', 703)]
    report, total_steps, total_tests = [], 0, 0
    for label, source in targets:
        result = scan(source)
        need(result['status'] == 'ORDINARY_DESCENT' and not result['hits'],
             'targeted no-earlier-reroute result, not an all-source claim')
        report.append((label, result['steps'], result['guard_tests']))
        total_steps += result['steps']
        total_tests += result['guard_tests']
    print('TARGETED [source,firstdescent,guardtests]:', report)
    print('Targeted work:', total_steps, 'literal odd edges;', total_tests,
          'candidate guards;11 sources,horizon128;no new earlier hit in that fixed target list.')
    small = scan(155)
    need(any(hit['checkpoint'] == 10 and hit['ell'] == 4 and hit['k'] == 1 for hit in small['hits']),
         'complete checkpoint search discovers the new later diagram')
    truncated = scan(703, 8)
    need(truncated['status'] == 'PENDING' and truncated['steps'] == 8
         and len(truncated['word']) == 8 and truncated['frontier'] > 703,
         'budget exhaustion retains the checked prefix and actual frontier')
    previous_family = prior.family(10)
    recovered = compile_family((1, 2), 33, 10, 1)
    need(recovered is not None and recovered['rho'] == previous_family.slope,
         'generic compiler recovers inherited after-reset slope')
    need((previous_family.source-recovered['source']) % recovered['period'] == 0,
         'inherited narrower four-letter guard lies in compiled cell')
    print('Positive controls:155 later checkpoint;prior(1,2),ell33,k10 family;703 cap8 staysPENDING.')

    phase_controls = 0
    for word in ((), (1,), (1, 2), w):
        r, _, _ = stats(word)
        kappa, modulus = phase(word)
        for k in range(1, 25):
            data = guard_data(word, r+1, k)
            need((data is not None) == ((k-kappa) % modulus == 0), 'full-depth phase iff')
            phase_controls += 1
    need(phase((1, 2)) == (1, 9), 'previous after-reset address recovered')
    root_boundary = compile_family((4,), 1, 1, 2)
    need(root_boundary['trimmed_root_member'] and root_boundary['source'] == 133
         and root_boundary['child'] == 33, 'padded5/1 member removed before first-hit export')
    print('Independent full-phase controls:', phase_controls, ';ell<r retained as separate truncated guard.')
    rejected = 0
    for thunk in (lambda: scan(True), lambda: scan(7, -1), lambda: stats((1, True)),
                  lambda: guard_data((1,), 0, 1), lambda: guard_data((1,), 2, False),
                  lambda: cut_family_hit(156)):
        try:
            thunk()
        except ValueError:
            rejected += 1
        else:
            raise ValueError('invalid exact domain accepted')
    print('Guard/type rejections:', rejected)
    print('PASS:', CHECKS, 'always-active compiler checks; no supplied or oracle child used by search.')


if __name__ == '__main__':
    main()
