"""Exact exhausted-H to expanding-12 switches with immutable-source payment.

Standalone reproduction; no import-time experiment or file writes. Existing
kernel/codec supply the inherited identities. No orbit-derived proof is hidden
inside evaluate(), certificate(), or the all-height completed-plan constructor.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
import collatz_recursive_dependency_kernel_20261004 as kernel

H = (729, 1024, 669)
G = (9, 8, 5)


def need(condition, message):
    if not condition:
        raise ValueError(message)


def positive(value, message):
    need(type(value) is int and value >= 1, message)


def rank(n):
    kernel.odd(n)
    if n == 1:
        return 0, 0
    k = kernel.v2(n-1)
    return 3**k*((n-1) >> k)**2, k


@dataclass(frozen=True)
class Switch:
    source: int
    h_repeats: int
    g_blocks: int


def metadata(m, q):
    positive(m, 'positive H count')
    positive(q, 'positive G count')
    first = kernel.power(H, m)[0]
    second = kernel.power(G, q)[0]
    return kernel.compose(first, second)


def source_word(m, q):
    positive(m, 'positive H count')
    positive(q, 'positive G count')
    return kernel.V+kernel.W*(m-1)+(3, 2)+(1, 2)*(q-1)


def evaluate(record):
    """Prefix guard and original-source payment; no terminal home inference."""
    need(type(record) is Switch, 'typed switch record')
    kernel.odd(record.source)
    positive(record.h_repeats, 'positive exact H count')
    positive(record.g_blocks, 'positive exact G count')
    n, m, q = record.source, record.h_repeats, record.g_blocks
    c = kernel.repeat_endpoint(n, m)
    need(kernel.v2(c+5) >= 3*q+1, 'retained actual-12 repeat guard')
    y = (9**q*(c+5) >> (3*q))-5
    need(y > c >= 111 and y % 2 == 1, 'actual G growth and positive endpoint')
    p, denominator, carry = metadata(m, q)
    need(p*n+carry == denominator*y, 'immutable-source affine summary')
    paid = y < n
    need(paid == ((denominator-p)*n > carry), 'exact source-relative payment iff')
    if q <= 2*m:
        need(paid, 'uniform two-block-per-H payment theorem')
    if q >= 3*m:
        need(not paid, 'three-block-per-H coefficient hostile')
    if paid:
        need(rank(y) < rank(n), 'original-source graph rank paid')
    return dict(terminal=c, endpoint=y, paid=paid, metadata=(p, denominator, carry),
                rank_paid=rank(y) < rank(n),
                maximal_h=(kernel.fuel(n)-1)//10 == m,
                maximal_g=(kernel.v2(c+5)-1)//3 == q)


def family(m, q, residual=1):
    """All-height exact exhaustion; residual=1,2,3 partitions the G exit fuel."""
    positive(m, 'positive H count')
    positive(q, 'positive G count')
    need(q >= 2, 'exhaustion family uses at least two G blocks')
    need(type(residual) is int and residual in (1, 2, 3), 'three residual fuel classes')
    p = 729**m
    shift = 3*q+residual
    t0 = 2144*pow(295*(1 << shift), -1, p) % p
    if t0 % 2 == 0:
        t0 += p
    need(t0 > 0 and t0 % 2 == 1, 'least positive odd ternary address')
    c0 = (t0 << shift)-5
    delta = 295*c0-669
    need(delta % p == 0, 'inverse-H terminal ternary guard')
    numerator = 669+1024**m*(delta//p)
    need(numerator % 295 == 0, 'integer source from inverse-H expression')
    n0 = numerator//295
    y0 = (9**q*t0 << residual)-5
    row = dict(m=m, q=q, residual=residual, t0=t0, t_period=2*p,
               source=n0, source_period=1 << (10*m+shift+1),
               terminal=c0, terminal_period=(1 << (shift+1))*p,
               endpoint=y0, endpoint_period=(1 << (residual+1))*3**(6*m+2*q))
    state = evaluate(Switch(n0, m, q))
    need(state['terminal'] == c0 and state['endpoint'] == y0,
         'least member identities')
    need(state['maximal_h'] and state['maximal_g'], 'both patterns exactly exhausted')
    return row


def member(m, q, residual=1, index=0):
    need(type(index) is int and index >= 0, 'nonnegative family parameter')
    row = family(m, q, residual)
    return Switch(row['source']+index*row['source_period'], m, q)


def certificate(record, endpoint_word):
    """Supplied exact endpoint proof -> original-source ROOT certificate."""
    state = evaluate(record)
    need(state['paid'], 'this adapter requires a paid switch')
    kernel.literal(state['endpoint'], endpoint_word, require_root=True)
    word = source_word(record.h_repeats, record.g_blocks)+endpoint_word
    kernel.literal(record.source, word, require_root=True)
    cert = kernel.codec_from_word(word)
    need(kernel.codec.expand(cert, bit_cap=kernel.codec.bit_bounds(cert)[1]) == record.source,
         'ROOT export preserves the supplied source')
    return cert


@dataclass(frozen=True)
class CompletedPlan:
    repeats: int
    root_exponent: int


def audit_plan(plan):
    need(type(plan) is CompletedPlan, 'typed completed switch plan')
    positive(plan.repeats, 'positive completed depth')
    need(type(plan.root_exponent) is int and plan.root_exponent >= 4
         and plan.root_exponent % 2 == 0, 'positive-height even root exponent')
    m = plan.repeats
    t0 = family(m, 2*m)['t0']
    modulus = 3**(10*m+1)
    target = (2*3**(4*m+1)*t0-14) % modulus
    need(pow(2, plan.root_exponent, modulus) == target, 'full completed ternary phase')
    return True


def completed_plan(m, lift=0):
    positive(m, 'positive completed depth')
    need(type(lift) is int and lift >= 0, 'nonnegative completed-family parameter')
    row = family(m, 2*m)
    depth = 10*m+1
    target = (2*3**(4*m+1)*row['t0']-14) % 3**depth
    s, period = kernel.principal_log4(target, depth)
    if s < 2:
        s += period
    plan = CompletedPlan(m, 2*(s+lift*period))
    audit_plan(plan)
    return plan


def materialize(plan, bit_cap=1000000):
    audit_plan(plan)
    positive(bit_cap, 'positive expansion cap')
    m, b = plan.repeats, plan.root_exponent
    need(b+16*m+2 <= bit_cap, 'declared expansion cap')
    denominator = 2*3**(4*m+1)
    numerator = (1 << b)+14
    need(numerator % denominator == 0, 'integer maximal-run parameter')
    t = numerator//denominator
    row = family(m, 2*m)
    need(t >= row['t0'] and (t-row['t0']) % row['t_period'] == 0,
         'completed plan lies on its exact exhaustion row')
    record = member(m, 2*m, index=(t-row['t0'])//row['t_period'])
    state = evaluate(record)
    need(3*state['endpoint']+1 == 1 << b, 'explicit one-edge ROOT exit')
    return record


def plan_source_mod(plan, modulus):
    audit_plan(plan)
    positive(modulus, 'positive modular precision')
    m, b = plan.repeats, plan.root_exponent
    p, q, carry = metadata(m, 2*m)
    denominator = 3*p
    enlarged = denominator*modulus
    numerator = (q*(pow(2, b, enlarged)-1)-3*carry) % enlarged
    need(numerator % denominator == 0, 'closed-form modular source division')
    return numerator//denominator


def main():
    need(Fraction(104781, 6487) < 17, 'one-cell domination threshold')
    need(metadata(1, 2) == (59049, 65536, 141229), 'two-credit composition')
    need(metadata(1, 3) == (531441, 524288, 1598741), 'three-credit hostile')
    positive_cases = hostile_cases = 0
    for m in range(1, 13):
        for q in (2*m, 3*m):
            word = source_word(m, q)
            need(kernel.stats(word) == metadata(m, q), 'actual versus formal ordered carry')
            for residual in (1, 2, 3):
                row = family(m, q, residual)
                for index in (0, 1, 7, 31):
                    record = member(m, q, residual, index)
                    state = evaluate(record)
                    need(state['maximal_h'] and state['maximal_g'], 'both maximal repeat counts')
                    need(state['terminal'] == row['terminal']+index*row['terminal_period'],
                         'all-height terminal affine parameter')
                    need(state['endpoint'] == row['endpoint']+index*row['endpoint_period'],
                         'all-height endpoint affine parameter')
                    need(kernel.fuel(record.source) == 10*m+5, 'exact H fuel at the switch')
                    need(kernel.v2(state['terminal']+5) == 3*q+residual, 'exact new-pattern fuel')
                    need(kernel.literal(record.source, word) == state['endpoint'],
                         'independent actual source-prefix replay')
                    need(kernel.literal(state['terminal'], (1, 2)*q) == state['endpoint'],
                         'independent actual switched replay')
                    if state['paid']:
                        positive_cases += 1
                    else:
                        hostile_cases += 1
                        if residual in (2, 3):
                            need(record.source % 4 == state['endpoint'] % 4 == 3
                                 and rank(state['endpoint']) > rank(record.source),
                                 'size increase also fails graph-rank payment on these two strata')
    print('Uniform payment: H^m then any legal G^q with q<=2m pays original source; G=actual12.')
    print('Switch universe: m=1..12, q=2m/3m, residual fuel1/2/3, row parameters0/1/7/31.')
    print('Literal exact exhausted-family controls:', positive_cases,
          'pay integer size and rank;', hostile_cases, 'fail integer size payment.')
    print('The96 designated residual2/3 hostile controls also strictly increase graph rank.')
    for q in (2, 3):
        record = member(1, q)
        state = evaluate(record)
        print('Sharp control:', record, 'terminal', state['terminal'],
              'endpoint', state['endpoint'], 'paid', state['paid'])
    rank_hostile = member(1, 3, 2)
    state = evaluate(rank_hostile)
    print('Three-credit graph-rank hostile:', rank_hostile.source, 'via', state['terminal'],
          'to', state['endpoint'], 'ranks', rank(rank_hostile.source), rank(state['endpoint']))
    observations = 0
    for m in range(1, 5):
        record = member(m, 2*m)
        state = evaluate(record)
        supplied = kernel.observed_word(state['endpoint'])
        observations += len(supplied)
        cert = certificate(record, supplied)
        need(kernel.codec.ranks(cert) == (10*m+len(supplied),
                                          26*m+sum(supplied)+len(supplied)),
             'first-hit splice ranks')
    print('Four independently supplied endpoint proofs: observed odd edges', observations,
          '; these observations are outside the switch verifier.')
    plan = completed_plan(1)
    record = materialize(plan)
    cert = certificate(record, (plan.root_exponent,))
    print('Constructed completed control: m1, root exponent', plan.root_exponent,
          'source bits', record.source.bit_length(), 'ranks', kernel.codec.ranks(cert))
    modular = 0
    for m in range(1, 13):
        for lift in (0, 1, 7):
            plan = completed_plan(m, lift)
            word = source_word(m, 2*m)+(plan.root_exponent,)
            row = family(m, 2*m)
            need(plan_source_mod(plan, row['source_period']) == row['source'],
                 'symbolic completed source lies on maximal-exhaustion row')
            for modulus in (2**20, 3**8, 19**3, 295):
                need(plan_source_mod(plan, modulus) == kernel.inverse_word_mod(word, modulus),
                     'independent closed-form and inverse-letter symbolic readers')
                modular += 1
            cert = kernel.codec_from_word(word)
            need(kernel.codec.ranks(cert) == (10*m+1, 26*m+plan.root_exponent+1),
                 'symbolic completed first-hit ranks')
    print('Completed plans m1..12, lifts0/1/7:', modular,
          'independent modular comparisons; only m1 source expanded.')
    good = member(1, 2)
    hostiles = (lambda: evaluate(replace(good, source=good.source+2)),
                lambda: evaluate(replace(good, h_repeats=True)),
                lambda: evaluate(replace(good, source=float(good.source))),
                lambda: member(1, 2, residual=True), lambda: member(1, 2, index=-1),
                lambda: certificate(good, (2,)),
                lambda: audit_plan(CompletedPlan(1, 964)),
                lambda: materialize(completed_plan(2), bit_cap=1000000))
    for hostile in hostiles:
        kernel.rejected(hostile)
    print('Type/source/terminal/phase/expansion-cap hostiles rejected:', len(hostiles))
    print('Scope: conditional payment and constructed ROOT families; no new source coverage beyond H.')
    print('Local expansion is permitted only against retained original-source debt; exhaustion alone pays nothing.')


if __name__ == '__main__':
    main()
