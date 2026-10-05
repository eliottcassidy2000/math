"""Repair supplied common-future families by cutting the certified child path.

No API searches for a missing root certificate. The finite main supplies them
explicitly, and separates a small hostile census from the elementary proofs.
"""
from dataclasses import dataclass
from fractions import Fraction
import json

import inverse_ray_ternary_addresses_20261004 as codec


def need(ok, message):
    if not ok:
        raise ValueError(message)


def data(word):
    P = Q = 1
    B = 0
    for a in word:
        need(type(a) is int and a > 0, 'positive exact valuation')
        P, Q, B = 3*P, Q*2**a, 3*B+Q
    return P, Q, B


def replay(source, word):
    need(type(source) is int and source > 0 and source % 2, 'positive odd source')
    states = [source]
    for a in word:
        need(source != 1, 'strict first-hit prefix')
        value = 3*source+1
        actual = (value & -value).bit_length()-1
        need(type(a) is int and a == actual, 'actual valuation guard')
        source = value >> actual
        states.append(source)
    return tuple(states)


def certificate_route(certificate):
    codec.audit_certificate(certificate)
    source = codec.expand(certificate, bit_cap=codec.bit_bounds(certificate)[1])
    word = tuple(codec.exponent(node) for node in codec.chain(certificate))
    states = replay(source, word)
    need(states[-1] == 1, 'supplied rooted first-hit path')
    return states, word


def potential(certificate):
    _, word = certificate_route(certificate)
    P, Q, _ = data(word)
    return Fraction(P, Q)


@dataclass(frozen=True)
class Join:
    source: int
    source_word: tuple
    child_certificate: object
    child_steps: int

    def validate(self):
        need(type(self.child_steps) is int and self.child_steps >= 0, 'exact child prefix length')
        child_states, child_word = certificate_route(self.child_certificate)
        need(self.child_steps <= len(child_word), 'child prefix lies within the supplied certificate')
        source_states = replay(self.source, self.source_word)
        need(0 < child_states[0] < self.source, 'actual smaller child identity')
        need(source_states[-1] == child_states[self.child_steps], 'same checked endpoint')
        return source_states, child_states, child_word

    def family(self):
        _, child_states, child_word = self.validate()
        r, s = len(self.source_word), self.child_steps
        A, D = sum(self.source_word), sum(child_word[:s])
        P = 2**(A+1)*3**max(s-r, 0)
        Q = 2**(D+1)*3**max(r-s, 0)
        return dict(source=self.source, child=child_states[0], source_period=P,
                    child_period=Q, slope=Fraction(Q, P),
                    source_word=self.source_word, child_word=child_word[:s])


def child_suffix_cuts(join):
    """Retain every smaller certified suffix; return its own exact guard and rank."""
    _, child_states, child_word = join.validate()
    old = join.family()
    suffix = join.child_certificate
    output = []
    for cut in range(join.child_steps+1):
        child = child_states[cut]
        if child < join.source:
            new_join = Join(join.source, join.source_word, suffix, join.child_steps-cut)
            new = new_join.family()
            p, q, _ = data(child_word[:cut])
            multiplier = Fraction(p, q)
            factor = 3**min(cut, max(join.child_steps-len(join.source_word), 0))
            need(new['slope'] == old['slope']*multiplier, 'child-cut slope law')
            need(old['source_period'] == factor*new['source_period'], 'ternary source-guard relaxation')
            if child < child_states[0]:
                need(multiplier < 1 and new['slope'] < old['slope'], 'positive carry improves a descending cut')
            if old['slope'] > 1 and child*old['slope'] <= child_states[0]:
                need(new['slope'] < 1, 'sufficient child-side repayment threshold')
            output.append(dict(cut=cut, join=new_join, family=new, coefficient=multiplier,
                               source_guard_factor=factor, all_height=new['slope'] <= 1))
        if cut < join.child_steps:
            suffix = suffix.parent
    return tuple(output)


def attach_source_prefix(join):
    """Build the original source AST from the supplied child; no orbit search."""
    source_states, _, _ = join.validate()
    suffix = join.child_certificate
    for _ in range(join.child_steps):
        suffix = suffix.parent
    for index in range(len(join.source_word)-1, -1, -1):
        source, target = source_states[index:index+2]
        row = source % 3
        least = codec.kappa(target % 9, row)
        exponent = join.source_word[index]
        need(exponent >= least and (exponent-least) % 6 == 0, 'rooted inverse prefix')
        suffix = codec.extend(suffix, row, (exponent-least)//6)
    need(codec.expand(suffix, bit_cap=codec.bit_bounds(suffix)[1]) == join.source,
         'constructed certificate keeps the supplied source')
    return suffix


def audit_lifts(family, parameters):
    for t in parameters:
        n = family['source']+family['source_period']*t
        m = family['child']+family['child_period']*t
        left = replay(n, family['source_word'])
        right = replay(m, family['child_word'])
        need(left[-1] == right[-1], 'actual lifted common endpoint')
        if family['slope'] <= 1:
            need(0 < m < n, 'uniform original-source rank')


def explicit_route(n):
    # Used only to supply the finite premises / independent hostile controls.
    return certificate_route(codec.encode_source(n, step_cap=5000))


def main():
    print('SUPPLIED STATE JOIN EXTENSIONS: PROVED cut/potential laws; FINITE-EXACT declared controls')
    bad_child = codec.encode_source(231)
    child_states, child_word = certificate_route(bad_child)
    bad_join = Join(233, (), bad_child, child_states.index(233))
    bad = bad_join.family()
    cuts = child_suffix_cuts(bad_join)
    repair = next(row for row in cuts if row['family']['child'] == 31)
    good = repair['family']
    need(bad['slope'] == Fraction(2**27, 3**17) > 1, 'inherited expanding pair')
    need(good['slope'] == Fraction(8192, 59049), 'repaired child31 slope')
    need((good['source_period'], good['child_period']) == (118098, 16384), 'larger guarded family')
    need(repair['source_guard_factor'] == 2187, 'exact ternary precision released')
    audit_lifts(bad, (0, 1, 17))
    audit_lifts(good, (0, 1, 2, 17, 10**6))
    need(bad['child']+bad['child_period'] > bad['source']+bad['source_period'], 'original lift fails at first positive parameter')
    for t in (0, 1, 13):
        old_source = bad['source']+bad['source_period']*t
        new_source = good['source']+good['source_period']*(2187*t)
        need(old_source == new_source, 'repair preserves and enlarges source applicability')
    print('Expanding pair233/231:', str(bad['slope']), '; source period', bad['source_period'])
    print('Cut at child31:', str(good['slope']), '; family233+118098t <-31+16384t; guard factor2187')
    print('All retained smaller suffixes:', [(row['cut'], row['family']['child'], str(row['family']['slope'])) for row in cuts])

    source_cert = attach_source_prefix(bad_join)
    need(source_cert == codec.encode_source(233), 'point proof before family improvement')
    need(potential(source_cert)/potential(bad_child) == bad['slope'], 'root potential recovers the join slope')
    nstates, nword = certificate_route(source_cert)
    for further in range(len(nword)+1):
        continued = Join(233, nword[:further], bad_child, bad_join.child_steps+further)
        need(continued.family()['slope'] == bad['slope'], 'same pair cannot repair by moving down its common future')
    print('Same-pair continuation obstruction:', len(nword)+1, 'strict common endpoints retain expanding slope')

    states75, word75 = explicit_route(75)
    cert73 = codec.encode_source(73)
    states73, word73 = certificate_route(cert73)
    hostile = Join(75, word75, cert73, len(word73))
    partial = next(row for row in child_suffix_cuts(hostile) if row['cut'] == 5)
    need(partial['family']['child'] == 71 < 73 and partial['family']['slope'] > 1,
         'a smaller child improves but need not finish slope repayment')
    print('Partial-repair hostile75/73:cut to71 still has slope', str(partial['family']['slope']), '>1')

    cases = []
    for n, m, source_steps in ((7, 3, 4), (27, 15, 37), (703, 123, 51), (515, 257, 2)):
        nstates, nword = explicit_route(n)
        child_cert = codec.encode_source(m)
        mstates, mword = certificate_route(child_cert)
        endpoint = nstates[source_steps]
        child_steps = mstates.index(endpoint)
        join = Join(n, nword[:source_steps], child_cert, child_steps)
        family = join.family()
        result = attach_source_prefix(join)
        need(result == codec.encode_source(n), 'transported supplied original source')
        audit_lifts(family, (0, 1, 17))
        candidates = child_suffix_cuts(join)
        need(all(row['family']['source'] == n for row in candidates), 'candidate cuts do not replace the request')
        cases.append(dict(source=n, child=m, endpoint=endpoint, forward_steps=source_steps,
                          reverse_steps=child_steps, slope=str(family['slope']),
                          retained_cuts=len(candidates), source_period=family['source_period']))
    print('Inherited stubborn-source / grounded257 controls:', json.dumps(cases))
    p257, _ = explicit_route(257)
    need(p257[:5] == (257, 193, 145, 109, 41), '257 uses the formerly pending27 suffix')
    need(explicit_route(171)[0][:2] == (171, 257), 'smaller ancestor is not automatically grounded')
    print('257 scope:171->257, but257->193->145->109->41 needs the shared grounded suffix; no circular premise')

    # Focused hostile probe, not a discovery source or a global theorem.
    routes = {n: explicit_route(n) for n in range(1, 1002, 2)}
    indexed = {n: {x:i for i,x in enumerate(states)} for n,(states,_) in routes.items()}
    expanding = repaired = endpoint_below = endpoint_not_below = 0
    for n in range(3, 1002, 2):
        ns, nw = routes[n]
        for m in range(1, n, 2):
            ms, mw = routes[m]
            r = next(i for i,x in enumerate(ns) if x in indexed[m])
            s = indexed[m][ns[r]]
            lam = Fraction(2**sum(mw[:s])*3**r, 2**sum(nw[:r])*3**s)
            if lam <= 1:
                continue
            expanding += 1
            if ns[r] < n:
                endpoint_below += 1
            else:
                endpoint_not_below += 1
            candidates = [lam*Fraction(3**p, 2**sum(mw[:p])) for p,x in enumerate(ms[:s+1]) if x<n]
            if min(candidates) <= 1:
                repaired += 1
    need((expanding, repaired) == (2914, 2914), 'declared finite child-cut signal')
    need((endpoint_below, endpoint_not_below) == (2822, 92), 'distinguish direct-descent endpoints from nontrivial repairs')
    print('FINITE-EXACT pair probe:all odd1<=m<n<=1001;', expanding, 'expanding earliest-pair joins,', repaired, 'have a usable child cut')
    print('Of those:2822 endpoints are already below the source;92 endpoints are not, and require an earlier child cut')
    print('OPEN: no proof that every expanding pair admits such a cut; no source coverage inferred')

    rejected = 0
    for thunk in (lambda: Join(233, (), codec.ROOT, 0).validate(),
                  lambda: Join(233, (), bad_child, True).validate(),
                  lambda: Join(233, (), bad_child, 0).validate(),
                  lambda: replay(1, (2,)),
                  lambda: Join(27, (1, 2), codec.encode_source(13), 1).validate()):
        try:
            thunk()
        except ValueError:
            rejected += 1
        else:
            raise ValueError('hostile unexpectedly accepted')
    need(rejected == 5, 'strict root/source/endpoint guards')
    print('Rejected type, wrong endpoint and root-padding controls:', rejected)
    large_certificate = codec.extend(codec.ROOT, 2, 2000)
    large_states, large_word = certificate_route(large_certificate)
    need(len(large_word) == 1 and large_states[0].bit_length() > 10000 and large_states[-1] == 1,
         'valid supplied certificate is not rejected by an incidental default bit cap')
    print('Supplied one-edge certificate beyond10000bits:passed using its explicit codec bit bound')
    print('PASS: APIs consume supplied certificates; symbolic family conclusions remain conditional on their child premises')


if __name__ == '__main__':
    main()
