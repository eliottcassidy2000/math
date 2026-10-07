"""Signed sibling ladders with marked affine states and exact tail guards.

The existing225 head bank is read only. Production exports actual paid
receipts or symbolic exponent phases, never discovers a ROOT certificate.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
from collections import Counter

import collatz_twoanchor_head_decoder_20261007b as old
import collatz_completion_anchor_20261007b as base
import collatz_uncovered_join_routes_20261007 as routes


@dataclass(frozen=True)
class Relation:
    debt: int
    offset: int


PLUS_ONE = Relation(3, -26)
TWO_TWOS = Relation(4, 10)


def mark(relation):
    if type(relation) is not Relation:
        raise ValueError("exact marked Relation required")
    old.integer(relation.debt, 1)
    old.integer(relation.offset)
    if relation.offset % 2:
        raise ValueError("offset must preserve odd source type")
    return relation


@dataclass(frozen=True)
class Ladder:
    relation: Relation
    source: tuple
    partner: tuple
    gap: int


def compile_ladder(source, gap, relation=PLUS_ONE):
    old.word_type(source)
    old.integer(gap)
    mark(relation)
    P, Q, B = old.carrier(source)
    L, C = len(source)+relation.debt, sum(source)-2*gap
    if C < L:
        return None
    if gap >= 0:
        target = B+P*relation.offset+(1 << C)*((4**gap-1)//3)
    else:
        target = B+P*relation.offset-Q*((4**(-gap)-1)//3)
    partner = old.decode(L, C, target)
    return None if partner is None else Ladder(relation, source, partner, gap)


def audit(packet):
    if type(packet) is not Ladder:
        raise ValueError("exact Ladder required")
    mark(packet.relation)
    old.word_type(packet.source)
    old.word_type(packet.partner)
    old.integer(packet.gap)
    if compile_ladder(packet.source, packet.gap, packet.relation) != packet:
        raise ValueError("head packet changed its marked affine identity")
    return packet


def all_ladders(source, relation=PLUS_ONE):
    old.word_type(source)
    mark(relation)
    P, Q, B = old.carrier(source)
    limit = (sum(source)-len(source)-relation.debt)//2
    gaps = list(range(1, limit+1))
    i = 0
    while B+P*relation.offset-Q*((4**i-1)//3) >= 1:
        gaps.append(-i)
        i += 1
    return tuple(p for r in sorted(gaps)
                 if (p := compile_ladder(source, r, relation)) is not None)


def reduce_zero(packet):
    """Strip shared final letters, then expose the nonzero signed ladder."""
    audit(packet)
    if packet.gap != 0:
        raise ValueError("zero-gap endpoint equality required")
    u, v, suffix = packet.source, packet.partner, ()
    while u and u[-1] == v[-1]:
        suffix = (u[-1],)+suffix
        u, v = u[:-1], v[:-1]
    if not u or not v or (v[-1]-u[-1]) % 2:
        raise ArithmeticError("zero-gap normal form violated its parity invariant")
    a, b = u[-1], v[-1]
    core = compile_ladder(u[:-1], (b-a)//2, packet.relation)
    if core is None or core.gap == 0 or core.partner != v[:-1]:
        raise ArithmeticError("zero-gap stripping lost the full carrier")
    return core, a, suffix


def coarse_cell(packet):
    """Formal valuation head/tail guard; strict exports also exclude ROOT padding."""
    audit(packet)
    P, Q, B = old.carrier(packet.source)
    extra = max(0, -2*packet.gap)
    effective = sum(packet.source)+extra
    modulus = 1 << (effective+1)
    residue = (-(3*B+Q)*pow(3*P, -1, modulus)) % modulus
    return residue, modulus, effective


def native_receipt(packet, parameter):
    audit(packet)
    old.integer(parameter, 0)
    c = 1+max(0, -2*packet.gap)
    left = packet.source+(c,)
    right = packet.partner+(c+2*packet.gap,)
    residue, modulus = old.native_cell(right)
    factor, d = 3**packet.relation.debt, packet.relation.offset
    least = max(3, (-d)//(factor-1)+1)
    first = max(0, (least-residue+modulus-1)//modulus)
    y = residue+modulus*(first+parameter)
    x = factor*y+d
    P, Q, B = old.carrier(right)
    endpoint = (P*y+B)//Q
    return routes.audit(routes.Receipt(x, y, left, right, endpoint))


def discharge(receipt, supplied_child_word):
    """Authenticate a supplied first-hit child word; never search its future."""
    return routes.discharge(receipt, supplied_child_word)


def state_type(packet, twos):
    audit(packet)
    old.integer(twos, 2)
    if not packet.source or packet.source[0] == 2:
        raise ValueError("nonempty head after the maximal two-run required")
    if not ((packet.relation == PLUS_ONE and twos >= 3)
            or (packet.relation == TWO_TWOS and twos == 2)):
        raise ValueError("clearing state is not authenticated at this two-run length")


@dataclass(frozen=True)
class Phase:
    ladder: Ladder
    twos: int
    residue: int
    period: int
    first_parameter: int


def phase(packet, twos):
    state_type(packet, twos)
    x, _, A = coarse_cell(packet)
    bits = 2*twos+A
    modulus = 1 << bits
    target = (1+(1 << (2*twos-1))*(x-1)*pow(3, -twos, modulus)) % modulus
    e, period = base.log_three(target, bits)
    minimum = 2*A+twos+1
    first = max(0, (minimum-(e+1)+period-1)//period)
    return Phase(packet, twos, e+1, period, first)


def audit_phase(record):
    if type(record) is not Phase:
        raise ValueError("exact Phase required")
    for x in (record.twos, record.residue, record.period, record.first_parameter):
        old.integer(x, 0)
    if phase(record.ladder, record.twos) != record:
        raise ValueError("phase lost its source guard or height cut")
    return record


def native_guard(packet, twos, K):
    state_type(packet, twos)
    old.integer(K, 1)
    x, M, A = coarse_cell(packet)
    divisor = 1 << (2*twos-1)
    modulus = divisor*M
    num = (pow(3, twos, modulus)*(pow(3, K-1, modulus)-1)) % modulus
    return num % divisor == 0 and (1+num//divisor) % M == x


def residual_parameter(packet, K, twos):
    state_type(packet, twos)
    old.integer(K, 5)
    x, _, A = coarse_cell(packet)
    modulus = 1 << (2*twos+A)
    target = (1+(1 << (2*twos-1))*(x-1)*pow(3, -twos, modulus)) % modulus
    return target*pow(3, -(K-1), modulus) % modulus, modulus


def apply_residual(packet, source, K, twos, deletion=4):
    state_type(packet, twos)
    routes.odd(source)
    old.integer(K, 5)
    if type(deletion) is not int or deletion not in (3,4):
        raise ValueError("clearing deletion must be3 or4")
    if base.v2(source+1) != K:
        raise ValueError("source changed its exact initial ones-run")
    _, _, A = coarse_cell(packet)
    if K-1 < 2*A+twos:
        raise ValueError("retain the sufficient immutable-source growth cut")
    left_prefix = (1,)*(K-1)+(2,)*twos+packet.source
    Z = routes.replay(source, left_prefix)[-1]
    if Z <= source:
        raise ValueError("source head has not retained proved growth")
    c = base.v2(3*Z+1)
    if c+2*packet.gap < 1:
        raise ValueError("source-side ladder lacks its terminal valuation reserve")
    child = ((source+1) >> deletion)-1
    if twos == 2:
        clearing = (4,) if deletion == 3 else (2,2)
        tail_twos = ()
    else:
        clearing = (4,1,1) if deletion == 3 else (2,2,1,1)
        tail_twos = (2,)*(twos-3)
    right = ((1,)*(K-deletion-1)+clearing+tail_twos
             +packet.partner+(c+2*packet.gap,))
    return routes.audit(routes.Receipt(source, child, left_prefix+(c,), right,
                                      (3*Z+1) >> c))


def main():
    checks = 0
    def check(ok):
        nonlocal checks
        checks += 1
        if not ok: raise ArithmeticError("exact control failed")
    def rejects(fn):
        try: fn()
        except (ValueError, TypeError): check(True)
        else: check(False)

    # Change the rule grammar on a fixed bounded head universe; no enlarged census.
    packets = []
    for ell in range(1,4):
        for u in product(range(1,13), repeat=ell):
            if u[0] == 2: continue
            rows = all_ladders(u)
            packets.extend(rows)
            P,Q,B = old.carrier(u)
            for p in rows:
                vP,vQ,vB = old.carrier(p.partner)
                scale = F(4)**p.gap
                check(F(vP+vB,vQ) == scale*F(P+B,Q)+(scale-1)/3)
                if p.gap == 0:
                    core,a,suffix = reduce_zero(p)
                    check(p.source == core.source+(a,)+suffix)
                    check(p.partner == core.partner+(a+2*core.gap,)+suffix)
                else:
                    for t in (0,1):
                        receipt = native_receipt(p,t)
                        check(receipt.source > receipt.child > 1 and receipt.endpoint > 1)
                for gap in (-7,-6,7,8):
                    candidate = compile_ladder(u,gap)
                    check(candidate is None or candidate in rows)
    counts = Counter(p.gap for p in packets)

    # Independent fixed-carrier comparison against all small positive compositions.
    table = {}
    for L in range(1,8):
        for C in range(L,14):
            for v in old.compositions(C,L):
                table[(L,C,old.carrier(v)[2])] = v
    comparison_count = 0
    for rel in (PLUS_ONE,TWO_TWOS):
        for ell in range(1,3):
            for u in product(range(1,7),repeat=ell):
                P,Q,B=old.carrier(u)
                for gap in range(-4,5):
                    C=sum(u)-2*gap;L=ell+rel.debt
                    if not L <= C <= 13: continue
                    scale=F(4)**gap
                    target=F(B+P*rel.offset)+(1<<C)*(scale-1)/3
                    check(target.denominator == 1)
                    answer=compile_ladder(u,gap,rel)
                    check((answer.partner if answer else None) == table.get((L,C,int(target))))
                    comparison_count += 1

    source_ladder=compile_ladder((12,1,1),-1)
    check(source_ladder.partner == (3,1,1,3,2,6))
    check(coarse_cell(source_ladder) == (75093,131072,16))
    bank=old.finite_head_bank()
    check(len(bank)==225)
    check(all(not (p.head[:3]==(12,1,1) or (12,1,1)[:len(p.head)]==p.head) for p in bank))
    check(F(4,131072//2)==F(1,16384))
    p3=phase(source_ladder,3)
    check((p3.residue,p3.period)==(485217,1048576))
    # A missing reserve is an actual positive source, not an impossible formal value.
    u=(12,1,1);res,mod=old.native_cell(u+(1,))
    X=res+mod
    Z=routes.replay(X,u)[-1]
    check(base.v2(3*Z+1)==1 and X % 131072 != 75093)
    check(compile_ladder((12,1,1),1) is None)

    j2=[compile_ladder((10,),1,TWO_TWOS),
        compile_ladder((1,12),1,TWO_TWOS),
        compile_ladder((1,1,14),1,TWO_TWOS)]
    check(j2[0].partner==(1,2,1,1,3))
    check(j2[1].partner==(2,1,1,2,4,1))
    check(j2[2].partner==(3,1,1,1,4,3,1))
    check([phase(p,2).residue % 16 for p in j2]==[9,5,13])
    check((phase(j2[2],2).residue,phase(j2[2],2).period)==(177949,262144))
    positive=compile_ladder((10,),1)
    coarse_residue,coarse_modulus,_=coarse_cell(positive)
    check(341 % coarse_modulus==coarse_residue)
    check(routes.replay(341,(10,))[-1]==1)
    rejects(lambda: routes.replay(341,(10,2)))
    grounded=native_receipt(positive,0)
    child_root=(3,1,1,3,3,2,1,3,1,1,3,4,1,3,1,2,3,4)
    check((grounded.source,grounded.child,grounded.endpoint)==(63829,2365,281))
    check(routes.replay(63829,discharge(grounded,child_root))[-1]==1)
    rejects(lambda: discharge(grounded,()))
    rejects(lambda: discharge(grounded,(4,)))
    controls=[(source_ladder,J) for J in range(3,9)]+[(positive,3)]+[(p,2) for p in j2]
    residual_receipts=0
    for packet,J in controls:
        record=phase(packet,J)
        audit_phase(record)
        for t in range(8):
            K=record.residue+record.period*(record.first_parameter+t)
            check(native_guard(packet,J,K))
            check(not native_guard(packet,J,K+2))
        _,_,A=coarse_cell(packet)
        K=2*A+J+1
        residue,modulus=residual_parameter(packet,K,J)
        for t in range(8):
            n=(1<<K)*(residue+modulus*t)-1
            for D in (3,4):
                receipt=apply_residual(packet,n,K,J,D)
                check(receipt.child==((n+1)>>D)-1)
                check(routes.replay(n,receipt.source_word)[-1]==receipt.endpoint)
                check(routes.replay(receipt.child,receipt.child_word)[-1]==receipt.endpoint)
                residual_receipts+=1

    # Equal endpoint heads carry a signed primitive predecessor; gap0 is not new by itself.
    zero=compile_ladder((10,1,4),0)
    core,a,suffix=reduce_zero(zero)
    check(core==positive and a==1 and suffix==(4,))
    check(mark(PLUS_ONE) != mark(TWO_TWOS))
    rejects(lambda: audit(replace(source_ladder,gap=True)))
    rejects(lambda: audit(replace(source_ladder,relation=Relation(3.0,-26))))
    rejects(lambda: audit(replace(source_ladder,partner=(3,1,1,3,2,6.0))))
    rejects(lambda: compile_ladder((1,True),1))
    rejects(lambda: mark(Relation(3,-25)))
    rejects(lambda: phase(positive,2))
    rejects(lambda: phase(j2[0],3))
    rejects(lambda: audit_phase(replace(p3,period=p3.period*2)))
    rejects(lambda: native_receipt(source_ladder,True))
    rejects(lambda: reduce_zero(source_ladder))
    check(routes.replay(native_receipt(source_ladder,0).source,
                        native_receipt(source_ladder,0).source_word)[-1]>1)

    print("PROVED: complete signed-gap decoder at a marked affine state; finite bounds per supplied head.")
    print("Fixed nonpadding head universe: letters1..12,lengths1..3; gap counts",sorted(counts.items()))
    print("Zero-gap matches reduce to signed ladder plus terminal/common suffix; no claimed new coverage from their count.")
    print("Independent carrier comparisons",comparison_count,"; actual residual receipts",residual_receipts)
    print("Inherited source ladder(12,1,1) gap-1; new exact coarse reserve X=75093 mod131072.")
    print("Entire coarse guard misses the unchanged225-bank; added relative >=3 head mass1/16384.")
    print("J3 Mersenne phaseK=485217 mod1048576; symbolic only, no huge ROOT discovery.")
    print("Changed J2 stateX=81Y+10: heads10,1-12,1-1-14 give exponent residues9,5,13 mod16.")
    print("Requested J2 branch positive control: head1114,phase177949 mod262144; memory owns grounding.")
    print("All receipts retain exact source, signed final valuation reserve and first-ROOT boundaries.")
    print("Formal coarse guard hostile341 reachesROOT under head10; padded continuation is rejected by strict exports.")
    print("One inherited frozen child ROOT word discharges63829; wrong-source and empty suffixes rejected.")
    print("No universal complement coverage or unsupplied child ROOT proof.")
    print("Exact checks",checks)


if __name__=='__main__':main()
