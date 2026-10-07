"""Unique positive-word decoder for the specified two-anchor head template.

Exact affine matching is separate from a supplied-source guard. The native
export below constructs actual common-future receipts; it never finds ROOT.
"""
from dataclasses import dataclass, replace
from fractions import Fraction
from itertools import combinations, product
from collections import Counter

import collatz_uncovered_join_routes_20261007 as routes


def integer(n, lower=None):
    if type(n) is not int or (lower is not None and n < lower):
        raise ValueError("exact integer with stated lower bound required")
    return n


def word_type(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("tuple of exact positive valuation letters required")
    return word


def carrier(word):
    word_type(word)
    P, Q, B = 1, 1, 0
    for a in word:
        P, Q, B = 3*P, (1 << a)*Q, 3*B+Q
    return P, Q, B


def decode(length, cost, carry):
    """Return the unique positive word, or None for impossible exact data."""
    integer(length, 1)
    integer(cost, 0)
    integer(carry)
    if cost < length or carry < 1:
        return None
    initial = (3**length, 1 << cost, carry)
    power = initial[0]//3
    letters = []
    for remaining in range(length, 1, -1):
        difference = carry-power
        if difference <= 0:
            return None
        a = (difference & -difference).bit_length()-1
        if a < 1 or cost-a < remaining-1:
            return None
        letters.append(a)
        cost -= a
        carry = difference >> a
        power //= 3
    if carry != 1 or cost < 1:
        return None
    result = tuple(letters+[cost])
    if carrier(result) != initial:
        raise ArithmeticError("decoded word did not reconstruct exact data")
    return result


@dataclass(frozen=True)
class HeadPartner:
    head: tuple
    gap: int
    partner: tuple


def compile_partner(head, gap=1):
    word_type(head)
    integer(gap, 1)
    P, Q, B = carrier(head)
    length = len(head)+3
    cost = sum(head)-2*gap
    if cost < length:
        return None
    target_carry = B-26*P+(1 << cost)*((4**gap-1)//3)
    partner = decode(length, cost, target_carry)
    return None if partner is None else HeadPartner(head, gap, partner)


def all_partners(head):
    word_type(head)
    limit = (sum(head)-len(head)-3)//2
    return tuple(p for r in range(1, limit+1)
                 if (p := compile_partner(head, r)) is not None)


def audit_partner(packet):
    if type(packet) is not HeadPartner:
        raise ValueError("exact HeadPartner required")
    word_type(packet.head)
    word_type(packet.partner)
    integer(packet.gap, 1)
    if compile_partner(packet.head, packet.gap) != packet:
        raise ValueError("partner packet does not authenticate its full word")
    return packet


def append_terminal(packet, terminal=1):
    audit_partner(packet)
    integer(terminal, 1)
    return packet.head+(terminal,), packet.partner+(terminal+2*packet.gap,)


def native_cell(word):
    P, Q, B = carrier(word)
    return ((Q-B)*pow(P, -1, 2*Q)) % (2*Q), 2*Q


def native_receipt(packet, parameter):
    """All positive native children, with odd terminal letters and no ROOT hit."""
    audit_partner(packet)
    integer(parameter, 0)
    left, right = append_terminal(packet, 1)
    residue, modulus = native_cell(right)
    child = residue+modulus*parameter
    if child < 3:
        raise ArithmeticError("odd-terminal native child cannot be ROOT")
    source = 27*child-26
    P, Q, B = carrier(right)
    endpoint = (P*child+B)//Q
    return routes.audit(routes.Receipt(source, child, left, right, endpoint))


def discharge(packet, parameter, supplied_child_word):
    """Consume a supplied first-hit child ROOT word, never discover one."""
    return routes.discharge(native_receipt(packet, parameter), supplied_child_word)


def prefix_minimal(words):
    selected = []
    for w in sorted(set(words), key=lambda w: (len(w), w)):
        word_type(w)
        if not any(w[:len(v)] == v for v in selected):
            selected.append(w)
    return tuple(selected)


def antichain_report(packets):
    for packet in packets:
        audit_partner(packet)
    heads = prefix_minimal(p.head for p in packets if p.head and p.head[0] != 2)
    out = []
    for group in (1, 3):
        words = tuple(w for w in heads if (w[0] == 1 if group == 1 else w[0] >= 3))
        shift = 1 if group == 1 else 2
        mass = sum((Fraction(1, 1 << (sum(w)-shift)) for w in words), Fraction())
        out.append((group, len(words), tuple(sorted(Counter(map(sum, words)).items())), mass))
    return heads, tuple(out)


def finite_head_bank():
    """Declared finite universe only; one canonical gap per minimal head."""
    packets = tuple(p for length in range(1, 5)
                    for u in product(range(1, 13), repeat=length) if u[0] != 2
                    for p in all_partners(u))
    heads = prefix_minimal(p.head for p in packets)
    by_head = {}
    for packet in packets:
        previous = by_head.get(packet.head)
        if previous is None or packet.gap < previous.gap:
            by_head[packet.head] = packet
    return tuple(by_head[u] for u in heads)


def compositions(cost, length):
    if cost < length:
        return
    for cuts in combinations(range(1, cost), length-1):
        points = (0,)+cuts+(cost,)
        yield tuple(points[i+1]-points[i] for i in range(length))


def main():
    checks = 0

    def check(ok):
        nonlocal checks
        checks += 1
        if not ok:
            raise ArithmeticError("exact control failed")

    def rejects(fn):
        try:
            fn()
        except (ValueError, TypeError):
            check(True)
        else:
            check(False)

    # Independent exhaustive carrier table, built without using the decoder.
    table = {}
    composition_count = 0
    for length in range(1, 8):
        for cost in range(length, 17):
            for w in compositions(cost, length):
                P, Q, B = carrier(w)
                key = (length, cost, B)
                check(key not in table)
                table[key] = w
                check(decode(*key) == w)
                composition_count += 1
    small_data_count = 0
    for length in range(1, 8):
        for cost in range(17):
            for B in range(256):
                check(decode(length, cost, B) == table.get((length, cost, B)))
                small_data_count += 1

    packets = []
    head_count = 0
    for length in range(1, 5):
        for u in product(range(1, 13), repeat=length):
            head_count += 1
            P, Q, B = carrier(u)
            A = sum(u)
            answers = all_partners(u)
            packets.extend(answers)
            by_gap = {p.gap: p for p in answers}
            for gap in range(1, max(1, (A-length-3)//2+1)):
                C, L = A-2*gap, length+3
                target = B-26*P+(1 << C)*((4**gap-1)//3)
                if C <= 16:
                    expected = table.get((L, C, target))
                    actual = by_gap[gap].partner if gap in by_gap else None
                    check(actual == expected)
            for packet in answers:
                p, q, b = carrier(packet.partner)
                check(len(packet.partner) == length+3)
                check(sum(packet.partner) == A-2*packet.gap)
                check(Fraction(p+b, q) == 4**packet.gap*Fraction(P+B, Q)
                      + Fraction(4**packet.gap-1, 3))
                for terminal in (1, 2, 5):
                    left, right = append_terminal(packet, terminal)
                    lp, lq, lb = carrier(left)
                    rp, rq, rb = carrier(right)
                    check(lq == rq and rp == 27*lp and rb == lb-26*lp)
                for t in (0, 1):
                    receipt = native_receipt(packet, t)
                    check(receipt.child < receipt.source and receipt.endpoint > 1)
                    check(routes.replay(receipt.source, receipt.source_word)[-1] == receipt.endpoint)
                    check(routes.replay(receipt.child, receipt.child_word)[-1] == receipt.endpoint)

    counts = Counter((len(p.head), p.gap) for p in packets)
    check(head_count == 22620 and len(packets) == 274)
    check(counts == {(1,1):1, (2,1):7, (3,1):40, (3,3):1,
                     (4,1):221, (4,2):3, (4,3):1})
    primary = tuple(p for p in packets if p.gap == 1)
    r1_heads, r1_report = antichain_report(primary)
    all_heads, all_report = antichain_report(packets)
    bank = finite_head_bank()
    check(tuple(p.head for p in bank) == all_heads and len(bank) == 225)
    check(len({p.head for p in packets if p.head[0] != 2}) == len(bank))
    for words, report in ((r1_heads, r1_report), (all_heads, all_report)):
        check(all(not (a != b and b[:len(a)] == a) for a in words for b in words))
        check(all(row[-1] <= 1 for row in report))
    # Leading2 is redundant only for anchor evaluation, not for actual source guards.
    for p in packets:
        if p.head[0] == 2:
            smaller = compile_partner(p.head[1:], p.gap)
            check(smaller is not None and p.partner == (2,)+smaller.partner)

    short_checks = 0
    short_answers = []
    for L in (4, 5):
        for C in range(L, 25):
            a = C+2 if L == 4 else C+1
            u = (a,) if L == 4 else (1, a)
            P, Q, B = carrier(u)
            for v in compositions(C, L):
                p, q, b = carrier(v)
                if (p+b)*Q == (4*(P+B)+Q)*q:
                    short_answers.append((u, v))
                short_checks += 1
    check(short_checks == 53130)
    check(short_answers == [((10,), (3,1,1,3)), ((1,10), (1,1,2,2,3))])
    check([((n & -n).bit_length()-1) for n in (104,176,224,256)] == [3,4,5,8])
    check([((n & -n).bit_length()-1) for n in (310,364,400,448,512)] == [1,2,4,6,9])

    p10 = compile_partner((10,))
    receipt = native_receipt(p10, 0)
    check((receipt.source, receipt.child, receipt.endpoint) == (63829,2365,281))
    child_root = (3,1,1,3,3,2,1,3,1,1,3,4,1,3,1,2,3,4)
    result = discharge(p10, 0, child_root)
    check(len(result) == 15 and routes.replay(63829, result)[-1] == 1)
    check(Fraction(sum(carrier((10,))[::2]), carrier((10,))[1]) == Fraction(1,256))
    check(compile_partner((11,)) is None and compile_partner((1,11)) is None)
    check(compile_partner(()) is None and all_partners(()) == ())
    for args in ((2,2,3),(2,2,2),(2,2,4),(2,1,5),(1,1,3),(3,9,-1)):
        check(decode(*args) is None)
    for bad in (True, 4.0, -1):
        rejects(lambda bad=bad: decode(bad, 8, 175))
    rejects(lambda: decode(4, 8.0, 175))
    rejects(lambda: decode(4, 8, True))
    for bad in ((True,), (10.0,), [10], (0,)):
        rejects(lambda bad=bad: compile_partner(bad))
    rejects(lambda: compile_partner((10,), True))
    rejects(lambda: audit_partner(replace(p10, gap=1.0)))
    rejects(lambda: audit_partner(replace(p10, partner=(3,1,1,True))))
    rejects(lambda: audit_partner(replace(p10, partner=(3,1,2,2))))
    rejects(lambda: native_receipt(p10, True))
    rejects(lambda: native_receipt(p10, -1))
    rejects(lambda: routes.replay(1, receipt.source_word))
    rejects(lambda: discharge(p10, 0, ()))
    rejects(lambda: discharge(p10, 0, child_root+(2,)))
    rejects(lambda: discharge(p10, 1, child_root))

    print("PROVED: unique positive-word decoder and exact finite sibling-gap partner list.")
    print("Decoder independent composition universe", composition_count,
          "; complete small data triples", small_data_count)
    print("Head universe: letters1..12, lengths1..4; heads", head_count,
          "; partners", len(packets), "; r1", len(primary))
    print("Counts(length,gap)", sorted(counts.items()))
    print("r1 prefix-minimal nonpadding antichain", len(r1_heads), "; rows", r1_report)
    print("all-gap prefix-minimal nonpadding antichain", len(all_heads), "; rows", all_report)
    print("Conditional weights are native head-cylinder masses; no Mersenne/payment coverage inference here.")
    print("New r1 examples:", [(p.head,p.partner) for p in primary if len(p.head)==2 and p.head[0]!=2])
    print("New nonpadding higher gaps:", [(p.head,p.gap,p.partner) for p in packets if p.gap>1 and p.head[0]!=2])
    print("Short-head independent child compositions", short_checks, "; unique pairs", short_answers)
    print("All-height native exports", 2*len(packets), "checked instances; x=27y-26; terminal oddness excludes ROOT padding.")
    print("Supplied-child example63829 ->281 <-2365; authenticated child ROOT word yields15 source edges.")
    print("Hostiles: formal value at1 is not actual; missing/forged partner, source guard and ROOT suffix are rejected.")
    print("No universal head success, no new Mersenne guard claim, no ROOT discovery in production.")
    print("Exact checks", checks)


if __name__ == '__main__':
    main()
