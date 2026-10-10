"""Authenticated Mersenne joins, their two clock defects, and ROOT transport.

Production consumes actual finite parity words. Only main() discovers the
explicitly declared small control universe. Giant fan data are inherited,
not independently re-computed by this checker.
"""
from dataclasses import dataclass, replace
from collections import Counter, defaultdict

CHECKS = 0
EXPONENT_CAP = 4096


def need(ok, message):
    global CHECKS
    CHECKS += 1
    if not ok:
        raise ValueError(message)


def integer(n, minimum=0):
    need(type(n) is int and n >= minimum, 'exact integer domain')
    return n


def exponent(e):
    integer(e, 2)
    need(e <= EXPONENT_CAP, 'declared literal exponent cap')
    return e


def parities(word):
    need(type(word) is tuple, 'exact tuple of parity bits')
    need(all(type(b) is int and b in (0, 1) for b in word), 'exact binary letters')
    return word


def step(n):
    integer(n, 1)
    return (3*n+1)//2 if n & 1 else n//2


def replay(n, word, root=False):
    integer(n, 1); parities(word)
    need(type(root) is bool, 'exact ROOT flag')
    states = [n]
    for b in word:
        need(n != 1, 'no ROOT padding')
        need(n & 1 == b, 'actual source parity')
        n = (3*n+1)//2 if b else n//2
        states.append(n)
    if root:
        need(n == 1, 'supplied word reaches ROOT')
    return tuple(states)


def odd_word_bits(word):
    need(type(word) is tuple and all(type(a) is int and a >= 1 for a in word),
         'positive exact valuation word')
    return tuple(b for a in word for b in (1,)+(0,)*(a-1))


@dataclass(frozen=True)
class Receipt:
    source_exponent: int
    child_exponent: int
    source_word: tuple
    child_word: tuple


def audit(receipt):
    need(type(receipt) is Receipt, 'exact Receipt')
    e = exponent(receipt.source_exponent)
    k = exponent(receipt.child_exponent)
    need(e > k, 'strictly smaller Mersenne child')
    x = replay((1 << e)-1, receipt.source_word)[-1]
    y = replay((1 << k)-1, receipt.child_word)[-1]
    need(x == y and x >= 3, 'same actual pre-ROOT endpoint')
    return x


def defect(receipt):
    audit(receipt)
    return (len(receipt.source_word)-len(receipt.child_word)
            -receipt.source_exponent+receipt.child_exponent,
            sum(receipt.source_word)-sum(receipt.child_word))


def compose(left, right):
    """Align the actual shared middle orbit, retaining the entry barriers."""
    audit(left); audit(right)
    need(left.child_exponent == right.source_exponent, 'same middle source')
    b, c = len(left.child_word), len(right.source_word)
    common = min(b, c)
    need(left.child_word[:common] == right.source_word[:common],
         'same authenticated middle prefix')
    result = Receipt(left.source_exponent, right.child_exponent,
                     left.source_word+right.source_word[b:],
                     right.child_word+left.child_word[c:])
    audit(result)
    dl, dr, dc = defect(left), defect(right), defect(result)
    need(dc == (dl[0]+dr[0], dl[1]+dr[1]), 'both defects add under alignment')
    return result


def discharge(receipt, supplied_child_root):
    """Use the supplied first-hit child proof; do not search for one."""
    audit(receipt)
    replay((1 << receipt.child_exponent)-1, supplied_child_root, root=True)
    b = len(receipt.child_word)
    need(supplied_child_root[:b] == receipt.child_word, 'actual child suffix seam')
    word = receipt.source_word+supplied_child_root[b:]
    replay((1 << receipt.source_exponent)-1, word, root=True)
    return word


def rooted_fibre(e, supplied_root):
    exponent(e)
    replay((1 << e)-1, supplied_root, root=True)
    return len(supplied_root)-e+1, sum(supplied_root)


def from_supplied_roots(e, k, left_root, right_root):
    """Extract a receipt from proofs already supplied, not new ROOT evidence."""
    exponent(e); exponent(k)
    xs = replay((1 << e)-1, left_root, root=True)
    ys = replay((1 << k)-1, right_root, root=True)
    index = {y:j for j, y in enumerate(ys) if y >= 3}
    candidates = [(i, index[x]) for i, x in enumerate(xs) if x in index]
    need(bool(candidates), 'first-hit positive Mersennes share at least 8')
    i, j = candidates[0]
    result = Receipt(e, k, left_root[:i], right_root[:j])
    audit(result)
    return result


def required_time_defect(e, bank_max_residual):
    """Necessary net time defect, conditional on a route to that ROOT bank.

    No source expansion: 3^12 > 2^19 implies s(e) > 19(e-1)/12.
    The bank maximum is an explicit supplied premise, not inferred here.
    """
    integer(e, 2); integer(bank_max_residual)
    return 19*(e-1)//12+1-bank_max_residual


def fan_members():
    """Exact label accounting for inherited THM-4605(4,6), not its replay."""
    top = 99708993713
    residual = top-8
    old_d = {3, 4, *range(7, 29)}
    bottom = residual-28
    escape_d = {*range(1, 1911), 1929, 1930, 1931, 1932, 1935, 1936}
    members = set(range(residual, top+1))
    members.update(residual-d for d in old_d)
    members.update(bottom-d for d in escape_d)
    return frozenset(members)


def reject(function, *args):
    try:
        function(*args)
    except (TypeError, ValueError):
        need(True, 'malformed receipt rejected')
    else:
        need(False, 'malformed receipt accepted')


def discover_control(e, cap=100000):
    """Explicit bounded discovery ONLY for the finite main() controls."""
    exponent(e); integer(cap, 1)
    n = (1 << e)-1
    bits = []
    for _ in range(cap):
        if n == 1:
            return tuple(bits)
        bits.append(n & 1)
        n = step(n)
    raise ValueError('finite control discovery cap; no nonconvergence inference')


def main():
    # Strict source authentication: 7=M_3, 3=M_2, endpoint 5.
    crossing = Receipt(3, 2, odd_word_bits((1, 1, 2, 3)), odd_word_bits((1,)))
    need(audit(crossing) == 5 and defect(crossing) == (5, 3),
         'smaller Mersenne join can cross the normalized deletion fibres')
    root3 = odd_word_bits((1, 4))
    root7 = discharge(crossing, root3)
    need(root7 == odd_word_bits((1, 1, 2, 3, 4)), 'literal ROOT suffix transfer')

    # This finite oracle is declared; production only consumes its words.
    proofs = {e:discover_control(e) for e in range(2, 97)}
    fibres = defaultdict(list)
    invariants = {}
    for e, word in proofs.items():
        invariants[e] = rooted_fibre(e, word)
        fibres[invariants[e]].append(e)
    pair_counts = Counter()
    for e in range(3, 97):
        for k in range(2, e):
            r = from_supplied_roots(e, k, proofs[e], proofs[k])
            d = defect(r)
            expected = tuple(invariants[e][i]-invariants[k][i] for i in range(2))
            need(d == expected, 'authenticated receipt has exact fibre difference')
            need(discharge(r, proofs[k]) == proofs[e], 'transport is canonical first-hit proof')
            pair_counts['zero' if d == (0, 0) else 'crossing'] += 1
            # THM-4605 converse: the shared 8 endpoint has equal normalized
            # time and accumulated odd count exactly in a rooted fibre.
            at8_e, at8_k = len(proofs[e])-3, len(proofs[k])-3
            zero8 = (at8_e-(e-1) == at8_k-(k-1)
                     and sum(proofs[e][:at8_e]) == sum(proofs[k][:at8_k]))
            need(zero8 == (d == (0, 0)), 'rooted converse retains both clocks')
    composition_count = 0
    for e in range(4, 33):
        for k in range(3, e):
            j = 2
            left = from_supplied_roots(e, k, proofs[e], proofs[k])
            right = from_supplied_roots(k, j, proofs[k], proofs[j])
            joined = compose(left, right)
            need(discharge(joined, proofs[j]) == proofs[e], 'composed receipt discharges exactly')
            composition_count += 1

    # Collapsing a zero-defect fibre reduces duplicated obligations, not
    # the number of different fibres requiring an authenticated terminal.
    need(len(fibres) == 22 and max(map(len, fibres.values())) == 14,
         'declared finite fibre census')
    same = from_supplied_roots(4, 3, proofs[4], proofs[3])
    need(defect(same) == (0, 0), '15=M_4 and 7=M_3 share a fibre')
    need(defect(compose(same, crossing)) == (5, 3), 'zero edge does not pay crossing debt')
    for e in range(2, 97):
        need(invariants[e][0] >= required_time_defect(e, 0), 'exact rational lower bound')
    need(3**12 > 2**19, 'rational logarithmic bound proved by integers')
    members = fan_members()
    need(len(members) == 1949 and min(members) == 99708991741,
         'inherited fan labels and deepest member')
    need(99708993677-1911 not in members and 99708993677-1936 in members,
         'absorbed set is not a contiguous interval')
    lower = required_time_defect(min(members), 97982)
    need(lower == 157872472274, 'necessary time charge to the inherited small bank')

    reject(audit, replace(crossing, source_exponent=3.0))
    reject(audit, replace(crossing, child_exponent=True))
    reject(audit, replace(crossing, source_word=(True,)+crossing.source_word[1:]))
    reject(audit, replace(crossing, source_word=crossing.source_word+(0,)))
    reject(audit, Receipt(3, 2, proofs[3], proofs[2])) # ROOT is not a pre-ROOT join
    reject(replay, 1, (1, 0), True)
    reject(discharge, crossing, proofs[3])
    reject(compose, crossing, same)
    reject(discover_control, 7, 1)
    reject(from_supplied_roots, 97, 2, proofs[96], proofs[2])
    reject(audit, replace(crossing, source_exponent=EXPONENT_CAP+1))
    reject(required_time_defect, True, 97982)

    print('PROVED: authenticated pre-ROOT Mersenne receipt defects add under exact clock alignment.')
    print('Hostile: M_3=7 ->5 <-3=M_2 has defect (5,3); normalized M_4 ->M_3 has (0,0).')
    print('FINITE-EXACT controls: exponents2..96; 95 supplied ROOT words; 22 fibres; largest14.')
    print('Pairs:', dict(sorted(pair_counts.items())), '; composed receipts:', composition_count)
    print('Inherited fan accounting only:1949 labels, bottom M_99708991741; giant ROOT remains OPEN.')
    print('Necessary net time defect to bank K<=12800 with max residual97982:', lower)
    print('No giant orbit replay; finite literal exponent cap:', EXPONENT_CAP)
    print('Exact checks:', CHECKS)


if __name__ == '__main__':
    main()
