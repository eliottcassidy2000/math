"""Native cells, source-specific equality, and debt-observer losses.

The production compiler inspects supplied finite parity words, never discovers
an orbit or assumes that a smaller child is rooted.  T is the half-step map.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product

CHECKS = 0


def check(condition, message):
    global CHECKS
    CHECKS += 1
    if not condition:
        raise ValueError(message)


def integer(n, minimum=None):
    if type(n) is not int or (minimum is not None and n < minimum):
        raise ValueError('exact integer in the declared domain required')
    return n


def bits(word):
    if type(word) is not tuple or any(type(a) is not int or a not in (0, 1) for a in word):
        raise ValueError('tuple of exact parity bits required')
    return word


def valuations(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError('tuple of positive exact valuations required')
    return word


def T(n):
    integer(n, 1)
    return n // 2 if n % 2 == 0 else (3*n+1)//2


def carrier(word):
    bits(word)
    p, q, b = 1, 1, 0
    for a in word:
        m = 1+2*a
        p, q, b = m*p, 2*q, m*b+a*q
    return p, q, b


def native(word):
    p, q, b = carrier(word)
    return (-b*pow(p, -1, q)) % q, q


def literal(n, length):
    integer(n, 1); integer(length, 0)
    path, word = [n], []
    for _ in range(length):
        word.append(n % 2)
        n = T(n)
        path.append(n)
    return tuple(word), tuple(path)


@dataclass(frozen=True)
class Cell:
    offset: int
    source_bits: tuple
    child_bits: tuple
    kind: str
    residue: int
    modulus: int
    least: int
    slope_numerator: int
    carry_numerator: int
    point: object


def compile_cell(offset, source_bits, child_bits):
    """Complete positive-source collision set for one equal-length transcript.

    y is the child coordinate and x=y+offset.  T may formally continue at 1;
    the separate export adapter enforces strict first-ROOT boundaries.
    """
    integer(offset); bits(source_bits); bits(child_bits)
    if len(source_bits) != len(child_bits):
        raise ValueError('equal half-step word lengths required')
    p, q, b = carrier(source_bits)
    pp, _, bb = carrier(child_bits)
    rx, _ = native(source_bits)
    ry, _ = native(child_bits)
    cut = max(1, 1-offset)
    least = ry+q*(-((ry-cut)//q))
    d, c = p-pp, p*offset+b-bb
    kind, point = 'empty', None
    if (ry+offset-rx) % q == 0:
        if d == 0:
            if c == 0:
                kind = 'whole'
        elif (-c) % d == 0:
            z = (-c)//d
            if z >= cut and z % q == ry:
                kind, point = 'singleton', z
    return Cell(offset, source_bits, child_bits, kind, ry, q, least, d, c, point)


def audit_cell(cell):
    if type(cell) is not Cell:
        raise ValueError('typed Cell required')
    for n in (cell.offset, cell.residue, cell.modulus, cell.least,
              cell.slope_numerator, cell.carry_numerator):
        integer(n)
    if cell.point is not None:
        integer(cell.point, 1)
    if type(cell.kind) is not str or cell != compile_cell(
            cell.offset, cell.source_bits, cell.child_bits):
        raise ValueError('cell differs from its authenticated compiler output')
    return cell


def accepts(cell, child):
    audit_cell(cell); integer(child, 1)
    if child+cell.offset < 1 or child % cell.modulus != cell.residue:
        return False
    return cell.kind == 'whole' or (cell.kind == 'singleton' and child == cell.point)


def terminal_frame(cell):
    """Affine relation x_K=M*y_K+e, without evaluating y_K."""
    audit_cell(cell)
    p, q, b = carrier(cell.source_bits)
    pp, _, bb = carrier(cell.child_bits)
    m = F(p, pp)
    return m, (p*cell.offset+b-m*bb)/q


def replay_u(source, word, require_root=False):
    integer(source, 1); valuations(word)
    if source % 2 == 0 or type(require_root) is not bool:
        raise ValueError('positive odd source and boolean ROOT flag required')
    n = source
    for a in word:
        if n == 1:
            raise ValueError('earlier ROOT is not a strict receipt')
        z = 3*n+1
        actual = (z & -z).bit_length()-1
        if a != actual:
            raise ValueError('supplied valuation is not actual')
        n = z >> a
    if require_root and n != 1:
        raise ValueError('receipt does not end at ROOT')
    return n


@dataclass(frozen=True)
class Join:
    source: int
    child: int
    source_bits: tuple
    child_bits: tuple
    source_word: tuple
    child_word: tuple
    endpoint: int
    first_equal_time: int
    odd_meeting_time: int


def export_join(offset, child, source_bits, child_bits):
    """Actual strictly smaller odd-child join; no ROOT suffix is assumed.

    odd_meeting_time is the chosen odd endpoint after the supplied horizon,
    not generally the earliest common odd state. first_equal_time is earliest.
    """
    integer(offset, 1); integer(child, 1)
    source = child+offset
    if source % 2 == 0 or child % 2 == 0:
        raise ValueError('export requires both starting sources odd')
    cell = compile_cell(offset, source_bits, child_bits)
    if not accepts(cell, child):
        raise ValueError('supplied source does not satisfy the collision cell')
    wx, px = literal(source, len(source_bits))
    wy, py = literal(child, len(child_bits))
    if wx != source_bits or wy != child_bits or px[-1] != py[-1]:
        raise ValueError('literal transcript mismatch')
    if 1 in px[:-1] or 1 in py[:-1]:
        raise ValueError('formal ROOT padding cannot be exported')
    first = next(i for i, (x, y) in enumerate(zip(px, py)) if x == y)
    px, py = list(px), list(py)
    while px[-1] % 2 == 0:
        px.append(px[-1]//2)
        py.append(py[-1]//2)
    def odd_word(path):
        indices = [i for i, n in enumerate(path) if n % 2]
        return tuple(b-a for a, b in zip(indices, indices[1:]))
    uw, vw = odd_word(px), odd_word(py)
    if replay_u(source, uw) != px[-1] or replay_u(child, vw) != px[-1]:
        raise ValueError('odd-clock conversion mismatch')
    return Join(source, child, source_bits, child_bits, uw, vw, px[-1], first, len(px)-1)


def audit_join(join):
    if type(join) is not Join:
        raise ValueError('typed Join required')
    for n in (join.source, join.child, join.endpoint, join.first_equal_time, join.odd_meeting_time):
        integer(n, 1)
    valuations(join.source_word); valuations(join.child_word)
    if join != export_join(join.source-join.child, join.child, join.source_bits, join.child_bits):
        raise ValueError('join does not match its source and literal words')
    return join


def discharge(join, supplied_child_root):
    audit_join(join)
    replay_u(join.child, supplied_child_root, True)
    k = len(join.child_word)
    if supplied_child_root[:k] != join.child_word:
        raise ValueError('supplied child ROOT certificate does not contain the common future')
    result = join.source_word+supplied_child_root[k:]
    replay_u(join.source, result, True)
    return result


def mod2(x):
    if type(x) not in (int, F):
        raise ValueError('exact rational hidden coordinate required')
    x = F(x)
    if x.denominator % 2 == 0:
        raise ValueError('hidden coordinate is not 2-adically integral')
    return x.numerator % 2


def frame_step(m, e, j):
    """Standard Collatz pair update for one supplied reference digit."""
    if type(j) is not int or j not in (0, 1):
        raise ValueError('exact driver digit required')
    if mod2(m) != 1 or m <= 0:
        raise ValueError('positive 2-adic unit slope required')
    i = (mod2(m)*j+mod2(e)) % 2
    m, e = F(m), F(e)
    mm = m*F(1+2*i, 1+2*j)
    return mm, ((1+2*i)*e+i-mm*j)/2, i


def expect_bad(call):
    try:
        call()
    except (ValueError, TypeError):
        check(True, 'rejected hostile')
    else:
        check(False, 'malformed or ungrounded input was accepted')


def main():
    global CHECKS
    CHECKS = 0
    packets = 0
    for k in range(7):
        words = list(product((0, 1), repeat=k))
        for u in words:
            for v in words:
                for e in (-3, -2, 0, 1, 2, 16):
                    cell = compile_cell(e, u, v)
                    packets += 1
                    # Independent actual T replay at three complete native lifts.
                    for y in (cell.least, cell.least+cell.modulus, cell.least+2*cell.modulus):
                        wu, pu = literal(y+e, k)
                        wv, pv = literal(y, k)
                        truth = wu == u and wv == v and pu[-1] == pv[-1]
                        check(accepts(cell, y) == truth, 'native collision versus literal path')
                    if cell.kind == 'singleton':
                        wu, pu = literal(cell.point+e, k)
                        wv, pv = literal(cell.point, k)
                        check((wu, wv) == (u, v) and pu[-1] == pv[-1], 'exact singleton outside sample lifts')
    print('Complete transcript universe: lengths 0..6, all ordered parity-word pairs, offsets -3,-2,0,1,2,16.')
    print('Compiled packets:', packets)

    ux, uy = (1, 1, 1, 0, 0), (1, 1, 1, 0, 1)
    singleton = compile_cell(16, ux, uy)
    check((singleton.kind, singleton.point, singleton.residue, singleton.modulus) == ('singleton', 7, 7, 32), '23/7 native singleton')
    check(terminal_frame(singleton) == (F(1, 3), F(40, 3)), 'nonidentity terminal frame')
    join = export_join(16, 7, ux, uy)
    check((join.source_word, join.child_word, join.endpoint, join.first_equal_time, join.odd_meeting_time) ==
          ((1, 1, 5), (1, 1, 2, 3), 5, 5, 7), 'odd-clock singleton export')
    check(discharge(join, (1, 1, 2, 3, 4)) == (1, 1, 5, 4), 'supplied child grounding')
    check(not accepts(singleton, 39), 'same native cells need the original source')
    for y in (7, 39, 71):
        # Independent homogeneous-column determinant, with the odd-clock cost 7.
        det = (27*(y+16)+19)*128-(81*y+73)*128
        check(det == 128*(-54*y+378), 'marked source incidence determinant')
        check((det == 0) == (y == 7), 'incidence needs the supplied source')
    print('Singleton 23/7: first equal T-time 5 at 20; first common odd time 7 at 5; frame (1/3,40/3).')
    print('Supplied ROOT(7) transports to ROOT(23)=(1,1,5,4). Native lift 55/39 does not meet at this horizon.')

    u, v = (1, 1, 0, 0, 1, 0, 0), (1, 0, 1, 0, 0, 0, 1)
    whole = compile_cell(2, u, v)
    check((whole.kind, whole.residue, whole.modulus, whole.least) == ('whole', 49, 128, 49), 'whole native cylinder')
    check(terminal_frame(whole) == (1, 0), 'whole-cell structural identity')
    for t in range(64):
        j = export_join(2, 49+128*t, u, v)
        check(j.source == 51+128*t and j.first_equal_time <= 7, 'all-height positive control')
    first = export_join(2, 49, u, v)
    check((first.source_word, first.child_word, first.endpoint) == ((1, 3, 3), (2, 4, 1), 11), 'structural odd words')
    check(discharge(first, (2, 4, 1, 1, 2, 3, 4)) == (1, 3, 3, 1, 2, 3, 4), 'whole cell supplied ROOT control')
    lifted = export_join(2, 177, u, v)
    check((lifted.source_word, lifted.child_word, lifted.endpoint) == ((1, 3, 4), (2, 4, 2), 19), 'fixed odd words need the extra endpoint bit')
    extended = export_join(2, 49, u+(1,), v+(1,))
    check((extended.first_equal_time, extended.odd_meeting_time, extended.endpoint) == (7, 8, 17), 'chosen endpoint need not be the first common odd state')
    print('Whole cell: y=49+128t, x=y+2, common endpoint 11+27t after 7 half steps, t>=0.')

    # Exact finite partition rather than independent stream sampling.
    for k in range(1, 7):
        q = 1 << k
        for e in (2, 4, 16):
            exceptions = set()
            for r in range(1, q, 2):
                u = literal(r+e, k)[0]
                v = literal(r, k)[0]
                cell = compile_cell(e, u, v)
                if cell.kind == 'singleton':
                    exceptions.add(cell.point)
                for t in range(5):
                    y = r+t*q
                    px, py = literal(y+e, k)[1], literal(y, k)[1]
                    check(accepts(cell, y) == (px[-1] == py[-1]), 'partition leaf classification')
            check(len(exceptions) <= 2**(k-1), 'finite positive exceptional set bound')

    patterns = 0
    for k in range(1, 9):
        for word in product((0, 1), repeat=k):
            m, e, multiplier = F(1), F(2**k), 1
            for j in word:
                m, e, i = frame_step(m, e, j)
                multiplier *= 1+2*j
                check(i == j and m == 1, 'zero debt throughout finite agreement')
            check(e == multiplier and mod2(e) == 1, 'one missing source digit changes the next coupling')
            check(frame_step(m, e, 0)[2] == 1, 'next digits split')
            patterns += 1
    for k in range(1, 49):
        y, x = 2**(k+1)-1, 3*2**k-1
        wy, py = literal(y, k)
        wx, px = literal(x, k)
        check(wx == wy == (1,)*k, 'positive initial agreement run')
        check((py[-1], px[-1]) == (2*3**k-1, 3**(k+1)-1), 'sharp positive precision witness')
        check(py[-1] % 2 == 1 and px[-1] % 2 == 0 and 1 not in px+py, 'post-block split is pre-ROOT')
    print('Precision hostile: all', patterns, 'binary driver words through length 8; positive growing witnesses K=1..48.')
    print('Frames (1,0) and (1,2^K) have identical first K debt increments and zero K-block covariance, but different next couplings.')

    vectors = ((0, 0, 0), (0, 0, 0), (1, 0, 0), (0, 1, 0), (0, 0, 1))
    matrices, thirds = [], []
    for sign in (1, -1):
        increments = [tuple(vectors[(j+sign) % 5][a]-vectors[j][a] for a in range(3)) for j in range(5)]
        matrix = tuple(tuple(sum(z[a]*z[b] for z in increments) for b in range(3)) for a in range(3))
        matrices.append(matrix)
        thirds.append(sum(z[0]**2*z[1] for z in increments))
        check(all(sum(z[a] for z in increments) == 0 for a in range(3)), 'root-law martingale mean')
    check(matrices[0] == matrices[1] == ((2, -1, 0), (-1, 2, -1), (0, -1, 2)), 'orientation-blind covariance')
    check(thirds == [1, -1], 'orientation retained by a third mixed moment')
    print('Z5 translation couplings +/-1: identical 3x3 path covariance; E[xi1^2 xi2]=+/-1/5.')

    expect_bad(lambda: compile_cell(True, (), ()))
    expect_bad(lambda: compile_cell(2, (True,), (1,)))
    expect_bad(lambda: compile_cell(2, (1.0,), (1,)))
    expect_bad(lambda: compile_cell(2, [1], (1,)))
    expect_bad(lambda: compile_cell(2, (1,), ()))
    expect_bad(lambda: accepts(replace(singleton, point=7.0), 7))
    expect_bad(lambda: accepts(replace(singleton, kind='whole'), 7))
    expect_bad(lambda: accepts(singleton, True))
    expect_bad(lambda: export_join(16, 39, ux, uy))
    expect_bad(lambda: export_join(1, 5, (0, 1, 1, 0, 0), (1, 0, 0, 0, 1)))
    expect_bad(lambda: export_join(4, 1, literal(5, 4)[0], literal(1, 4)[0]))
    expect_bad(lambda: discharge(join, (4,)))
    expect_bad(lambda: discharge(replace(join, source=23.0), (1, 1, 2, 3, 4)))
    expect_bad(lambda: replay_u(True, (), True))
    expect_bad(lambda: replay_u(1.0, (), True))
    expect_bad(lambda: replay_u(1, (2,), True))
    expect_bad(lambda: frame_step(F(1, 2), 0, 0))
    expect_bad(lambda: frame_step(1, 0, True))
    expect_bad(lambda: frame_step(0, 0, 0))
    expect_bad(lambda: frame_step(-1, 0, 0))
    expect_bad(lambda: frame_step(2, 0, 0))
    check(replay_u(1, (), True) == 1, 'exact empty ROOT receipt')
    print('Hostile controls: 21 type, source, grounding, precision-domain, and ROOT-padding rejections.')
    print('CHECKS', CHECKS)
    print('Scope: exact finite transcript compiler and observer losses; no universal ROOT or Haar-to-integer inference.')


if __name__ == '__main__':
    main()
