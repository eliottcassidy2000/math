"""Marked four/six-point carry windows and an unexpanded D8 receipt template.

Tournament signs and motif counts are hostile quotients. Exact weighted
windows, terminal scale, ordered ports and an OPEN child are retained.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import combinations, product
from collections import Counter
from pathlib import Path
import sys

import collatz_eight_bit_completion_20261007c as eight
import tournament_recursive_four_state_20261004 as four


def require(ok, why):
    if not ok:
        raise ValueError(why)


def natural(n, least=0):
    require(type(n) is int and n >= least, "exact integer in declared domain")


@dataclass(frozen=True)
class Tail:
    offset: int


def terminal(value):
    if type(value) is Tail:
        natural(value.offset)
    else:
        natural(value, 1)


@dataclass(frozen=True)
class Window:
    points: tuple
    last: object


def encode_window(letters):
    require(type(letters) is tuple and len(letters) in (3, 5), "three/five-letter window")
    for a in letters[:-1]:
        natural(a, 1)
    terminal(letters[-1])
    p, q, b = 1, 1, 0
    points = [F(0)]
    for i, a in enumerate(letters):
        p, b = 3*p, 3*b+q
        points.append(F(b, p))
        if i+1 < len(letters):
            q <<= a
    return Window(tuple(points), letters[-1])


def decode_window(window):
    require(type(window) is Window and type(window.points) is tuple, "exact weighted window")
    z = window.points
    require(len(z) in (4, 6) and all(type(x) is F for x in z), "four/six rational vertices")
    terminal(window.last)
    require(z[0] == 0 and z[1] == F(1, 3), "marked translation and first gap")
    gaps = tuple(b-a for a, b in zip(z, z[1:]))
    require(all(d > 0 for d in gaps), "strict chronological order, without ties")
    letters = []
    for a, b in zip(gaps, gaps[1:]):
        power = 3*b/a
        require(power.denominator == 1 and power >= 2, "positive dyadic gap ratio")
        n = power.numerator
        require(n & (n-1) == 0, "gap ratio is an exact power of two")
        letters.append(n.bit_length()-1)
    answer = tuple(letters)+(window.last,)
    require(encode_window(answer) == window, "complete marked reconstruction")
    return answer


# Sparse exact polynomials in T=(2/3)^e and C=2^c.
# Keys are two nonnegative exponents; coefficients are Fractions.
def poly(terms):
    out = {}
    for (a, b), c in terms.items():
        natural(a); natural(b)
        require(type(c) is F, "exact rational coefficient")
        if c:
            out[a, b] = c
    return out


ZERO = {}
ONE = {(0, 0): F(1)}
T = {(1, 0): F(1)}
C = {(0, 1): F(1)}


def add(a, b):
    out = dict(a)
    for key, c in b.items():
        out[key] = out.get(key, F(0))+c
    return poly(out)


def scale(a, c):
    require(type(c) is F, "exact scale")
    return poly({key: c*v for key, v in a.items()})


def mul(a, b):
    out = {}
    for (i, j), x in a.items():
        for (k, l), y in b.items():
            key = i+k, j+l
            out[key] = out.get(key, F(0))+x*y
    return poly(out)


def glue(left, right):
    """(z,s)*(z',s')=(z+s*z',s*s'); order is chronological."""
    z, s = left
    zz, ss = right
    return add(z, mul(s, zz)), mul(s, ss)


def window_summary(window):
    word = decode_window(window)
    r = len(word)
    last = word[-1]
    a = sum(word[:-1])+(last.offset if type(last) is Tail else last)
    return ({(0, 0): window.points[-1]},
            {(0, int(type(last) is Tail)): F(1 << a, 3**r)})


def run_summary(drop=0):
    """The e-drop one-run, expressed using the same source-owned T."""
    natural(drop)
    s = scale(T, F(3, 2)**drop)
    return add(ONE, scale(s, F(-1))), s


def evaluate(p, t, c):
    require(type(t) is F and type(c) is F, "exact formal-variable specialization")
    return sum((v*t**i*c**j for (i, j), v in p.items()), F(0))


@dataclass(frozen=True)
class Packet:
    run: int
    oddpart: int
    left: tuple
    right: tuple
    terminal_status: str = "OPEN"


def make_packet(run, oddpart):
    natural(run, eight.MIN_RUN)
    natural(oddpart, 1)
    require(oddpart % 2 == 1, "source address uses its actual odd cofactor")
    a = eight.LEFT+(Tail(0),)
    b = eight.RIGHT+(Tail(2),)
    left = (encode_window(a[:3]), encode_window(a[3:]))
    right = tuple(encode_window(part) for part in (b[:5], b[5:10], b[10:13], b[13:]))
    return audit_packet(Packet(run, oddpart, left, right))


def audit_packet(packet):
    require(type(packet) is Packet, "exact compressed receipt packet")
    natural(packet.run, eight.MIN_RUN)
    natural(packet.oddpart, 1)
    require(packet.oddpart % 2 == 1, "canonical odd cofactor")
    require(type(packet.left) is tuple and type(packet.right) is tuple, "ordered window tuples")
    left = tuple(a for w in packet.left for a in decode_window(w))
    right = tuple(a for w in packet.right for a in decode_window(w))
    require(left == eight.LEFT+(Tail(0),) and right == eight.RIGHT+(Tail(2),),
            "exact ordered D8 heads and the linked final valuation")
    require(type(packet.terminal_status) is str and packet.terminal_status == "OPEN",
            "a word template cannot supply a child ROOT proof")
    residue, q = eight.audit_rule()
    require((2*pow(3, packet.run, 2*q)*packet.oddpart-1) % (2*q) == residue,
            "same source address, exact pulled-back head guard")
    z_left, s_left = run_summary()
    for window in packet.left:
        z_left, s_left = glue((z_left, s_left), window_summary(window))
    z_right, s_right = run_summary(eight.DROP)
    for window in packet.right:
        z_right, s_right = glue((z_right, s_right), window_summary(window))
    require(s_right == scale(s_left, F(1, 256)), "retained terminal clock ratio")
    require(z_right == scale(add(z_left, {(0, 0): F(255)}), F(1, 256)),
            "full carry relation, not only total valuation")
    return packet


def materialize(packet, bit_cap=4096):
    """Bind c at the same modest source; literal export still has a child obligation."""
    audit_packet(packet)
    natural(bit_cap, 1)
    bits = packet.run+1 if packet.oddpart == 1 else packet.run+1+packet.oddpart.bit_length()
    require(bits <= bit_cap, "explicit source-materialization cap")
    source = (packet.oddpart << (packet.run+1))-1
    receipt = eight.receipt(source, bit_cap)
    c = receipt.source_word[-1]
    def bind(windows):
        return tuple(c+a.offset if type(a) is Tail else a
                     for w in windows for a in decode_window(w))
    require(receipt.source_word == (1,)*packet.run+bind(packet.left), "ordered source export")
    require(receipt.child_word == (1,)*(packet.run-8)+bind(packet.right), "ordered child export")
    return receipt


def main():
    checks = Counter()
    def check(ok, key):
        checks[key] += 1
        if not ok: raise ArithmeticError(f"{key}: {checks[key]}")
    def rejects(fn):
        try: fn()
        except (ValueError, TypeError): check(True, "hostile")
        else: check(False, "hostile")

    # Intrinsic sign patterns remain transitive, while weights decode every word.
    for r in (3, 5):
        for word in product(range(1, 5), repeat=r):
            window = encode_window(word)
            check(decode_window(window) == word, "window_decode")
            z, s = window_summary(window)
            p, q, b = eight.routes.carrier(word)
            check(z == {(0, 0): F(b, p)} and s == {(0, 0): F(q, p)}, "summary")
            pts = window.points
            check(all(pts[j]-pts[i] > 0 for i, j in combinations(range(r+1), 2)), "transitive")
    words = [(1, 2, 3), (3, 1, 2), (1, 1, 2, 2, 3), (3, 2, 1, 2, 1)]
    for a, b, c in product(words, repeat=3):
        sa, sb, sc = [window_summary(encode_window(w)) for w in (a, b, c)]
        check(glue(glue(sa, sb), sc) == glue(sa, glue(sb, sc)), "composition")
        p, q, bb = eight.routes.carrier(a+b+c)
        check(glue(glue(sa, sb), sc) == ({(0, 0): F(bb, p)}, {(0, 0): F(q, p)}), "composition")
    a, b = words[:2]
    check(glue(window_summary(encode_window(a)), window_summary(encode_window(b))) !=
          glue(window_summary(encode_window(b)), window_summary(encode_window(a))), "order_loss")
    check(encode_window((1, 1, 1)).points == encode_window((1, 1, 4)).points,
          "terminal_loss")
    check(window_summary(encode_window((1, 1, 1))) != window_summary(encode_window((1, 1, 4))),
          "terminal_loss")
    check(four.class4(four.fixed_path(4, 3)) == four.class4(four.fixed_path(4, 4)) == 'S',
          "four_flip")
    check(four.class4(four.fixed_path(4, 3 ^ 1)) != four.class4(four.fixed_path(4, 4 ^ 1)),
          "four_flip")

    # Reproduce the inherited six-vertex collision; add the full four-motif census.
    sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
    import a000568_edge_perspective_extension_codex_s213 as six
    motif = []
    for mask in (344, 345):
        t = tuple(tuple(0 if i == j else int(six.tour.edge(mask, 6, i, j))
                        for j in range(6)) for i in range(6))
        motif.append(Counter(four.class4(tuple(tuple(t[i][j] for j in v) for i in v))
                             for v in combinations(range(6), 4)))
        check(six.tour.canonical(mask, 6) == mask, "six_collision")
        check(six.hamiltonian_paths(mask, 6) == 43, "six_collision")
        check(six.tour.score_sequence(mask, 6) == (2, 2, 2, 3, 3, 3), "six_collision")
    check(motif[0] == motif[1] == {'T': 1, '+': 2, '-': 2, 'S': 10}, "six_collision")
    for mode in ('size', 'internal'):
        check(six.class_sector_deck(344, 6, mode) == six.class_sector_deck(345, 6, mode), "six_collision")
    check(six.class_sector_deck(344, 6, 'cross') != six.class_sector_deck(345, 6, 'cross'), "six_collision")

    literals = 0
    for run in (68, 69, 96, 128):
        for t in range(4):
            source = eight.fixed_run_source(run, t)
            oddpart = (source+1) >> (run+1)
            packet = make_packet(run, oddpart)
            receipt = materialize(packet)
            c = receipt.source_word[-1]
            for windows, drop, word in ((packet.left, 0, receipt.source_word),
                                        (packet.right, 8, receipt.child_word)):
                summary = run_summary(drop)
                for w in windows: summary = glue(summary, window_summary(w))
                actual = tuple(evaluate(p, F(2, 3)**run, F(1 << c)) for p in summary)
                p, q, b = eight.routes.carrier(word)
                check(actual == (F(b, p), F(q, p)), "literal")
            check(receipt.child == ((source+1) >> 8)-1 and packet.terminal_status == "OPEN", "literal")
            literals += 1

    astronomical = make_packet(eight.PHASE-1, 1)
    check(astronomical.run+1-8 == 924745897, "symbolic")
    check(tuple(len(w.points) for w in astronomical.left) == (4, 6), "symbolic")
    check(tuple(len(w.points) for w in astronomical.right) == (6, 6, 4, 4), "symbolic")
    rejects(lambda: materialize(astronomical, 4096))
    rejects(lambda: audit_packet(replace(astronomical, terminal_status="ROOT")))
    rejects(lambda: audit_packet(replace(astronomical, run=astronomical.run+2)))
    rejects(lambda: audit_packet(replace(astronomical, left=tuple(reversed(astronomical.left)))))
    rejects(lambda: encode_window((1, True, 2)))
    rejects(lambda: decode_window(Window((F(0), F(1,3), F(1,3), F(2,3)), 1)))
    rejects(lambda: decode_window(Window((0, F(1,3), F(1), F(2)), 1)))
    rejects(lambda: encode_window((1, 2, Tail(-1))))
    rejects(lambda: make_packet(True, 1))
    print("PROVED: weighted four/six-point windows glue by (z,s)*(z',s')=(z+s*z',s*s').")
    print("Window universe: 1088 words, lengths3/5 and letters1..4; sign patterns all transitive.")
    print("Six-vertex masks344/345: same full4-motif census1,2,2,10 and H43; cross-sector arcs differ.")
    print("D8 fixed windows: source4+6 vertices; child6+6+4+4; formal T=(2/3)^e,C=2^c retained.")
    print("Exact identity: z_child=(z_source+255)/256, s_child=s_source/256; ordered words and guards retained.")
    print("Literal receipt controls", literals, "; symbolic child exponent924745897 remains OPEN and unexpanded.")
    print("Checks", dict(sorted(checks.items())))
    print("Total checks", sum(checks.values()))


if __name__ == '__main__':
    main()
