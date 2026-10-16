"""Source-marked elliptic lifts of actual Collatz words; no discovery oracle.

Run normally and with -O.  The JSON report is deterministic; this program never
writes files.  All claimed all-height statements are proved in the companion.
"""
from dataclasses import dataclass, replace
from fractions import Fraction as F
from itertools import product
import json


def need(ok, message):
    if not ok:
        raise ValueError(message)


def integer(x, lower=0):
    need(type(x) is int and x >= lower, 'exact integer in range')
    return x


def rational(x):
    need(type(x) in (int, F), 'exact rational')
    return F(x)


def odd(x, lower=1):
    integer(x, lower)
    need(x % 2 == 1, 'odd integer')
    return x


def word(w):
    need(type(w) is tuple, 'valuation tuple')
    for a in w:
        integer(a, 1)
    return w


def v2(n):
    integer(n, 1)
    return (n & -n).bit_length()-1


@dataclass(frozen=True)
class Frame:
    P: int
    Q: int
    B: int


def carrier(w):
    word(w)
    P = Q = 1
    B = 0
    for a in w:
        P, Q, B = 3*P, (1 << a)*Q, 3*B+Q
    return Frame(P, Q, B)


def decode_frame(frame):
    """Total membership reader: recover a positive word, or reject."""
    need(type(frame) is Frame, 'frame type')
    P, Q, B = frame.P, frame.Q, frame.B
    integer(P, 1)
    integer(Q, 1)
    integer(B)
    need(Q & (Q-1) == 0, 'power-of-two denominator')
    A = Q.bit_length()-1
    m = 0
    p = P
    while p > 1 and p % 3 == 0:
        p //= 3
        m += 1
    need(p == 1, 'power-of-three numerator')
    if m == 0:
        need((P, Q, B) == (1, 1, 0), 'empty frame')
        return ()
    N = B % Q
    need(A >= m and 0 < N < Q and N % 2 == 1, 'real dyadic phase')
    out = []
    while m > 1:
        Z = (N-pow(3, m-1, 1 << A)) % (1 << A)
        need(Z > 0, 'nonzero tail residue')
        a = v2(Z)
        need(a >= 1 and A-a >= m-1, 'positive valuation budget')
        out.append(a)
        N, A, m = Z >> a, A-a, m-1
    need(N == 1 and A >= 1, 'last valuation')
    out = tuple(out)+(A,)
    need(carrier(out) == frame, 'full carry, not only its real phase')
    return out


def native_cell(frame):
    decode_frame(frame)
    modulus = 2*frame.Q
    return ((frame.Q-frame.B)*pow(frame.P, -1, modulus)) % modulus, modulus


def replay(source, w):
    """Strict ROOT semantics: ROOT may occur only at the final endpoint."""
    odd(source)
    word(w)
    n = source
    for a in w:
        need(n != 1, 'no ROOT padding')
        z = 3*n+1
        need(v2(z) == a, 'actual valuation required')
        n = z >> a
    return n


def model(source, frame):
    """Cubic coefficients (a,b,c,d) for Y²=aX³+bX²+cX+d."""
    integer(source, 2)
    decode_frame(frame)
    u, v = F(frame.P, frame.Q), F(frame.B, frame.Q)
    A = source*(1-source)
    return (u, -3*u*v, 3*u*v*v+A*u**3, -u*v**3-A*u**3*v)


def on_curve(point, coefficients):
    need(type(point) is tuple and len(point) == 2, 'point tuple')
    need(type(coefficients) is tuple and len(coefficients) == 4, 'cubic tuple')
    x, y = map(rational, point)
    a, b, c, d = map(rational, coefficients)
    return y*y == ((a*x+b)*x+c)*x+d


def double(point, coefficients):
    need(on_curve(point, coefficients), 'point lies on supplied curve')
    x, y = map(rational, point)
    a, b, c, _ = map(rational, coefficients)
    need(a != 0 and y != 0, 'finite nonvertical doubling')
    slope = (3*a*x*x+2*b*x+c)/(2*y)
    z = (slope*slope-b)/a-2*x
    return z, -y+slope*(x-z)


@dataclass(frozen=True)
class Lift:
    source: int
    frame: Frame
    endpoint: int
    clock: int
    point: tuple


def compile_lift(source, w):
    odd(source, 3)
    endpoint = replay(source, w)
    frame = carrier(w)
    u = F(frame.P, frame.Q)
    return Lift(source, frame, endpoint, sum(w), (F(endpoint), u*u*source))


def audit_lift(receipt):
    need(type(receipt) is Lift, 'lift type')
    odd(receipt.source, 3)
    odd(receipt.endpoint)
    integer(receipt.clock)
    need(type(receipt.point) is tuple and len(receipt.point) == 2, 'point type')
    for coordinate in receipt.point:
        rational(coordinate)
    w = decode_frame(receipt.frame)
    need(replay(receipt.source, w) == receipt.endpoint, 'source-bound endpoint')
    need(receipt.clock == sum(w), 'literal shortcut clock')
    u = F(receipt.frame.P, receipt.frame.Q)
    need(receipt.point == (F(receipt.endpoint), u*u*receipt.source), 'marked point')
    need(on_curve(receipt.point, model(receipt.source, receipt.frame)), 'curve equation')
    return w


def compose(receipt, continuation):
    w = audit_lift(receipt)
    replay(receipt.endpoint, continuation)
    r = carrier(continuation)
    f = receipt.frame
    composed = Frame(r.P*f.P, r.Q*f.Q, r.P*f.B+r.B*f.Q)
    out = compile_lift(receipt.source, w+continuation)
    need(out.frame == composed, 'chronological affine composition')
    c, d = F(r.P, r.Q), F(r.B, r.Q)
    need(out.point == (c*receipt.point[0]+d, c*c*receipt.point[1]), 'point transport')
    return out


def discharge(receipt, supplied_child, supplied_root_word):
    audit_lift(receipt)
    odd(supplied_child)
    need(supplied_child == receipt.endpoint, 'unchanged supplied child')
    need(replay(supplied_child, supplied_root_word) == 1, 'supplied first-hit ROOT word')
    return compose(receipt, supplied_root_word)


def terminal_column(receipt):
    audit_lift(receipt)
    return (-1 if receipt.clock % 2 else 1), receipt.endpoint


def paid(receipt):
    """Strict dependency on a smaller integer; does not ground that integer."""
    audit_lift(receipt)
    f = receipt.frame
    return (f.P-f.Q)*receipt.source+f.B < 0


def golden_mul(x, y):
    """Q(phi), phi²=phi+1; pairs encode a+b*phi."""
    a, b = x
    c, d = y
    return a*c+b*d, a*d+b*c+b*d


def golden_inv(x):
    a, b = x
    norm = a*a+a*b-b*b
    need(norm != 0, 'nonzero golden element')
    return F(a+b, norm), F(-b, norm)


def golden_sub(x, y):
    return x[0]-y[0], x[1]-y[1]


def golden_chart(z):
    # h(z)=(z-phi')/(z-phi), phi'=1-phi.
    return golden_mul(golden_sub(z, (1, -1)), golden_inv(golden_sub(z, (0, 1))))


def main():
    checks = 0
    def check(ok, message):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(message)
    def rejects(call, message):
        try:
            call()
        except ValueError:
            check(True, message)
        else:
            check(False, message)

    frames = actual = split_controls = 0
    for length in range(6):
        for w in product(range(1, 5), repeat=length):
            f = carrier(w)
            check(decode_frame(f) == w, 'full frame round trip')
            frames += 1
            r, modulus = native_cell(f)
            for k in range(3):
                s = r+k*modulus
                if s < 3:
                    continue
                # Independently replay exact valuations, including ROOT padding.
                n = s
                first_hit_safe = True
                for a in w:
                    first_hit_safe = first_hit_safe and n != 1
                    z = 3*n+1
                    check(z % (1 << (a+1)) == 1 << a, 'native congruence supplies exact valuation')
                    n = z >> a
                if not first_hit_safe:
                    rejects(lambda s=s, w=w: compile_lift(s, w), 'ROOT-padded native cylinder')
                    continue
                receipt = compile_lift(s, w)
                check(audit_lift(receipt) == w, 'strict lift audit')
                check(receipt.endpoint == n, 'independent native replay')
                check(paid(receipt) == (n < s), 'immutable-source payment')
                u, v = F(f.P, f.Q), F(f.B, f.Q)
                check((receipt.point[0]-v)/u == s and receipt.point[1]/u**2 == s, 'inverse chart')
                base = (F(s), F(s))
                doubled = double(base, model(s, Frame(1, 1, 0)))
                check(doubled[0] == F(4*s*s-4*s+1, 4), 'nonintegral double coordinate')
                check(doubled[0].denominator == 4, 'Nagell-Lutz obstruction premise')
                transported_double = (u*doubled[0]+v, u*u*doubled[1])
                check(double(receipt.point, model(s, f)) == transported_double, 'group law commutes with chart')
                for cut in range(length+1):
                    check(compose(compile_lift(s, w[:cut]), w[cut:]) == receipt, 'all source-preserving seams')
                    split_controls += 1
                actual += 1

    root_ray = []
    for k in range(2, 65):
        s = (4**k-1)//3
        lift = compile_lift(s, (2*k,))
        check(lift.endpoint == 1 and audit_lift(lift) == (2*k,), 'all-height ROOT ray control')
        check(lift.point[1] == F(9*s, 16**k), 'non-torsion point at encoded ROOT')
        if k <= 5:
            root_ray.append({'k': k, 'source': s, 'point': [str(x) for x in lift.point]})
    prefix = compile_lift(3, (1,))
    rooted = discharge(prefix, 5, (4,))
    check(rooted == compile_lift(3, (1, 4)), 'supplied strict ROOT splice')
    check(rooted.point == (F(1), F(243, 1024)), 'ROOT is not torsion')
    check(discharge(rooted, 1, ()) == rooted, 'empty terminal suffix')
    check(not paid(prefix) and paid(rooted), 'endpoint coordinates, not curve rank, pay')
    left, right = terminal_column(prefix), terminal_column(compile_lift(5, ()))
    check((left, right) == ((-1, 5), (1, 5)), 'integer-terminal signed cancellation')
    check(model(3, prefix.frame) != model(5, Frame(1, 1, 0)), 'coalescence is not curve-model equality')
    check(prefix.point != (F(5), F(5)), 'coalescence is not marked-point equality')
    check(native_cell(carrier((1, 2))) == (11, 16), 'first same-slope source cell')
    check(native_cell(carrier((2, 1))) == (9, 16), 'second same-slope source cell')
    check(carrier((1, 2)) == Frame(9, 8, 5) and carrier((2, 1)) == Frame(9, 8, 7), 'slope does not determine carry')

    # Exact conjugacy; no decimal approximation to phi or logarithms is used.
    golden_count = 0
    multiplier = (-2, 1)  # -phi^-2 = phi-2.
    for denominator in range(1, 8):
        for numerator in range(-30, 31):
            z = F(numerator, denominator)
            if z == 1:
                continue
            lhs = golden_chart(golden_inv((z-1, 0)))
            rhs = golden_mul(multiplier, golden_chart((z, 0)))
            check(lhs == rhs, 'golden Mobius conjugacy')
            golden_count += 1
    check(golden_mul(multiplier, golden_mul((0, 1), (0, 1))) == (-1, 0), 'golden multiplier exact inverse square')

    hostiles = [
        lambda: carrier((True,)), lambda: carrier([1]),
        lambda: decode_frame(Frame(True, 1, 0)),
        lambda: decode_frame(Frame(3, 2.0, 1)),
        lambda: decode_frame(Frame(9, 8, 13)),
        lambda: decode_frame(Frame(6, 8, 5)),
        lambda: decode_frame(Frame(1, 2, 0)),
        lambda: compile_lift(1, ()), lambda: compile_lift(2, ()),
        lambda: compile_lift(3.0, (1,)), lambda: compile_lift(3, (2,)),
        lambda: compile_lift(3, (1, 4, 2)),
        lambda: audit_lift(replace(rooted, source=3.0)),
        lambda: audit_lift(replace(rooted, endpoint=True)),
        lambda: audit_lift(replace(rooted, clock=5.0)),
        lambda: audit_lift(replace(rooted, point=(1.0, rooted.point[1]))),
        lambda: audit_lift(replace(rooted, point=(F(1), F(1)))),
        lambda: audit_lift(replace(rooted, frame=Frame(9, 32, 37))),
        lambda: discharge(prefix, 7, (1, 1, 2, 3, 4)),
        lambda: discharge(prefix, 5, (4, 2)),
        lambda: discharge(rooted, 1, (2,)),
        lambda: compose(rooted, (2,)),
    ]
    for call in hostiles:
        rejects(call, 'typed/native/ROOT hostile')
    report = {
        'checks': checks,
        'scope': 'Exact finite controls plus proved source-marked lift; no new universal or source coverage.',
        'word_universe': {'lengths': [0, 5], 'valuations': [1, 4], 'frames': frames,
                          'native_first_hit_lifts': actual, 'all_splits': split_controls},
        'root_ray_controls_k': [2, 64], 'root_ray_samples': root_ray,
        'supplied_root_3': {'word': [1, 4], 'frame': [9, 32, 5],
                            'point': ['1', '243/1024'], 'shortcut_clock': 5},
        'coalescence_3_5': {'terminal': 5, 'clocks': [1, 0], 'signs': [-1, 1]},
        'same_slope_different_carry': [[9, 8, 5], [9, 8, 7]],
        'golden_exact_conjugacy_controls': golden_count,
        'typed_and_semantic_hostiles': len(hostiles),
        'open_obligation': 'A source-preserving strict rank or uniform signed-energy budget beyond this isomorphic decoration.'
    }
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
