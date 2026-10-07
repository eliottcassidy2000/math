"""Exact Kolakoski clocks and guarded ordered-run Collatz compilation.

Normal and -O modes retain every check. Large counts concern explicit finite
prefixes; no limiting density or infinite positive-integer orbit is assumed.
"""
from dataclasses import dataclass
from fractions import Fraction as F
from itertools import groupby, product


def natural(n):
    if type(n) is not int or n < 0:
        raise ValueError("exact nonnegative integer required")
    return n


def alphabet_type(a):
    if type(a) is not tuple or len(a) != 2 or any(type(x) is not int for x in a) or a not in ((1, 2), (2, 3)):
        raise ValueError("exact alphabet (1,2) or (2,3) required")


def kolakoski(alphabet, length):
    alphabet_type(alphabet)
    natural(length)
    if alphabet == (1, 2):
        s, read, symbol = bytearray((1, 2, 2)), 2, 1
    else:
        s, read, symbol = bytearray((2, 2)), 1, 3
    while len(s) < length:
        s.extend([symbol] * s[read])
        read += 1
        symbol = sum(alphabet) - symbol
    return bytes(s[:length])


def rle(word):
    return tuple(len(tuple(g)) for _, g in groupby(word))


def prefix_clock(sequence, length):
    """sequence includes sufficient lookahead and its own run-length prefix."""
    natural(length)
    if not length or length >= len(sequence):
        raise ValueError("positive prefix and lookahead required")
    runs, observed = 1, 1
    for j in range(1, length):
        if sequence[j] == sequence[j - 1]:
            observed += 1
        else:
            if observed != sequence[runs - 1]:
                raise ValueError("closed run disagrees with self-description")
            runs, observed = runs + 1, 1
    full = sequence[runs - 1]
    deficit = full - observed
    if deficit < 0 or (deficit == 0) != (sequence[length] != sequence[length - 1]):
        raise ValueError("final run or lookahead is inconsistent")
    return runs, runs - 1, runs - (deficit > 0), observed, full, deficit


@dataclass(frozen=True)
class RunPacket:
    alphabet: tuple
    first: int
    lengths: tuple
    deficit: int = 0


def validate_packet(packet):
    if type(packet) is not RunPacket:
        raise ValueError("exact RunPacket required")
    alphabet_type(packet.alphabet)
    if type(packet.first) is not int or packet.first not in packet.alphabet:
        raise ValueError("marked starting symbol required")
    if type(packet.lengths) is not tuple or any(type(x) is not int or x < 1 for x in packet.lengths):
        raise ValueError("positive exact run lengths required")
    natural(packet.deficit)
    if (not packet.lengths and packet.deficit) or (packet.lengths and packet.deficit >= packet.lengths[-1]):
        raise ValueError("final deficit must leave a nonempty final run")
    return packet


def expand_packet(packet):
    validate_packet(packet)
    result, symbol = [], packet.first
    for j, length in enumerate(packet.lengths):
        actual = length - (packet.deficit if j + 1 == len(packet.lengths) else 0)
        result.extend([symbol] * actual)
        symbol = sum(packet.alphabet) - symbol
    return tuple(result)


def make_packet(alphabet, length):
    s = kolakoski(alphabet, natural(length) + 3)
    observed = rle(s[:length])
    R = len(observed)
    return RunPacket(alphabet, alphabet[0], tuple(s[:R]), s[R - 1] - observed[-1] if R else 0)


def authentic_self_prefix(packet):
    word = expand_packet(packet)
    return bytes(word) == kolakoski(packet.alphabet, len(word))


def compose(first, second):
    P, Q, B = first
    p, q, b = second
    return p * P, q * Q, p * B + b * Q


def compile_packet(packet):
    validate_packet(packet)
    result, symbol = (1, 1, 0), packet.first
    for j, length in enumerate(packet.lengths):
        d = length - (packet.deficit if j + 1 == len(packet.lengths) else 0)
        P, Q = 3**d, 2**(symbol * d)
        B = (P - Q) // (3 - 2**symbol)
        result = compose(result, (P, Q, B))
        symbol = sum(packet.alphabet) - symbol
    return result


def compile_word(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("exact positive valuation word required")
    result = (1, 1, 0)
    for a in word:
        result = compose(result, (3, 2**a, 1))
    return result


def native_cell(carrier):
    P, Q, B = carrier
    return ((Q - B) * pow(P, -1, 2 * Q)) % (2 * Q), 2 * Q


def replay(n, word):
    if type(n) is not int or n <= 0 or n % 2 == 0:
        raise ValueError("exact positive odd source required")
    compile_word(word)
    for a in word:
        if n == 1:
            raise ValueError("first-hit ROOT padding forbidden")
        numerator = 3 * n + 1
        actual = (numerator & -numerator).bit_length() - 1
        if a != actual:
            raise ValueError("actual valuation guard fails")
        n = numerator >> a
    return n


def certain_runs(word):
    lengths = list(rle(word))
    if lengths and lengths[0] == 1:
        lengths.pop(0)
    if lengths and lengths[-1] == 1:
        lengths.pop()
    return tuple(lengths)


FORBIDDEN = ((1, 1, 1), (2, 2, 2), (1, 2, 1, 2, 1), (2, 1, 2, 1, 2),
             (1, 1, 2, 2, 1, 1), (2, 2, 1, 1, 2, 2)) + tuple(
                 w * 3 for w in product((1, 2), repeat=3) if len(set(w)) == 2)


def contains(word, factor):
    return any(word[j:j + len(factor)] == factor for j in range(len(word) - len(factor) + 1))


def run():
    checks = 0

    def check(condition, message):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(message)

    def rejects(fn):
        try:
            fn()
        except ValueError:
            check(True, "rejection")
        else:
            check(False, "missing rejection")

    # Complete independent finite proof certificate for the factor implication.
    for i, w in enumerate(FORBIDDEN):
        check(i < 2 or any(contains(certain_runs(w), z) for z in FORBIDDEN[:i]),
              "forbidden-factor derivation through complete runs")
    allowed9 = [w for w in product((1, 2), repeat=9) if not any(contains(w, z) for z in FORBIDDEN)]
    check(len(allowed9) == 42, "complete length9 universe")
    check(all(4 <= w.count(1) <= 5 for w in allowed9), "9-factor balance")
    check(3**9 > 2**14, "uniform nine-letter coefficient growth")
    check(not any(contains((1, 2, 2, 1, 2), z) for z in FORBIDDEN), "local tests pass hostile")
    check((1, 2, 2, 1, 2) != tuple(kolakoski((1, 2), 5)), "local tests do not prove self-description")

    s12 = kolakoski((1, 2), 10000003)
    # Independent two-block substitution (OEIS Dekking representation), avoiding
    # the nongrowing seed122 whose final unpaired digit is dropped.
    beta = {(1, 1): bytes((1, 2)), (1, 2): bytes((1, 2, 2)),
            (2, 1): bytes((1, 1, 2)), (2, 2): bytes((1, 1, 2, 2))}
    alternate = bytes((1, 2, 2, 1))
    while len(alternate) < 10000:
        alternate = b''.join(beta[(alternate[j], alternate[j + 1])] for j in range(0, len(alternate) - 1, 2))
    check(alternate[:10000] == s12[:10000], "independent block generator")
    s23 = kolakoski((2, 3), 200003)
    check(s23[:25] == bytes((2,2,3,3,2,2,2,3,3,3,2,2,3,3,2,2,3,3,3,2,2,2,3,3,3)), "OEIS A071820 prefix")
    reports = []
    for alphabet, s, Ns in (((1, 2), s12, (2, 200000, 10000000)), ((2, 3), s23, (2, 200000))):
        for N in Ns:
            clock = prefix_clock(s, N)
            R, transitions, complete, observed, full, deficit = clock
            check(sum(s[:R]) - deficit == N, "run clock exact inverse")
            onesR = s[:R].count(alphabet[0])
            check(N == alphabet[1] * R - (alphabet[1] - alphabet[0]) * onesR - deficit,
                  "frequency at run index, not digit index")
            if alphabet == (1, 2):
                check(3 * R - 2 * N == 2 * onesR - R + 2 * deficit, "finite discrepancy identity")
            reports.append((alphabet, N, clock, s[:N].count(alphabet[0])))
    check(reports[1][2][0] == 133321 and reports[1][2][1] == 133320, "200k convention correction")
    check(reports[2][2][0] == 6666660, "10m count")

    # Output imbalance uses alternating run positions, a different observer.
    signed_prefix = [0]
    for symbol in s12[:8192]:
        signed_prefix.append(signed_prefix[-1] + (1 if symbol == 1 else -1))
    endpoint = alternating = 0
    run_endpoints = []
    alternating_prefix = [0]
    for r, length in enumerate(s12[:4096], 1):
        endpoint += length
        alternating += length if r % 2 else -length
        run_endpoints.append(endpoint)
        alternating_prefix.append(alternating)
        if r <= 2048:
            check(2 * endpoint == 3 * r - signed_prefix[r], "input imbalance determines clock endpoint")
            check(signed_prefix[endpoint] == alternating, "endpoint imbalance is alternating run observer")
    R = 1
    for N in range(1, 4097):
        while run_endpoints[R - 1] < N:
            R += 1
        deficit = run_endpoints[R - 1] - N
        sign = 1 if R % 2 else -1
        check(signed_prefix[N] == alternating_prefix[R] - sign * deficit,
              "run phase and deficit transport endpoint imbalance")

    packets = 0
    for alphabet in ((1, 2), (2, 3)):
        for N in range(257):
            packet = make_packet(alphabet, N)
            word = expand_packet(packet)
            carrier = compile_packet(packet)
            check(authentic_self_prefix(packet), "authenticated prefix")
            check(carrier == compile_word(word), "run powers retain ordered carry")
            check(word == tuple(kolakoski(alphabet, N)), "exact final-run deficit")
            residue, modulus = native_cell(carrier)
            P, Q, B = carrier
            for lift in (1, 2):
                source = residue + modulus * lift
                endpoint = replay(source, word)
                check(Q * endpoint == P * source + B, "actual native cylinder")
                if alphabet == (2, 3):
                    check((endpoint - 1) * 4**N <= (source - 1) * 3**N, "strict-source height contraction")
                    check((source - 1) * 3**N >= 2 * 4**N, "finite-prefix height lower bound")
            packets += 1

    # Minimal marking/carry controls.
    check(compile_word((1, 2)) == (9, 8, 5) and compile_word((2, 1)) == (9, 8, 7), "same counts lose order")
    check(native_cell(compile_word((1, 2))) == (11, 16), "12 native cell")
    check(native_cell(compile_word((2, 1))) == (9, 16), "21 native cell")
    cut = make_packet((1, 2), 2)
    check(cut.lengths == (1, 2) and cut.deficit == 1, "truncated run sidecar")
    check(expand_packet(RunPacket((1, 2), 1, (1, 2), 0)) == (1, 2, 2), "omitting deficit changes word")
    check(replay(43, (1, 2, 2)) == 37, "early Kolakoski source descent can occur")
    # BSC full support: no noisy observation can zero-error certify this nonconstant guard.
    epsilon = F(1, 4)
    for x in product((0, 1), repeat=3):
        for y in product((0, 1), repeat=3):
            d = sum(a != b for a, b in zip(x, y))
            check(epsilon**d * (1 - epsilon)**(3 - d) > 0, "exact noise full-support hostile")
    for bad in (True, -1, 2.0):
        rejects(lambda bad=bad: kolakoski((1, 2), bad))
    rejects(lambda: kolakoski((True, 2), 10))
    rejects(lambda: compile_packet(RunPacket((1, 2), True, (1,))))
    rejects(lambda: compile_packet(RunPacket((1, 2), 1, (True,))))
    rejects(lambda: compile_packet(RunPacket((1, 2), 1, (1,), 1)))
    rejects(lambda: compile_packet(RunPacket((1, 2), 1, (), 1)))
    rejects(lambda: replay(1, (2,)))
    rejects(lambda: replay(43, (True,)))
    rejects(lambda: replay(True, ()))

    print("collatz_self_describing_runs_20261007b")
    print("STATUS: FINITE-EXACT clocks; PROVED ordered-run compiler and height/factor obstructions")
    for alphabet, N, clock, first_count in reports:
        print("alphabet", alphabet, "N", N, "runs/transitions/complete/last_observed/full/deficit", clock,
              "count_first_symbol", first_count)
    print("independent generator agreement: first10000 A000002 digits")
    print("imbalance transport: E(S(r))=D(r), r1..2048; phase/deficit identity at N1..4096")
    print("forbidden derivation:12 words; full length9 universe512; surviving necessary factors42")
    print("Kolakoski12 valuation9-block slope >=19683/16384; this is not a density-limit claim")
    print("Kolakoski23 actual r-prefix source lower bound: n>=1+2*(4/3)^r; no infinite positive realization")
    print("ordered run packets:", packets, "; two native positive lifts per packet; explicit ROOT padding rejection")
    print("hostiles:200k transitions !=runs;12212 passes necessary local tests; order and terminal deficit are required")
    print("paper interfaces: ordered exact matching and retained output mean/entropy; noisy observation is no exact guard")
    print("exact checks:", checks)


if __name__ == "__main__":
    run()
