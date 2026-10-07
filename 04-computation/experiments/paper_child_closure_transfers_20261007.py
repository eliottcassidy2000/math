"""Exact 12-run normal forms and guarded receiver transport.

Reproduce: python -B -O -X utf8 <this file>
No production routine discovers a ROOT word; splice/cut consume supplied words.
"""
from dataclasses import dataclass
from math import gcd


def nat(x, name="integer"):
    if type(x) is not int or x < 0:
        raise ValueError(name + " must be an exact nonnegative integer")
    return x


def odd(x):
    if type(x) is not int or x <= 0 or x % 2 == 0:
        raise ValueError("source must be an exact positive odd integer")
    return x


def valuation(n, p):
    if type(n) is not int or n <= 0 or type(p) is not int or p < 2:
        raise ValueError("positive integer and integer base >=2 required")
    v = 0
    while n % p == 0:
        n //= p
        v += 1
    return v


def step(n):
    odd(n)
    a = valuation(3 * n + 1, 2)
    return (3 * n + 1) >> a, a


def word_type(word):
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("word must be a tuple of exact positive integers")
    return word


def replay(n, word):
    """Validate a supplied actual word, forbidding any step after first ROOT."""
    odd(n)
    word_type(word)
    states = [n]
    for a in word:
        if n == 1:
            raise ValueError("ROOT padding forbidden")
        n, actual = step(n)
        if actual != a:
            raise ValueError("incorrect actual valuation")
        states.append(n)
    return tuple(states)


@dataclass(frozen=True)
class NormalForm:
    source: int
    depth: int
    core: int
    parameter: int


def normal_form(n):
    """Total exact maximal (1,2)-prefix decoder on positive odd integers."""
    odd(n)
    r = (valuation(n + 5, 2) - 1) // 3
    t = (n + 5) // (2 * 8**r)
    return NormalForm(n, r, 2 * 9**r * t - 5, t)


def prefix_pair(r, t):
    nat(r, "depth")
    nat(t, "parameter")
    if t < (3 if r == 0 else 1):
        raise ValueError("parameter gives a nonpositive source")
    return 2 * 8**r * t - 5, 2 * 9**r * t - 5


def expand(parent, r):
    """Inverse 12-prefix, defined iff 9^r divides parent+5."""
    odd(parent)
    nat(r, "depth")
    q, rem = divmod(parent + 5, 9**r)
    if rem:
        raise ValueError("ternary divisibility guard fails")
    return odd(8**r * q - 5)


def inverse_depth(core):
    odd(core)
    if normal_form(core).depth:
        raise ValueError("a canonical core is required")
    return valuation(core + 5, 3) // 2


def pull_receiver_guard(r, a, modulus):
    """Return a t congruence or None; prefix_pair retains the positivity cut."""
    nat(r, "depth")
    if type(a) is not int or type(modulus) is not int or modulus < 1:
        raise ValueError("exact residue and positive modulus required")
    coefficient = 2 * 9**r
    g = gcd(coefficient, modulus)
    if (a + 5) % g:
        return None
    period = modulus // g
    residue = ((a + 5) // g * pow(coefficient // g, -1, period)) % period
    return residue, period


def word_carrier(word):
    word_type(word)
    P, Q, B = 1, 1, 0
    for a in word:
        P, Q, B = 3 * P, Q * 2**a, 3 * B + Q
    return P, Q, B


def receiver_carrier(r, P, Q, B):
    """Transport an authenticated affine receiver; does not prove its guards."""
    nat(r, "depth")
    if type(P) is not int or P < 1 or type(Q) is not int or Q < 1 or type(B) is not int:
        raise ValueError("exact positive P,Q and integer B required")
    return P * 9**r, Q * 8**r, 5 * P * (9**r - 8**r) + B * 8**r


def forward_depth_budget(word):
    """Largest r with receiver slope*(9/8)^r<1; -1 means no such r."""
    P, Q, _ = word_carrier(word)
    if P >= Q:
        return -1
    r, p, q = 0, P, Q
    while 9 * p < 8 * q:
        r, p, q = r + 1, 9 * p, 8 * q
    return r


def mersenne_child(K, r):
    if type(K) is not int or K < 5 or K % 2 == 0:
        raise ValueError("odd K>=5 required")
    if type(r) is not int or r < 1:
        raise ValueError("positive exact depth required")
    return expand(2 ** (K - 2) - 1, r)


def decode_mersenne_child(n):
    f = normal_form(n)
    if f.depth == 0:
        return None
    z = f.core + 1
    if z & (z - 1):
        return None
    exponent = z.bit_length() - 1
    if exponent < 3 or exponent % 2 == 0:
        return None
    return exponent + 2, f.depth


@dataclass(frozen=True)
class RootReceipt:
    source: int
    word: tuple


def validate_root(receipt):
    if type(receipt) is not RootReceipt or replay(receipt.source, receipt.word)[-1] != 1:
        raise ValueError("a strict supplied ROOT receipt is required")
    return receipt


def expand_root(receipt, r):
    validate_root(receipt)
    n = expand(receipt.source, r)
    result = RootReceipt(n, (1, 2) * r + receipt.word)
    return validate_root(result)


def retract_root(receipt):
    validate_root(receipt)
    f = normal_form(receipt.source)
    if receipt.word[:2 * f.depth] != (1, 2) * f.depth:
        raise ValueError("supplied ROOT word misses the mandatory prefix")
    return validate_root(RootReceipt(f.core, receipt.word[2 * f.depth:]))


@dataclass(frozen=True)
class ReceiverReceipt:
    source: int
    depth: int
    core: int
    receiver_word: tuple
    endpoint: int
    paid: bool
    reached_root: bool


def attach_receiver(n, receiver_word):
    f = normal_form(n)
    endpoint = replay(f.core, receiver_word)[-1]
    replay(n, (1, 2) * f.depth + receiver_word)
    return ReceiverReceipt(n, f.depth, f.core, receiver_word, endpoint,
                           endpoint < n, endpoint == 1)


def run():
    checks = 0

    def check(condition, label):
        nonlocal checks
        checks += 1
        if not condition:
            raise RuntimeError(label)

    def rejects(fn, label):
        try:
            fn()
        except ValueError:
            check(True, label)
        else:
            check(False, label)

    # Independent literal decoder on the complete positive odd interval.
    for n in range(1, 20000, 2):
        f = normal_form(n)
        x, literal = n, 0
        while x != 1:
            y, a = step(x)
            z, b = step(y)
            if (a, b) != (1, 2):
                break
            x, literal = z, literal + 1
        check((literal, x) == (f.depth, f.core), "literal normal form")
        check(f.parameter % 8 != 0, "maximal parameter")
        check(prefix_pair(f.depth, f.parameter) == (n, f.core), "decode pair")
        check(expand(f.core, f.depth) == n, "inverse reconstruction")
        check(normal_form(f.core).depth == 0, "core idempotence")
        check(f.depth <= inverse_depth(f.core), "finite inverse chain")

    pairs = 0
    for r in range(9):
        for t in range(3 if r == 0 else 1, 65):
            n, m = prefix_pair(r, t)
            states = replay(n, (1, 2) * r)
            check(states[-1] == m, "prefix pair")
            check(all(x > n for x in states[1:]), "every proper state grows")
            check(normal_form(n).depth == r + valuation(t, 2) // 3, "depth composition")
            pairs += 1
    guard_cases = 0
    for r in range(5):
        for modulus in (1, 2, 3, 4, 8, 9, 16, 19, 27, 64, 81, 243):
            for a in range(modulus):
                cell = pull_receiver_guard(r, a, modulus)
                for t in range(3, 35):
                    actual = (2 * 9**r * t - 5 - a) % modulus == 0
                    predicted = cell is not None and (t - cell[0]) % cell[1] == 0
                    check(actual == predicted, "nonunit guard pullback")
                    guard_cases += 1

    family = []
    for K in range(5, 606, 2):
        m = 2**(K - 2) - 1
        R = valuation(m + 5, 3) // 2
        check(R == (1 + valuation(K - 4, 3)) // 2, "inherited LTE depth")
        for r in range(1, R + 1):
            n = mersenne_child(K, r)
            check(decode_mersenne_child(n) == (K, r), "tagged family decoder")
            check(replay(n, (1, 2) * r)[-1] == m, "family actual return")
            check(valuation(n + 5, 2) == 3 * r + 2, "child depth sidecar")
            family.append((K, r))

    # Finite supplied-certificate controls: this test helper alone discovers words.
    roots = []
    for b in range(1, 512, 2):
        if normal_form(b).depth:
            continue
        x, w = b, []
        for _ in range(2000):
            if x == 1:
                break
            x, a = step(x)
            w.append(a)
        check(x == 1, "finite ROOT control terminates within stated cap")
        receipt = RootReceipt(b, tuple(w))
        validate_root(receipt)
        for r in range(inverse_depth(b) + 1):
            lifted = expand_root(receipt, r)
            check(retract_root(lifted) == receipt, "ROOT expansion retraction")
            check(len(lifted.word) == len(w) + 2 * r, "strict odd-rank increment")
            roots.append(lifted)

    # Exact normalized-odd cylinder counts, with a retained unresolved tail.
    R = 4
    counts = [0] * (R + 1)
    for n in range(1, 2 * 8**R, 2):
        counts[min(R, normal_form(n).depth)] += 1
    expected = [7 * 8**(R - r - 1) for r in range(R)] + [1]
    check(counts == expected, "depth partition")

    hostile = attach_receiver(3067, (2,))
    check((hostile.depth, hostile.core, hostile.endpoint) == (3, 4369, 3277), "source-relative hostile")
    check(hostile.endpoint < hostile.core and not hostile.paid, "core descent is unpaid")
    check(not attach_receiver(27, ()).paid, "return alone pays nothing")
    check(decode_mersenne_child(11) is None, "12 child need not be Mersenne child")
    check(pull_receiver_guard(2, 0, 9) is None, "incompatible ternary target")
    check(pull_receiver_guard(2, 4, 9) == (0, 1), "forced ternary target")
    for receipt in roots:
        f = normal_form(receipt.source)
        tail = receipt.word[2 * f.depth:]
        P, Q, B = word_carrier(tail)
        p, q, b = receiver_carrier(f.depth, P, Q, B)
        check(p * receipt.source + b == q, "transported ROOT carrier")
        check(((q - p) * receipt.source > b) == (receipt.source > 1), "original-source payment iff")
    for r in range(9):
        for t in range(3, 20):
            n, m = prefix_pair(r, t)
            for P, Q, B in ((3, 4, 1), (1, 64, -63), (729, 1024, 669)):
                p, q, b = receiver_carrier(r, P, Q, B)
                check((p * n + b) * Q == (P * m + B) * q, "signed affine interface")
    check(receiver_carrier(2, 3, 4, 1) == (243, 256, 319), "r2 forward receiver threshold")
    check(receiver_carrier(3, 3, 4, 1) == (2187, 2048, 3767), "r3 forward receiver cannot pay")
    for word in ((), (1,), (2,), (3,), (4,), (1, 2), (1, 4), (2, 4), (6, 8, 10)):
        budget = forward_depth_budget(word)
        P, Q, _ = word_carrier(word)
        check((budget == -1) == (P >= Q), "empty forward budget")
        if budget >= 0:
            check(P * 9**budget < Q * 8**budget, "last possible receiver depth")
        check(P * 9**(budget + 1) >= Q * 8**(budget + 1), "first forbidden receiver depth")
    check(forward_depth_budget((2,)) == 2, "word2 exact budget")
    for bits in range(1, 25):
        r = bits // 3 + 1
        n, _ = prefix_pair(r, 1)
        n2, _ = prefix_pair(r + 1, 1)
        check(n % 2**bits == n2 % 2**bits and normal_form(n).depth != normal_form(n2).depth,
              "no fixed dyadic observation decodes all depths")
    # Signed fixed-point hostile is external to every positive-domain API.
    check((3 * -5 + 1) // 2 == -7 and (3 * -7 + 1) // 4 == -5, "signed 12 cycle")
    for x in (True, 1.0, -5, 0, 2):
        rejects(lambda x=x: normal_form(x), "typed positive domain")
    for r in (True, 1.0, -1):
        rejects(lambda r=r: expand(13, r), "typed depth")
    rejects(lambda: replay(1, (2,)), "ROOT padding")
    rejects(lambda: replay(3, (True,)), "boolean letter")
    rejects(lambda: replay(3, [1]), "noncanonical word")
    rejects(lambda: expand(7, 1), "illegal inverse")
    rejects(lambda: validate_root(RootReceipt(27, (1, 2))), "unclosed prefix")
    rejects(lambda: attach_receiver(27, (2,)), "receiver actual guard")
    rejects(lambda: pull_receiver_guard(1, False, 3), "typed residue")
    rejects(lambda: pull_receiver_guard(1, 1, True), "typed modulus")
    rejects(lambda: mersenne_child(4, 1), "typed family domain")
    rejects(lambda: inverse_depth(11), "canonical core guard")

    print("paper_child_closure_transfers_20261007")
    print("STATUS: PROVED elementary normal form / CONDITIONAL ROOT splicing / FINITE-EXACT controls")
    print("complete literal universe: 10000 positive odd sources below 20000")
    print("prefix-pair controls:", pairs)
    print("receiver guard comparisons:", guard_cases)
    print("Mersenne-child controls: odd K=5..605, legal depths:", len(family))
    print("supplied ROOT transport controls, canonical cores below512:", len(roots))
    print("depth counts q=0,1,2,3,>=4 on4096 normalized odd residues:", counts)
    print("hostile: 3067 --(12)^3--> 4369 --2--> 3277; core descends, original source unpaid")
    print("receiver word2: r2 payment iff13*n>319; r3 slope2187/2048 forbids payment")
    print("boundary: signed -5 --1--> -7 --2--> -5 is excluded; no finite dyadic depth decoder")
    print("paper transfer: exact boundary / fixed lattice / retained conditioning / fixed-word scope")
    print("exact checks:", checks)


if __name__ == "__main__":
    run()
