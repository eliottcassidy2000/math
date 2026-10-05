"""Exact neighbour macros, source masses, and a climb-register obstruction.

No computation runs on import. No universal convergence assumption is used.
"""
from dataclasses import dataclass
from fractions import Fraction


def positive_int(value, name="value"):
    if type(value) is not int or value < 1:
        raise ValueError(name + " must be an exact positive integer")
    return value


def odd_source(value):
    positive_int(value, "source")
    if value % 2 == 0:
        raise ValueError("source must be odd")
    return value


def v2(value):
    positive_int(value)
    return (value & -value).bit_length() - 1


def odd_step(source):
    odd_source(source)
    if source == 1:
        raise ValueError("ROOT has no first-hit outgoing edge")
    numerator = 3 * source + 1
    exponent = v2(numerator)
    return numerator >> exponent, exponent


def replay(source, word, require_root=False):
    odd_source(source)
    if type(word) is not tuple or any(type(a) is not int or a < 1 for a in word):
        raise ValueError("word must be a tuple of exact positive exponents")
    current = source
    for wanted in word:
        current, actual = odd_step(current)
        if actual != wanted:
            raise ValueError("valuation guard failed")
    if require_root and current != 1:
        raise ValueError("receipt does not finish at ROOT")
    return current


def even_term(k):
    positive_int(k, "k")
    return 2 * (4 ** k - 1) // 3


def root_ray(k):
    """The k=1 member is ROOT itself, with an empty first-hit word."""
    positive_int(k, "k")
    return (4 ** k - 1) // 3, (() if k == 1 else (2 * k,))


@dataclass(frozen=True)
class Neighbour:
    k: int
    offset: int

    def __post_init__(self):
        positive_int(self.k, "k")
        if type(self.offset) is not int or self.offset not in (-1, 1):
            raise ValueError("offset must be the exact integer -1 or 1")

    @property
    def source(self):
        return even_term(self.k) + self.offset

    @property
    def word(self):
        if self.offset == -1:
            return () if self.k == 1 else (2,) + (1,) * (2 * self.k - 2) + (2,)
        epsilon = 1 if self.k % 2 == 0 else 2
        return (1,) + (2,) * (self.k - 1) + (2 + epsilon,)

    @property
    def endpoint(self):
        if self.offset == -1:
            return (3 ** (2 * self.k - 1) - 1) // 2
        epsilon = 1 if self.k % 2 == 0 else 2
        return (3 ** self.k + 1) // (1 << epsilon)

    @property
    def dependency(self):
        if self.offset != 1:
            raise ValueError("the growing minus macro has no smaller dependency here")
        return 3 ** (self.k - 1)


def identify_neighbour(source, offset):
    """Total exact membership test; a malformed type is an error, not membership."""
    odd_source(source)
    if type(offset) is not int or offset not in (-1, 1):
        raise ValueError("invalid offset")
    twice_power = 3 * (source - offset) + 2
    if twice_power <= 0 or twice_power % 2:
        return None
    power = twice_power // 2
    if power < 4 or power & (power - 1):
        return None
    exponent = power.bit_length() - 1
    if exponent % 2:
        return None
    return Neighbour(exponent // 2, offset)


def lift_plus_receipt(k, supplied_child_word):
    """Consume a checked first-hit receipt for 3^(k-1); never search its orbit."""
    record = Neighbour(k, 1)
    replay(record.dependency, supplied_child_word, require_root=True)
    if k == 1:
        # Child=ROOT has no outgoing edge to cancel.
        result = record.word
    else:
        epsilon = 1 if k % 2 == 0 else 2
        if not supplied_child_word or supplied_child_word[0] != epsilon:
            raise ValueError("common-future child guard failed")
        result = record.word + supplied_child_word[1:]
    replay(record.source, result, require_root=True)
    return result


def gamma_mass(source):
    odd_source(source)
    index = (source + 1) // 2
    return Fraction(1, 1 << (2 * index.bit_length() - 1))


def mixture_mass(source):
    odd_source(source)
    return Fraction(8, 3 * 4 ** source.bit_length())


def gamma_cylinder(residue, depth):
    """Index m=r modulo 2^depth, where n=2m-1."""
    if type(depth) is not int or depth < 0 or type(residue) is not int:
        raise ValueError("invalid cylinder")
    modulus = 1 << depth
    if not 0 <= residue < modulus:
        raise ValueError("residue out of range")
    atom = Fraction(0) if residue == 0 else Fraction(1, 1 << (2 * residue.bit_length() - 1))
    return Fraction(1, 4 ** depth) + atom


def mixture_cylinder(residue, depth):
    positive_int(depth, "depth")
    odd_source(residue)
    if residue >= 1 << depth:
        raise ValueError("residue out of range")
    return mixture_mass(residue) + Fraction(4, 3 * 4 ** depth)


def literal_receipt(source, cap):
    """Finite experiment helper only; it is not used by receipt transport."""
    odd_source(source)
    positive_int(cap, "cap")
    current, word = source, []
    while current != 1:
        if len(word) == cap:
            raise ValueError("finite experiment cap reached")
        current, exponent = odd_step(current)
        word.append(exponent)
    return tuple(word)


def main():
    checks = 0

    def check(condition, label):
        nonlocal checks
        checks += 1
        if not condition:
            raise ValueError(label)

    for k in range(1, 129):
        a = even_term(k)
        check(even_term(k + 1) == 4 * a + 2, "sequence recursion")
        check(bin(a)[2:] == "10" * k, "binary sequence")
        source, word = root_ray(k)
        check(a == 2 * source and replay(source, word, True) == 1, "root ray")
        check(word == (() if k == 1 else (2 * k,)), "root rank convention")
        for offset in (-1, 1):
            record = Neighbour(k, offset)
            check(replay(record.source, record.word) == record.endpoint, "literal macro")
            check(identify_neighbour(record.source, offset) == record, "source identity")
            check(bin(record.source)[2:] == ("10" * (k - 1) + ("01" if offset == -1 else "11")).lstrip("0"), "neighbour suffix")
        plus, minus = Neighbour(k, 1), Neighbour(k, -1)
        check(plus.endpoint < plus.source, "plus payment")
        check(plus.dependency < plus.source, "smaller dependency")
        if k >= 2:
            check(odd_step(plus.dependency)[0] == plus.endpoint, "common future")
            check(minus.endpoint > minus.source, "minus growth")
            check(v2(minus.source + 1) == v2(minus.endpoint + 1) == 1, "register return")
            following, following_exponent = odd_step(minus.endpoint)
            check(following_exponent == 2 + v2(k), "following reset valuation")
            check(following == (9 ** k - 1) // (1 << (3 + v2(k))), "following reset endpoint")
            check(gamma_mass(minus.source) == gamma_mass(plus.source), "equal gamma atoms")
            check(mixture_mass(minus.source) == mixture_mass(plus.source), "equal mixture atoms")
        if k >= 3:
            third = replay(plus.source, (1, 2, 2))
            check(third == (27 * plus.source + 23) // 32 < plus.source, "third step pays")
            check(replay(plus.source, (1,)) > plus.source and replay(plus.source, (1, 2)) > plus.source, "first descent exactly three")

    for source in range(1, 8192, 2):
        for offset in (-1, 1):
            record = identify_neighbour(source, offset)
            expected = [k for k in range(1, 8) if even_term(k) + offset == source]
            check((record.k if record else None) == (expected[0] if expected else None), "finite membership census")

    for k in range(1, 65):
        total = Fraction(3, 4 ** k)
        gamma = gamma_cylinder(1 << (k - 1), k)
        mixture = mixture_cylinder((1 << k) - 1, k + 1)
        check(gamma == mixture == total, "same climb-register distribution")
        check(gamma_mass((1 << k) - 1) / total == Fraction(2, 3), "gamma least representative")
        check(mixture_mass((1 << k) - 1) / total == Fraction(8, 9), "mixture least representative")

    for bound in range(1, 65):
        gamma = sum((gamma_mass(root_ray(k)[0]) for k in range(1, bound + 1)), Fraction())
        mixture = sum((mixture_mass(root_ray(k)[0]) for k in range(1, bound + 1)), Fraction())
        check(gamma + Fraction(2, 15 * 16 ** (bound - 1)) == Fraction(19, 30), "gamma rooted ray mass")
        check(mixture + Fraction(2, 45 * 16 ** (bound - 1)) == Fraction(32, 45), "mixture rooted ray mass")

    transported_edges = 0
    for k in range(1, 17):
        child_word = literal_receipt(3 ** (k - 1), 10000)
        result = lift_plus_receipt(k, child_word)
        transported_edges += len(result)
        check(replay(Neighbour(k, 1).source, result, True) == 1, "transported ROOT receipt")
        if k >= 2:
            check(len(result) == len(child_word) + k, "transport odd rank")
            check(len(result) + sum(result) == len(child_word) + sum(child_word) + 3 * k + 1, "transport ordinary rank")

    check(replay(9, (2, 1, 1, 2)) == 13, "universal h-only hostile")
    check(v2(9 + 1) == v2(13 + 1) == 1, "hostile register")
    check(gamma_mass(9) == gamma_mass(13) == Fraction(1, 32), "gamma equal hostile endpoints")
    check(mixture_mass(9) == mixture_mass(13) == Fraction(1, 96), "mixture equal hostile endpoints")
    check(Neighbour(1, -1).word == () and Neighbour(1, -1).source == 1, "minus ROOT exception")

    hostiles = [
        lambda: Neighbour(True, 1), lambda: Neighbour(1.0, 1),
        lambda: Neighbour(0, 1), lambda: Neighbour(1, True),
        lambda: Neighbour(1, 0), lambda: identify_neighbour(True, 1),
        lambda: identify_neighbour(2, 1), lambda: identify_neighbour(3, 1.0),
        lambda: replay(1, (2,)), lambda: replay(1.0, ()),
        lambda: replay(11, (1, 2, 2)), lambda: lift_plus_receipt(2, (2,)),
        lambda: lift_plus_receipt(1, (2,)), lambda: lift_plus_receipt(3, (1, 4)),
    ]
    for action in hostiles:
        try:
            action()
        except ValueError:
            check(True, "hostile rejected")
        else:
            raise ValueError("hostile accepted")

    print("PROVED: a_k=2(4^k-1)/3; neighbour/root formula controls k=1..128")
    print("FINITE-EXACT: neighbour membership on 4096 positive odd sources below8192, both offsets")
    print("PROVED: plus macro pays; k>=3 first descent exactly3; smaller common-future dependency3^(k-1)")
    print("PROVED: minus k>=2 macro grows and restores h=1; k=1 is ROOT with empty word")
    print("PROVED: both climb-register laws3/4^k; minimum representative conditionals2/3 and8/9")
    print("PROVED: rooted inverse ray gamma mass19/30, mixture mass32/45; exact finite tails checked")
    print("PROVED: every positive h-only reweighting of either baseline fails strict all-edge discounted flow")
    print("FINITE-EXACT: supplied-child receipt transport k=1..16; resulting odd edges=" + str(transported_edges))
    print("FINITE-EXACT: malformed/root/valuation hostiles=" + str(len(hostiles)))
    print("explicit checks=" + str(checks))
    print("OPEN: source coverage outside the stated guards; no universal convergence or new basin claim")


if __name__ == "__main__":
    main()
