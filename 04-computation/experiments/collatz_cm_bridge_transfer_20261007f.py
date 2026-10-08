"""Protected-source sibling bridges; elementary mathematics, not modularity.

Production functions never discover a ROOT path. main() uses frozen supplied
words and finite independent arithmetic controls. Run with Python -B [-O].
"""
from dataclasses import dataclass, replace
from itertools import product


def natural(x, name="integer", minimum=0):
    if type(x) is not int or x < minimum:
        raise ValueError(name)
    return x


def odd(n):
    natural(n, "positive odd source", 1)
    if n % 2 == 0:
        raise ValueError("positive odd source")
    return n


def letters(w):
    if type(w) is not tuple or not w:
        raise ValueError("nonempty exact tuple")
    for a in w:
        natural(a, "positive valuation", 1)
    return w


def carrier(w):
    letters(w)
    p, q, b = 1, 1, 0
    for a in w:
        p, q, b = 3*p, q*(1 << a), 3*b+q
    return p, q, b


def step(n):
    odd(n)
    z = 3*n+1
    a = (z & -z).bit_length()-1
    return z >> a, a


def replay(n, w, root=False):
    odd(n)
    if type(w) is not tuple:
        raise ValueError("exact tuple")
    for a in w:
        natural(a, "positive valuation", 1)
        if n == 1:
            raise ValueError("earlier ROOT")
        n, actual = step(n)
        if actual != a:
            raise ValueError("wrong native valuation")
    if root and n != 1:
        raise ValueError("not ROOT")
    return n


def sibling(n, k):
    odd(n)
    natural(k, "sibling depth")
    return (1 << (2*k))*n+((1 << (2*k))-1)//3


def lift_phase(n, w):
    """Unique k mod3^len(w) making F_w^{-1}(S^k(n)) integral."""
    odd(n)
    p, q, b = carrier(w)
    modulus = 3*p
    target = ((q+3*b)*pow(q*(3*n+1), -1, modulus)) % modulus
    phase, period = 0, 1
    for _ in w:
        next_modulus = 9*period
        hits = [phase+j*period for j in range(3)
                if pow(4, phase+j*period, next_modulus)
                == target % next_modulus]
        if len(hits) != 1:
            raise ArithmeticError("principal-unit lift")
        phase = hits[0]
        period *= 3
    return phase, period


def budget_cap(n, w):
    """Largest nonnegative k that can pay h<n; -1 means none."""
    odd(n)
    p, q, b = carrier(w)
    cap = (3*p*n+q+3*b-1)//(q*(3*n+1))
    return (cap.bit_length()-1)//2 if cap >= 1 else -1


def paying_lifts(n, w):
    """Complete finite list of positive sibling depths for this fixed n,w."""
    odd(n)
    if n == 1:
        raise ValueError("source must exceed ROOT")
    phase, period = lift_phase(n, w)
    first = phase if phase else period
    return tuple(range(first, budget_cap(n, w)+1, period))


@dataclass(frozen=True)
class Receipt:
    source: int
    child: int
    lift: int
    inverse_word: tuple
    source_word: tuple
    child_word: tuple
    endpoint: int


def make_receipt(n, w, k):
    odd(n)
    if n == 1:
        raise ValueError("source must exceed ROOT")
    letters(w)
    natural(k, "positive sibling depth", 1)
    p, q, b = carrier(w)
    value = q*sibling(n, k)-b
    if value % p:
        raise ValueError("inverse guard")
    h = value//p
    if not 1 < h < n:
        raise ValueError("unpaid child")
    endpoint, a = step(n)
    child_word = w+(a+2*k,)
    if replay(h, w) != sibling(n, k):
        raise ArithmeticError("inverse identity")
    if replay(h, child_word) != endpoint:
        raise ArithmeticError("common future")
    return Receipt(n, h, k, w, (a,), child_word, endpoint)


def audit_receipt(r):
    if type(r) is not Receipt:
        raise ValueError("exact Receipt")
    for n in (r.source, r.child, r.endpoint):
        odd(n)
    natural(r.lift, "positive sibling depth", 1)
    letters(r.inverse_word)
    letters(r.source_word)
    letters(r.child_word)
    if r != make_receipt(r.source, r.inverse_word, r.lift):
        raise ValueError("noncanonical receipt")
    return r


def discharge(r, child_root_word):
    """Consume only a supplied, exact, first-hit child ROOT certificate."""
    audit_receipt(r)
    replay(r.child, child_root_word, root=True)
    if child_root_word[:len(r.child_word)] != r.child_word:
        raise ValueError("incompatible supplied child proof")
    result = r.source_word+child_root_word[len(r.child_word):]
    replay(r.source, result, root=True)
    return result


def ternary_fuel(n):
    odd(n)
    z, depth = 8*n+3, 0
    while z % 3 == 0:
        z //= 3
        depth += 1
    return depth


def family_receipt(n, depth=None):
    """Maximum depth is the smallest child within this declared family."""
    odd(n)
    fuel = ternary_fuel(n)
    if depth is None:
        depth = fuel
    natural(depth, "paid depth", 6)
    if depth > fuel:
        raise ValueError("missing ternary fuel")
    return make_receipt(n, (1,)*(depth-1)+(2,), 1)


ROOT191 = (1,1,1,1,1,2,4,3,3,3,1,2,3,4)
ROOT1215 = (1,1,1,1,1,2,3,4,1,1,2,1,1,3,1,1,3,1,2,6,1,1,1,1,2,2,1,
            2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4)


def main():
    checks = 0
    def check(condition):
        nonlocal checks
        checks += 1
        if not condition:
            raise AssertionError(checks)
    def rejects(f):
        try:
            f()
        except (ValueError, TypeError):
            check(True)
        else:
            check(False)

    # Exhaustive independent phase scan for lengths <=3, letters1..4.
    for length in range(1, 4):
        for w in product(range(1, 5), repeat=length):
            p, q, b = carrier(w)
            for n in range(3, 66, 2):
                phase, period = lift_phase(n, w)
                direct = tuple(k for k in range(p)
                               if (q*sibling(n, k)-b) % p == 0)
                check(direct == (phase,) and period == p)
                expected = []
                # Native depths in a complete period plus all potentially paid
                # depths; the independent test keeps the strict inequality.
                for k in range(1, max(p+1, budget_cap(n, w)+2)):
                    value = q*sibling(n, k)-b
                    if value % p == 0 and 0 < value//p < n:
                        expected.append(k)
                check(paying_lifts(n, w) == tuple(expected))
                for k in expected:
                    r = make_receipt(n, w, k)
                    check(replay(r.child, r.child_word) == replay(n, r.source_word))

    # Independent all-height progression controls; no proof search.
    for t in range(512):
        n = 273+1458*t
        r = family_receipt(n, 6)
        check(r.child == 191+1024*t)
        check(n % 3 == 0 and r.child < n)
        check(replay(r.child, r.inverse_word) == 4*n+1)
        maximum = family_receipt(n)
        check(maximum.child <= r.child)
        for depth in range(1, ternary_fuel(n)+1):
            w = (1,)*(depth-1)+(2,)
            p, q, b = carrier(w)
            h = (q*(4*n+1)-b)//p
            check(replay(h, w) == 4*n+1)
            check((h < n) == (depth >= 6))
        if t % 2:
            check(r.source_word == (1,) and r.endpoint > n)
            check(r.child_word == (1,1,1,1,1,2,3))
    for depth in range(1, 25):
        phase = (11*pow(3**depth, -1, 16)) % 16
        for j in range(12):
            t = phase+16*j
            n = (3**depth*t-3)//8
            h = 2**depth*t-1
            check(n > 0 and n % 2 == 1 and h % 2 == 1)
            check(replay(h, (1,)*(depth-1)+(2,)) == 4*n+1)
            check((h < n) == (depth >= 6))

    first = family_receipt(273, 6)
    check(discharge(first, ROOT191) == (2,3,3,3,1,2,3,4))
    growing = family_receipt(1731, 6)
    rooted = discharge(growing, ROOT1215)
    check(replay(1731, rooted, root=True) == 1 and len(rooted) == 53)
    check(family_receipt(273).child == 127)
    # A protected source may not be swapped for a neighbour satisfying a guard.
    rejects(lambda: family_receipt(275, 6))
    rejects(lambda: discharge(first, ROOT1215))
    rejects(lambda: audit_receipt(replace(first, source=True)))
    rejects(lambda: audit_receipt(replace(first, child=float(first.child))))
    rejects(lambda: audit_receipt(replace(first, lift=True)))
    rejects(lambda: audit_receipt(replace(first, source_word=(2.0,))))
    rejects(lambda: replay(True, (), root=True))
    rejects(lambda: replay(1.0, (), root=True))
    rejects(lambda: replay(1, (2,), root=True))
    rejects(lambda: paying_lifts(1, (2,)))
    rejects(lambda: family_receipt(273, True))
    rejects(lambda: lift_phase(273, (True,)))
    rejects(lambda: lift_phase(273, []))
    rejects(lambda: make_receipt(7, (1,), 1))
    rejects(lambda: make_receipt(273, (1,)*5+(2,), 2))
    # Exact all-depth obstruction witness: the whole small basin is closed.
    small = {1, 3, 5}
    for n in small:
        check(step(n)[0] in small)
    check(min(sibling(7, k) for k in range(10)) == 7)
    check(replay(7, (1,1,2,3,4), root=True) == 1)

    print("Protected-source CM-inspired bridge: elementary exact controls")
    print("Phase universe:84 words(length1..3,letters1..4),32 odd sources3..65")
    print("Native k phase: unique modulo3^length; source-budget search is finite and complete")
    print("Family: n=273+1458t -> h=191+1024t, t>=0; density among odds=1/729")
    print("Growing subfamily: n=1731+2916t -> h=1215+2048t")
    print("Exact payment iff depth>=6; maximal depth=v3(8n+3) selects smallest family child")
    print("Supplied ROOT discharges:273 (8 edges),1731 (53 edges)")
    print("All-depth hostile:7 has no smaller sibling-inverse child; its ROOT word is supplied")
    print("Paper headline trusted for architecture only; universal source coverage remains OPEN")
    print(f"PASS {checks} exact checks")


if __name__ == "__main__":
    main()
