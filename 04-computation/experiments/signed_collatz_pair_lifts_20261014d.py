"""Sheet-labelled Collatz inverse fibres and supplied terminal certificates.

Run with python -B [-O] this_file.py. No production API discovers an orbit.
"""
from dataclasses import dataclass, replace
from hashlib import sha256
from itertools import product
import json


def integer(x):
    if type(x) is not int:
        raise TypeError("exact integer required")
    return x


def sheet(b):
    integer(b)
    if b not in (-1, 1):
        raise ValueError("sheet must be -1 or +1")
    return b


def odd(n):
    integer(n)
    if not n or not n % 2:
        raise ValueError("nonzero odd source required")
    return n


def word(w):
    if type(w) is not tuple:
        raise TypeError("tuple word required")
    for a in w:
        integer(a)
        if a < 1:
            raise ValueError("positive valuations required")
    return w


def v2(n):
    integer(n)
    if not n:
        raise ValueError("valuation at zero undefined")
    n = abs(n)
    return (n & -n).bit_length() - 1


def shortcut(n, b):
    integer(n)
    sheet(b)
    return n // 2 if n % 2 == 0 else (3 * n + b) // 2


def step(n, b):
    odd(n)
    sheet(b)
    z = 3 * n + b
    a = v2(z)
    return z // (1 << a), a


def carrier(w):
    word(w)
    p, q, carry = 1, 1, 0
    for a in w:
        p, carry, q = 3 * p, 3 * carry + q, q << a
    return p, q, carry


def replay(n, b, w):
    odd(n)
    sheet(b)
    word(w)
    for a in w:
        n, actual = step(n, b)
        if a != actual:
            raise ValueError("non-native valuation")
    return n


@dataclass(frozen=True)
class Port:
    sheet: int
    source: int
    exponent: int
    endpoint: int


def audit_port(p):
    if type(p) is not Port:
        raise TypeError("exact Port required")
    sheet(p.sheet)
    odd(p.source)
    odd(p.endpoint)
    integer(p.exponent)
    if p.exponent < 1 or step(p.source, p.sheet) != (p.endpoint, p.exponent):
        raise ValueError("port fails actual edge authentication")
    return p


def port_at(q, a):
    """The unique port at valuation a over both sheets with endpoint q."""
    odd(q)
    integer(a)
    if q % 3 == 0 or a < 1:
        raise ValueError("endpoint must be a ternary unit; exponent positive")
    z = (1 << a) * q
    b = 1 if z % 3 == 1 else -1
    return audit_port(Port(b, (z - b) // 3, a, q))


def lift(p, count):
    audit_port(p)
    integer(count)
    if count < 0:
        raise ValueError("nonnegative lift count required")
    power = 1 << count
    parity = -1 if count % 2 else 1
    n = power * p.source + p.sheet * (power - parity) // 3
    result = Port(parity * p.sheet, n, p.exponent + count, p.endpoint)
    return audit_port(result)


def mersenne_ports(q):
    """Complete positive Mersenne-source intersection of the fibre over q.

    Return pairs (E, authenticated_port), where source=2**E-1 and E>=1.
    This is a finite exact classifier, not a bounded exponent search.
    """
    odd(q)
    if q < 1:
        raise ValueError("positive endpoint required")
    if q == 1:
        return ((1, port_at(q, 2)), (1, port_at(q, 1)), (2, port_at(q, 3)))
    if (q + 1) % 3:
        return ()
    power = (q + 1) // 3
    if power < 2 or power & (power - 1):
        return ()
    r = power.bit_length() - 1
    return ((r + 1, port_at(q, 1)), (r + 2, port_at(q, 2)))


def canonical_cycle(cycle):
    if type(cycle) is not tuple or not cycle:
        raise TypeError("nonempty tuple cycle required")
    for n in cycle:
        odd(n)
    if len(set(cycle)) != len(cycle):
        raise ValueError("cycle has repeated vertices")
    i = min(range(len(cycle)), key=lambda j: (abs(cycle[j]), cycle[j]))
    return cycle[i:] + cycle[:i]


def audit_cycle(b, cycle):
    sheet(b)
    if canonical_cycle(cycle) != cycle:
        raise ValueError("cycle is not canonically marked")
    for x, y in zip(cycle, cycle[1:] + cycle[:1]):
        if step(x, b)[0] != y:
            raise ValueError("cycle edge failed")
    return cycle


@dataclass(frozen=True)
class Terminal:
    sheet: int
    source: int
    tail: tuple
    cycle: tuple


def audit_terminal(t):
    if type(t) is not Terminal:
        raise TypeError("exact Terminal required")
    sheet(t.sheet)
    odd(t.source)
    word(t.tail)
    audit_cycle(t.sheet, t.cycle)
    x = t.source
    for a in t.tail:
        if x in t.cycle:
            raise ValueError("terminal padding before final entry")
        x, actual = step(x, t.sheet)
        if a != actual:
            raise ValueError("non-native terminal word")
    if x not in t.cycle:
        raise ValueError("tail did not reach supplied terminal cycle")
    return t


def discharge_port(p, endpoint_terminal):
    """Only a destination-sheet certificate may cross this seam."""
    audit_port(p)
    audit_terminal(endpoint_terminal)
    t = endpoint_terminal
    if (p.sheet, p.endpoint) != (t.sheet, t.source):
        raise ValueError("sheet or endpoint mismatch")
    tail = () if p.source in t.cycle else (p.exponent,) + t.tail
    return audit_terminal(Terminal(p.sheet, p.source, tail, t.cycle))


def reflected(t):
    audit_terminal(t)
    return audit_terminal(Terminal(-t.sheet, -t.source, t.tail,
                                   canonical_cycle(tuple(-x for x in t.cycle))))


# Exact discovery was bounded at 10,000 odd steps for these four records only.
# These supplied words, not a discovery loop, are the production evidence.
MINUS17 = (17, 25, 37, 55, 41, 61, 91)
FIXED = (
    Terminal(1, 1473, (2,2,3,1,1,3,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4), (1,)),
    Terminal(-1, 1473, (1,1,1,1,1,8,3,1,1,1,3), MINUS17),
    Terminal(1, 26861, (3,2,1,4,3,3,1,2,2,3,1,1,1,1,4,1,2,1,1,2,1,1,1,2,3,1,1,2,1,2,1,1,1,1,1,3,1,1,1,4,2,2,4,3,1,1,5,4), (1,)),
    Terminal(-1, 26861, (1,3,1,2,1,2,1,1,1,3,2,1,2,1,1,1,1,1,4,1,3,2,2,1,3,1,2,4,3,1,1,2,1,1,1,3,2,1,3,2,3,2,1,3,7), (1,)),
)


def supplied_terminal(source, b):
    odd(source)
    sheet(b)
    for t in FIXED:
        if (t.source, t.sheet) == (source, b):
            return audit_terminal(t)
    raise ValueError("no supplied terminal record for this marked source")


def _discover_bounded(source, b, cap):
    """Test-only independent bounded exploration, never called by production."""
    x, seen, states, vals = source, {}, [], []
    for _ in range(cap):
        if x in seen:
            j = seen[x]
            return tuple(vals[:j]), canonical_cycle(tuple(states[j:]))
        seen[x] = len(states)
        states.append(x)
        # Independent division loop rather than step()/v2().
        x = 3 * x + b
        a = 0
        while x % 2 == 0:
            a += 1
            x //= 2
        vals.append(a)
    return None


def main():
    checks = 0
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            raise RuntimeError(label)
    def rejects(f, label):
        try:
            f()
        except (TypeError, ValueError):
            check(True, label)
        else:
            check(False, label)

    for n in range(-999, 1000, 2):
        up, ap = step(n, 1)
        um, am = step(n, -1)
        check(shortcut(n, 1) - shortcut(n, -1) == 1, "consecutive shortcut")
        check((shortcut(n, 1) % 3, shortcut(n, -1) % 3) == (2, 1), "shortcut residues")
        check(min(ap, am) == 1 and max(ap, am) >= 2, "paired valuations")
        check((1 << ap) * up - (1 << am) * um == 2, "paired carry")
        check(up % 3 == (-1) ** ap % 3 and um % 3 == -(-1) ** am % 3, "odd residues")
        for b in (-1, 1):
            q, a = step(n, b)
            check(step(-n, -b) == (-q, a), "negation")
            p = audit_port(Port(b, n, a, q))
            for m in range(9):
                z = lift(p, m)
                check(z == port_at(q, a + m), "complete fibre lift")
            check(lift(p, 2).source == 4 * n + b, "two lifts are sibling")
    for q in range(-99, 100, 2):
        if q % 3:
            for a in range(1, 17):
                p = port_at(q, a)
                check(p.source * q > 0, "port sign")
                check(lift(p, 1).source == 2*p.source+p.sheet, "adjacent port")
    for length in range(1, 6):
        for w in product(range(1, 5), repeat=length):
            p, q, carry_ = carrier(w)
            check(carry_ % 3 == pow(2, sum(w[:-1]), 3), "carry unit")
            check(2 * carry_ % p != 0, "same word crosssheet obstruction")

    plus_cycles = ((1,), (-1,), (-5, -7), tuple(-x for x in MINUS17))
    cycle_rows = []
    for cycle in plus_cycles:
        audit_cycle(1, cycle)
        w = tuple(step(x, 1)[1] for x in cycle)
        p, q, carry_ = carrier(w)
        check((q - p) * cycle[0] == carry_, "cycle carry equation")
        t = Terminal(1, cycle[0], (), cycle)
        check(reflected(reflected(t)) == t, "terminal negation involution")
        expanded = []
        for x, a in zip(cycle, w):
            for _ in range(a):
                expanded.append(x)
                x = shortcut(x, 1)
        check(len(set(expanded)) == sum(w), "shortcut simple cycle")
        check(all(shortcut(x, 1) == y for x,y in zip(expanded, expanded[1:]+expanded[:1])), "shortcut cycle")
        cycle_rows.append(dict(odd_cycle=cycle, odd_word=w, shortcut_cycle=expanded,
                               denominator=q-p, carry=carry_))

    fixed_rows = []
    for t in FIXED:
        audit_terminal(t)
        check(_discover_bounded(t.source, t.sheet, 10000) == (t.tail, t.cycle), "independent supplied trace")
        check(reflected(reflected(t)) == t, "signed supplied reflection")
        p, q, carry_ = carrier(t.tail)
        endpoint = replay(t.source, t.sheet, t.tail)
        check(p*t.source+t.sheet*carry_ == q*endpoint, "full supplied affine carrier")
        fixed_rows.append(dict(source=t.source, sheet=t.sheet, odd_rank=len(t.tail),
                               halving_cost=sum(t.tail), endpoint=endpoint, cycle=t.cycle,
                               word=t.tail, word_sha256=sha256(bytes(t.tail)).hexdigest()))

    root_plus = Terminal(1, 1, (), (1,))
    root_minus = Terminal(-1, 1, (), (1,))
    five_plus = Terminal(1, 5, (4,), (1,))
    five_minus = Terminal(-1, 5, (), (5, 7))
    for j in range(32):
        n, m = (10*4**j-1)//3, (20*4**j+1)//3
        pp, pm = port_at(5, 2*j+1), port_at(5, 2*j+2)
        check((pp.sheet, pp.source, pm.sheet, pm.source) == (1,n,-1,m), "hostile family sources")
        check(discharge_port(pp, five_plus).cycle == (1,), "plus family ROOT")
        check(discharge_port(pm, five_minus).cycle == (5,7), "minus family nonROOT")
        rejects(lambda: discharge_port(pm, five_plus), "sheet-blind seam rejected")
        for a, t in ((2*j+2, root_plus), (2*j+1, root_minus)):
            z = discharge_port(port_at(1,a),t)
            check(z.tail == (() if z.source == 1 else (a,)), "ROOT padding removed")

    # A single numeric endpoint with independently supplied proofs on both sheets
    # grounds every member of its combined inverse fibre, without sign switches.
    for a in range(1, 65):
        p = port_at(26861, a)
        t = discharge_port(p, supplied_terminal(26861, p.sheet))
        check(t.cycle == (1,), "two-sheet grounded fibre")
        check(len(t.tail) == (49 if a % 2 else 46), "grounded fibre rank")
        check(sum(t.tail) == a + (91 if a % 2 else 86), "grounded fibre cost")

    check(mersenne_ports(26861) == (), "grounded fibre has no Mersenne source")
    check((26861 + 1) // 3 == 8954 and 8954 & (8954 - 1) != 0,
          "exact power-of-two obstruction")
    check(tuple((e,p.sheet,p.exponent) for e,p in mersenne_ports(1)) ==
          ((1,1,2),(1,-1,1),(2,-1,3)), "small Mersenne exceptions")
    for e in range(1, 129):
        for b in (-1, 1):
            n = (1 << e) - 1
            q, a = step(n, b)
            check((e,Port(b,n,a,q)) in mersenne_ports(q), "Mersenne classifier completeness")
    for q in range(1, 1000, 2):
        expected = []
        for e in range(1, 12):
            for b in (-1, 1):
                n = (1 << e) - 1
                if step(n,b)[0] == q:
                    expected.append((e,b,n,step(n,b)[1]))
        actual = [(e,p.sheet,p.source,p.exponent) for e,p in mersenne_ports(q)]
        check(sorted(expected) == sorted(actual), "independent bounded Mersenne control")

    # Same integer endpoint can carry incompatible futures on different sheets.
    check(step(3,1)[0] == step(7,-1)[0] == 5, "common integer endpoint")
    check(step(5,1)[0] == 1 and step(5,-1)[0] == 7, "different continuation")
    check(step(1,1)[1] == 2 and step(1,-1)[1] == 1, "ROOT shortcut clock distinction")
    good = port_at(5,1)
    for bad in (replace(good, sheet=True), replace(good, source=3.0),
                replace(good, exponent=True), replace(good, endpoint=5.0),
                replace(good, endpoint=7)):
        rejects(lambda bad=bad: audit_port(bad), "malformed port")
    for bad in (replace(five_plus, sheet=True), replace(five_plus, source=5.0),
                replace(five_plus, tail=(True,)), replace(five_plus, tail=[4]),
                replace(five_plus, cycle=(True,)), replace(five_plus, tail=(4,2))):
        rejects(lambda bad=bad: audit_terminal(bad), "malformed terminal")
    for args in ((3,1), (1,0), (True,2), (1,2.0)):
        rejects(lambda args=args: port_at(*args), "invalid fibre")
    rejects(lambda: lift(good,-1), "negative lift")
    rejects(lambda: discharge_port(port_at(1,2),five_plus), "wrong endpoint")
    for bad in (True, 1.0, -1, 0, 2):
        rejects(lambda bad=bad: mersenne_ports(bad), "invalid Mersenne endpoint")

    result = dict(status="PROVED elementary identities; FINITE-EXACT controls; no universal basin classification",
                  checks=checks, identity_universe="all odd -999..999, both sheets, lifts0..8",
                  word_universe="length1..5, letters1..4", known_plus_nonzero_cycles=cycle_rows,
                  supplied_terminal_records=fixed_rows,
                  fully_grounded_fibre=dict(endpoint=26861, exponent="a>=1",
                      sheet="(-1)^(a+1)", source="(26861*2^a-(-1)^(a+1))/3",
                      odd_rank="49 for odd a; 46 for even a",
                      halving_cost="a+91 for odd a; a+86 for even a",
                      positive_mersenne_source_intersection=[]),
                  infinite_hostile_family="plus (10*4^j-1)/3 reaches1; minus (20*4^j+1)/3 enters(5,7), j>=0",
                  production_orbit_discovery=False)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
