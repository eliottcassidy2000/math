"""Exact prime and Hamiltonian certificates for the nine-to-twenty-five window."""
from itertools import product
from math import isqrt
from pathlib import Path
import sys

RESULTS = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(RESULTS / "chanpak_prime_frontier_selector_20261118"))
sys.path.insert(0, str(RESULTS / "chanpak_tiling_viewpoint_20261123"))
from join_loads import analyze
from precedence_bank import count_and_order
sys.path.insert(0, str(RESULTS / "chanpak_gap_star_ports_20261122"))
from replay import count_from_profile, flipped, profile, star_rows, strong
from check_perfect_band import remove_common_extremes

P = (5, 7, 11, 13, 17, 19, 23)
A = (2, 4, 8, 10, 14)
ORIGINAL = (11, 13, 17, 19, 23)
LATE = {
    32397: [(23,35060114), (31,28983694), (3,216663367), (31,24964507),
            (5,143055378), (2,515990256), (7,122641422), (87,8551115), (3,216671559)],
    213237: [(2,19162820870), (3,13432939772), (2,21753521944), (59,550949287),
             (11,3936528912), (43,701908019), (2027,16493519), (2,21249748028), (2,16583567427)],
    254577: [(2,28694721470), (2,30923579112), (507,124629424), (25,2273623836),
             (3,14470826882), (2,32116025056), (2,30558865544), (2,26347876946), (2,26873988458)],
}


# (square root, lower offset): (lower core, flipped pair, tags, left and right addresses)
FLIPS = {
    (32397, -2): ((12, 81, 74, 162, 45, 135, 57, 87), (4, 7), 5, 48, 2952),
    (32397, 0): ((0, 5, 129, 7, 15, 31, 63, 123), (2, 4), 34, 0, 3126502),
    (32397, 4): ((0, 33, 67, 23, 7, 93, 155, 63), (4, 5), 9, 0, 2403024),
    (32397, 6): ((0, 33, 67, 23, 7, 93, 155, 63), (4, 5), 9, 0, 2646553),
    (32397, 8): ((2, 32, 19, 23, 67, 157, 47, 95), (4, 5), 36, 0, 1949061),
    (213237, -2): ((0, 5, 49, 39, 75, 147, 47, 95), (4, 5), 5, 0, 168808150),
    (213237, 0): ((0, 1, 19, 23, 195, 31, 47, 111), (4, 5), 68, 0, 130147542),
    (213237, 4): ((0, 1, 19, 71, 43, 207, 151, 31), (5, 7), 7, 0, 195134908),
    (213237, 8): ((12, 81, 74, 162, 45, 135, 57, 87), (3, 6), 36, 15, 404645),
    (254577, 0): ((4, 9, 146, 37, 11, 23, 63, 123), (2, 4), 32, 0, 160302552),
    (254577, 6): ((0, 1, 3, 23, 39, 143, 63, 95), (5, 6), 21, 0, 553127135),
    (254577, 12): ((0, 1, 19, 23, 195, 31, 47, 111), (4, 5), 6, 0, 399719416),
}


def crt(pairs):
    value, period = 0, 1
    for target, modulus in pairs:
        value += period * ((target-value) * pow(period, -1, modulus) % modulus)
        period *= modulus
    return value


def address(bits, old_bank=False):
    # q=129*u*37^i*751^j. Choose the independent odd parts of i by CRT.
    b11, b13, b17, b19, b23 = bits
    if old_bank and b13 and b17:
        return None
    u = 17 if old_bank and b13 else (41 if b13 else 1)
    # Sign-class logs at 11, 13, 17, 19, 23 of the fixed prefactor.
    logs = {1: (3,0,3,2,10), 41:(6,1,6,6,11), 17:(7,2,None,3,17)}[u]
    r11, r13, r17, r19, r23 = logs
    target13 = 4 if b13 else 1
    j9 = (r19-4*b19) % 9
    if u == 17:
        j = j9
        i_even = (0, 2)
    else:
        j = crt([(j9,9), ((b17-r17-1) % 2,2)])
        i_even = ((b17-r17-j) % 8,8)
    i = crt([i_even, ((target13-r13-j) % 3,3),
             ((1+2*b11-r11-3*j)*3 % 5,5), ((r23+6*j-5*b23) % 11,11)])
    return u, i, j


def verify_address(bits, u, i, j):
    for p in (2,3)+P:
        q = 129 * u * pow(37,i,p) * pow(751,j,p) % p
        blocked = [a for a in A if (q*q+a) % p == 0]
        expected = [A[ORIGINAL.index(p)]] if p in ORIGINAL and bits[ORIGINAL.index(p)] else []
        assert blocked == expected


def first_patterns(limit):
    # The largest tested value has sqrt less than limit+1. Trial division is exact.
    sieve = bytearray([1]) * (limit+2)
    sieve[:2] = bytes(2)
    for p in range(2,isqrt(limit+1)+1):
        if sieve[p]:
            sieve[p*p::p] = bytes((limit+1-p*p)//p + 1)
    primes = [p for p, yes in enumerate(sieve) if yes]
    def isprime(n):
        for p in primes:
            if p*p > n:
                return True
            if n % p == 0:
                return False
        return True
    return [q for q in range(3, limit+1, 2)
            if q % 30 in (3,27) and all(isprime(q*q+a) for a in A)]


def main():
    assert count_and_order((118,120,58,113,64,80,4)) == (41,list(range(7)))
    # Ordinary classes through six, using the existing complete small-class bank.
    possible = set()
    data = RESULTS.parents[1] / "session212/cycle-levels"
    for n in range(1,7):
        for line in (data / f"classes{n}.rows").read_text().splitlines():
            rows = tuple(map(int,line.split()))
            try:
                h, _ = count_and_order(rows)
            except AssertionError:
                continue
            possible.add(h)
    assert possible == {1,5,9,17,37}
    print("Heights with transitive precedence through order six", sorted(possible))
    for old, number in [(True,24),(False,32)]:
        success = 0
        for bits in product(range(2), repeat=5):
            entry = address(bits,old)
            if entry is not None:
                verify_address(bits, *entry)
                success += 1
        assert success == number
        print("Exact addresses constructed", success)
    base = 37**6*751
    for p in P:
        assert pow(base,720,p) == 1
    lifts = [(pow(base,720,p*p)-1)//p for p in ORIGINAL]
    assert lifts == [10,10,15,12,12]
    print("Independent prime-power lift digits", lifts)
    roots = first_patterns(max(LATE))
    assert roots == [3,*LATE]
    for q, joins in LATE.items():
        orders = []
        for a, (b,e) in zip(range(-2,15,2), joins):
            obj = analyze(b,e)
            assert obj["H"] == q*q+a and len(obj["rows"]) == obj["order"]
            orders.append(obj["order"])
        print("Prime constellation root",q,"nine literal tournament orders",orders)

    for (q, lower), (core, edge, tags, x, y) in FLIPS.items():
        cores = (core, flipped(core, *edge))
        counts = [count_from_profile(profile(c,tags), x, y) for c in cores]
        assert counts == [q*q+lower,q*q+lower+2]
        boards = remove_common_extremes(*(star_rows(c,tags,x,y) for c in cores))
        assert all(strong(t) for t in boards) and sum((a^b).bit_count() for a,b in zip(*boards)) == 2
    print("Strong single-arc witnesses at lower offsets",
          {q:[lo for root,lo in FLIPS if root==q] for q in LATE})


if __name__ == "__main__":
    main()
