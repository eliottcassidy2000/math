"""Typed prime-prefix and divisor-word orbit audit. Standard library only."""

from itertools import combinations, product
from math import gcd, isqrt, prod
from pathlib import Path


def check(condition, message="check failed"):
    if not condition:
        raise RuntimeError(message)


def prime(n):
    if n < 2:
        return False
    if n % 2 == 0:
        return n == 2
    return all(n % d for d in range(3, isqrt(n) + 1, 2))


def factor(n):
    out = {}
    d = 2
    while d*d <= n:
        while n % d == 0:
            out[d] = out.get(d, 0) + 1
            n //= d
        d += 1
    if n > 1:
        out[n] = out.get(n, 0) + 1
    return out


def divisors(n):
    out = set()
    for d in range(1, isqrt(n) + 1):
        if n % d == 0:
            out.update((d, n//d))
    return out


def divisor_sets(n):
    D = divisors(n) - {1, n}
    S = {d for d in D if all(a == 1 for a in factor(d).values())}
    U = {d for d in D if prime(d)}
    return D, S, U


def is_square(n):
    return isqrt(n)**2 == n


def action(p, q, r):
    """C3 action on proper divisors of p^2 q r, via subset/extra chart."""
    check(len({p, q, r}) == 3 and all(prime(x) for x in (p, q, r)))
    n = p*p*q*r
    D, S, U = divisor_sets(n)
    rename = {p: q, q: r, r: p}
    g = {d: prod(rename[z] for z in factor(d)) for d in S}
    lift = {p: p*p, q: p*p*q, r: p*p*r}
    g.update({lift[z]: lift[rename[z]] for z in (p, q, r)})
    check(set(g) == D and set(g.values()) == D)
    check(all(g[g[g[d]]] == d for d in D))
    check({d for d in D if g[d] == d} == {p*q*r})
    check({g[d] for d in S} == S)
    return n, D, S, U, g


def orbit_partition(D, S, g, k):
    out = {}
    for w in product(sorted(D), repeat=k):
        if len(set(w)) == 1 and w[0] in S:
            continue
        w1 = tuple(g[d] for d in w)
        w2 = tuple(g[d] for d in w1)
        check(len({w, w1, w2}) == 3)
        root = min(w, w1, w2)
        out.setdefault(root, set()).add(w)
    check(all(len(v) == 3 for v in out.values()))
    return out


def main():
    lines = []

    def report(value):
        lines.append(str(value))

    report("DUCK PRIME/PREFIX/DIVISOR AUDIT 2026-09-25")
    report("Universe: first30 retained prime prefixes for three specified filters; exact threshold196 witnesses; divisor examples and literal word orbits at k1..5 for N60 and k1..3 for N84; first decimal candidates k2..9. Standard library; explicit checks survive -O.")
    P = [p for p in range(2, 300) if prime(p)]
    filters = [set(), {11}, {2, 3, 11}]
    all_sums = {}
    for E in filters:
        retained = [p for p in P if p not in E]
        sums = [sum(retained[:k]) for k in range(1, 31)]
        squares = [(k, v) for k, v in enumerate(sums, 1) if is_square(v)]
        below = max((k, v) for k, v in enumerate(sums, 1) if v < 196)
        above = next((k, v) for k, v in enumerate(sums, 1) if v > 196)
        check(196 not in sums and above[0] == below[0] + 1)
        report(f"Filter omit{sorted(E)}: retained prefixes straddle196 at {below},{above}; square sums among first30={squares}")
        report(f"  k8..15: (k,cutoff-then-delete,first-k-retained)={[(k,sum(p for p in P[:k] if p not in E),sum(retained[:k])) for k in range(8,16)]}")
        all_sums[tuple(sorted(E))] = sums
    check(all_sums[()][10:12] == [160, 197])
    check(all_sums[(11,)][10:12] == [186, 227])
    check(all_sums[(2, 3, 11)][8:10] == [181, 222])
    # Finite difference-of-squares certificates for the two deletion sums.
    for c in (11, 16):
        roots = []
        for a in divisors(c):
            b = c//a
            if a <= b and (a+b) % 2 == 0:
                high, low = (a+b)//2, (b-a)//2
                check(high*high-low*low == c)
                roots.append((high*high, low*low))
        report(f"All nonnegative square differences high-low={c}: {sorted(roots)}")
    check(sum(P[:5]) == 28 and sum(P[:6]) == 41)
    report("Beyond cutoff11, subtracting11 or16 cannot leave BOTH original and filtered prime-prefix sums square: candidate highs36 or{16,25} miss the required prefix.")

    for n in (4, 8, 10, 60, 84, 196):
        D, S, U = divisor_sets(n)
        report(f"Value{n}: factors={factor(n)}, (F,S,U)=({len(D)},{len(S)},{len(U)}), balanced={len(D)==len(S)+len(U)}")
    check(tuple(len(x) for x in divisor_sets(10)) == (2, 2, 2))
    check(tuple(len(x) for x in divisor_sets(60)) == (10, 7, 3))
    additive = [(a, 10-a) for a in range(1, 6)]
    prime_additive = [(a, b) for a, b in additive if prime(a) and prime(b)]
    proper_factors = [(d, 60//d) for d in range(2, isqrt(60)+1) if 60 % d == 0]
    check(prime_additive == [(3, 7), (5, 5)])
    check(proper_factors == [(2,30), (3,20), (4,15), (5,12), (6,10)])
    report(f"10 additive pairs={additive}; prime pairs={prime_additive}, distinct-only leaves(3,7).")
    report(f"60 has one prime factorization but five unordered proper factor groupings={proper_factors}.")

    counts = []
    n, D, S, U, g = action(2, 3, 5)
    report(f"N60 literal divisor permutation={sorted(g.items())}; squarefree subset={sorted(S)}")
    previous = None
    for k in range(1, 6):
        partition = orbit_partition(D, S, g, k)
        count = len(partition)
        check(count == (10**k-7)//3)
        # Independent Burnside count using only fixed letters / forbidden words.
        fixed = sum(g[d] == d for d in D)
        check(fixed == 1 and (10**k + 2*fixed**k)//3 - 3 == count)
        if previous is not None:
            children = {root: 0 for root in previous}
            births = 0
            for root in partition:
                prefix = root[:-1]
                if len(set(prefix)) == 1 and prefix[0] in S:
                    births += 1
                else:
                    p1 = tuple(g[d] for d in prefix)
                    p2 = tuple(g[d] for d in p1)
                    children[min(prefix, p1, p2)] += 1
            check(set(children.values()) == {10} and births == 21)
            check(count == 10*len(previous)+21)
        counts.append(count)
        previous = partition
    check(counts == [1,31,331,3331,33331])
    report(f"Literal free-word orbit counts at k1..5={counts}; all orbits size3; each prior orbit has10 children, with21 new orbits at every next level.")
    _, D2, S2, U2, g2 = action(2, 3, 7)
    check([len(orbit_partition(D2,S2,g2,k)) for k in range(1,4)] == counts[:3])
    report("Independent value84, same exponent profile: orbit counts1,31,331 agree; prime values are discarded by this count.")
    # Group and graph information deliberately do not agree.
    check(12 % 4 == 0 and g[12] % g[4] != 0)
    check(is_square(3+6) and not is_square(g[3]+g[6]))
    order = sorted(D)
    rank = {d: i+1 for i,d in enumerate(order)}
    label_action = {rank[d]: rank[g[d]] for d in D}
    before = [abs(j-1) for j in range(2,11)]
    after = [abs(label_action[j]-label_action[1]) for j in range(2,11)]
    check(sorted(before) == list(range(1,10)) and len(set(after)) < 9)
    report("Hostiles:4|12 but g(4)=12 does not divide g(12)=20; 3+6=9 but g(3)+g(6)=5+15=20; divisor group action is not an arithmetic graph automorphism.")
    report(f"On rank labels1..10, the gracefully labeled star centered1 has edge differences{before}; after transported C3 permutation, center={label_action[1]}, differences={after}, no longer distinct.")
    # The other nonprime balanced family generates binary repunit cardinalities.
    D8, S8, U8 = divisor_sets(8)
    check(D8 == {2,4} and S8 == U8 == {2})
    for k in range(1, 9):
        words = [w for w in product(sorted(D8), repeat=k) if w != (2,)*k]
        check(len(words) == 2**k-1)
    report("Cube profilep^3: two-letter proper-divisor words minus the sole squarefree constant word count2^k-1; mixed profile gives(10^k-7)/3 after the free C3 quotient.")
    sequence_profiles = []
    for k in range(2, 10):
        n = (10**k-7)//3
        fsu = tuple(len(x) for x in divisor_sets(n))
        sequence_profiles.append((k, n, fsu, fsu[0] == fsu[1]+fsu[2]))
    check(all(row[3] for row in sequence_profiles[:-1]) and not sequence_profiles[-1][3])
    report(f"Orbit COUNT values themselves: (k,value,(F,S,U),balanced)={sequence_profiles}")
    report("ALL CHECKS PASSED. Infinite orbit formulas and no-prefix196 conclusion have proofs in the note; no primality or Collatz rank follows from the word quotient.")
    output = "\n".join(lines)+"\n"
    destination = Path(__file__).resolve().parents[2] / "05-knowledge/results/duck_primes_20260925.out"
    destination.write_text(output, encoding="utf-8")
    print(output, end="")


if __name__ == "__main__":
    main()
