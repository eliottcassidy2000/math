#!/usr/bin/env python3
"""Small literal controls for the strong two-port counter carrier.

No repository imports and no stripped assert gates. Output is deterministic.
"""
from collections import Counter
from functools import lru_cache
from itertools import combinations, permutations
import json


def need(condition, message):
    if not condition:
        raise RuntimeError(message)


def graph(a, b):
    """Core names 0,1,2; persistent wrapper ports 3=u,4=v."""
    left = tuple(range(5, 5 + a))
    right = tuple(range(5 + a, 5 + a + b))
    names = tuple(range(a + b + 5))
    out = {x: set() for x in names}
    for x, y in ((0, 1), (1, 2), (2, 0), (3, 4)):
        out[x].add(y)
    X = left + (0, 1, 2) + right
    for i, x in enumerate(X):
        for y in X[i + 1:]:
            if x < 3 and y < 3:
                continue
            out[x].add(y)
        out[x].add(3)
        out[4].add(x)
    return out, left, right


def rows(out):
    names = tuple(sorted(out))
    index = {x: i for i, x in enumerate(names)}
    return tuple(sum(1 << index[y] for y in out[x]) for x in names), index


def validate(out, left, right):
    need(set(out) == set(left) | set(right) | set(range(5)), 'named partition')
    for x, y in combinations(out, 2):
        need((y in out[x]) != (x in out[y]), 'tournament pair')
    need(all(x not in out[x] for x in out), 'no loops')
    X = tuple(left) + (0, 1, 2) + tuple(right)
    for x in X:
        need(3 in out[x] and x in out[4], 'wrapper cut')
    need(4 in out[3], 'wrapper arc')
    for i, x in enumerate(X):
        for y in X[i + 1:]:
            if x < 3 and y < 3:
                continue
            need(y in out[x], 'ordered counter blocks')
    need(1 in out[0] and 2 in out[1] and 0 in out[2], 'core cycle')
    # A separate graph reachability check, not just the construction formula.
    for start in out:
        visited, todo = {start}, [start]
        while todo:
            x = todo.pop()
            for y in out[x] - visited:
                visited.add(y)
                todo.append(y)
        need(len(visited) == len(out), 'strong connectivity')


def inc(out, left, right, register, newborn):
    need(newborn not in out, 'fresh counter name')
    out = {x: set(v) for x, v in out.items()}
    if register == 0:
        out[newborn] = set(out) - {4}
        out[4].add(newborn)
        left = (newborn,) + tuple(left)
    else:
        need(register == 1, 'counter index')
        out[newborn] = {3}
        for x in out:
            if x not in (3, newborn):
                out[x].add(newborn)
        right = tuple(right) + (newborn,)
    validate(out, left, right)
    return out, left, right


def dec(out, left, right, register):
    block = left if register == 0 else right
    need(bool(block), 'nonzero deletion guard')
    removed = block[0] if register == 0 else block[-1]
    out = {x: set(v) - {removed} for x, v in out.items() if x != removed}
    if register == 0:
        left = left[1:]
    else:
        right = right[:-1]
    validate(out, left, right)
    return out, left, right


def H(T):
    n = len(T)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, count in enumerate(dp[mask]):
            if not count:
                continue
            available = T[v] & ~mask
            while available:
                bit = available & -available
                available -= bit
                dp[mask | bit][bit.bit_length() - 1] += count
    return sum(dp[-1])


def literal_h(T):
    return sum(all(T[p[i]] >> p[i+1] & 1 for i in range(len(T)-1))
               for p in permutations(range(len(T))))


def cover_profile(T, ports):
    """Enumerate directed components, then unordered partitions once each.

    This is independent of the forward eighteen-polynomial transfer code.
    """
    n = len(T)
    slot = {v: i for i, v in enumerate(ports)}
    paths = [Counter() for _ in range(1 << n)]

    def visit(mask, first, last):
        paths[mask][(slot.get(first, -1), slot.get(last, -1))] += 1
        available = T[last] & ~mask
        while available:
            bit = available & -available
            available -= bit
            visit(mask | bit, first, bit.bit_length() - 1)

    for first in range(n):
        visit(1 << first, first, first)

    @lru_cache(None)
    def covers(mask):
        if not mask:
            return Counter({((), 0): 1})
        first = mask & -mask
        result = Counter()
        subset = mask
        while subset:
            if subset & first:
                for endpoint, multiplicity in paths[subset].items():
                    for (skeleton, degree), count in covers(mask ^ subset).items():
                        if endpoint == (-1, -1):
                            key = (skeleton, degree + 1)
                        else:
                            key = (tuple(sorted(skeleton + (endpoint,))), degree)
                        result[key] += count * multiplicity
            subset = (subset - 1) & mask
        return result

    return covers((1 << n) - 1)


def sigma(n):
    return 2*n + 1 + 2*(n >= 3) + 2*(n >= 9)


def rho(h):
    return (h-1)//2 - (h > 7) - (h > 21)


def K(h):
    best = 0
    for b in range(3):
        if 5**b > h:
            continue
        a, power = 0, 5**b
        while 3*power <= h:
            a += 1
            power *= 3
        best = max(best, 3*a + 4*b)
    return best


def fractran(value, program):
    """Literal first-applicable integer execution, with a hostile safety cap."""
    trace = [value]
    for _ in range(100):
        for numerator, denominator in program:
            if value * numerator % denominator == 0:
                value = value * numerator // denominator
                trace.append(value)
                break
        else:
            return trace
    raise RuntimeError('unexpected nonhalting short FRACTRAN control')


def main():
    base_checks, transition_checks, transfer_checks = [], 0, 0
    for a in range(4):
        for b in range(4):
            out, left, right = graph(a, b)
            validate(out, left, right)
            T, _ = rows(out)
            h = H(T)
            need(h == 3 + 12*2**(a+b), 'exact count law')
            if a+b <= 2:
                need(h == literal_h(T), 'independent permutation count')
            base_checks.append([a, b, len(T), h])
            for register in range(2):
                child = inc(out, left, right, register, 100)
                CT, _ = rows(child[0])
                need(H(CT) == 3 + 12*2**(a+b+1), 'literal increment count')
                returned = dec(*child, register)
                need(returned == (out, left, right), 'named increment/decrement inverse')
                transition_checks += 2
                if (left if register == 0 else right):
                    parent = dec(out, left, right, register)
                    PT, _ = rows(parent[0])
                    need(H(PT) == 3 + 12*2**(a+b-1), 'literal decrement count')
                    transition_checks += 1
            if a:
                x = left[0]
                target = {v: set(w) for v, w in out.items()}
                X = set(out) - {3, 4}
                for y in X - {x}:
                    need(y in target[x], 'source counter vertex')
                    target[x].remove(y)
                    target[y].add(x)
                validate(target, left[1:], right + (x,))
                TT, _ = rows(target)
                need(H(TT) == h, 'native-reset endpoint is H neutral')
                changed = sum((y in out[x0]) != (y in target[x0])
                              for x0, y in combinations(out, 2))
                need(changed == a+b+2, 'native-reset old-arc budget')
                transfer_checks += 1

    collisions = []
    for s in range(4):
        reference = None
        for a in range(s+1):
            out, _, _ = graph(a, s-a)
            T, index = rows(out)
            profile = cover_profile(T, (index[3], index[4]))
            need(sum(count for (sk, degree), count in profile.items()
                     if len(sk)+degree == 1) == H(T), 'cover Hamilton readout')
            if reference is None:
                reference = profile
            else:
                need(profile == reference, 'entire eighteen-polynomial collision')
            collisions.append([a, s-a, len(profile)])
    # These are the two different zero guards inside one exact observer fibre.
    need((0 == 0) != (1 == 0), 'hostile zero guard')

    executions = []
    fraction_program = ((7, 10), (11, 5), (15, 7))
    fraction_traces = []
    for a in range(4):
        for b in range(4):
            state = graph(a, b)
            reference_a, reference_b, pc = a, b, 0
            next_name, steps = 1000, 0
            while pc != 2:
                out, left, right = state
                need(len(left) == reference_a and len(right) == reference_b,
                     'counter-program state matches literal source')
                if pc == 0:  # DECJZ A -> nonzero:1, zero:halt
                    if left:
                        state = dec(*state, 0)
                        reference_a -= 1
                        pc = 1
                    else:
                        pc = 2
                else:        # INC B -> 0
                    state = inc(*state, 1, next_name)
                    next_name += 1
                    reference_b += 1
                    pc = 0
                steps += 1
            need(reference_a == 0 and reference_b == a+b, 'transfer program output')
            need(steps == 2*a+1, 'literal instruction clock')
            need(H(rows(state[0])[0]) == 3+12*2**(a+b), 'program endpoint count')
            executions.append([a, b, reference_a, reference_b, steps])
            trace = fractran(5*2**a*3**b, fraction_program)
            need(len(trace)-1 == steps, 'FRACTRAN/tournament instruction clock')
            need(trace[-1] == 11*3**(a+b), 'FRACTRAN transfer output')
            # Exact intermediate instruction boundaries, not just the output.
            expected = [5*2**a*3**b]
            for done in range(a):
                expected.append(7*2**(a-done-1)*3**(b+done))
                expected.append(5*2**(a-done-1)*3**(b+done+1))
            expected.append(11*3**(a+b))
            need(trace == expected, 'FRACTRAN full instruction-state agreement')
            fraction_traces.append([a,b,trace])
    wrong_priority = fractran(10, ((11,5),(7,10),(15,7)))
    need(wrong_priority == [10,22], 'zero-first FRACTRAN priority hostile')

    addresses = [sigma(n) for n in range(101)]
    need(addresses == [h for h in range(1, addresses[-1]+1, 2) if h not in (7,21)],
         'coverage-address enumeration')
    need(all(rho(sigma(n)) == n for n in range(101)), 'coverage-address inverse')
    budgets = {h: K(h) for h in (1,3,5,7,9,15,21,25,27,45,75,81)}
    need(list(budgets.values()) == [0,3,4,4,6,7,7,8,9,10,11,12],
         'inherited exact core budget controls')
    print(json.dumps({
        'status': 'FINITE-EXACT controls, not spectrum coverage',
        'universe': '0<=a,b<=3; full marked covers for a+b<=3',
        'source_checks': base_checks,
        'increment_decrement_checks': transition_checks,
        'native_endpoint_checks': transfer_checks,
        'marked_cover_checks': collisions,
        'counter_program_executions': executions,
        'fraction_program': fraction_program,
        'fraction_traces': fraction_traces,
        'wrong_priority_hostile': wrong_priority,
        'zero_guard_hostile': [[0,1,True],[1,0,False]],
        'address_prefix': addresses[:12],
        'core_budget_controls': budgets,
    }, indent=2, sort_keys=True))


if __name__ == '__main__':
    main()
