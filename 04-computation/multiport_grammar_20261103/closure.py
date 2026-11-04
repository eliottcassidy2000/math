#!/usr/bin/env python3
"""Exact controls for the ordered-end, old-port-preserving obstruction.

The unbounded theorem is proved in closure/PROOF.md. This finite control uses
all labeled tournaments through order 5, all actual port masks, and all labeled
new blocks through order 3. It never quotients by an observed marginal.
No assertions: all verification gates remain active under python -O.
"""
from collections import Counter
from functools import lru_cache
from itertools import permutations


GATES = 0


def check(test, label):
    global GATES
    GATES += 1
    if not test:
        raise RuntimeError(label)


def tournament(n, code):
    rows = [0] * n
    bit = 0
    for a in range(n):
        for b in range(a + 1, n):
            if (code >> bit) & 1:
                rows[a] |= 1 << b
            else:
                rows[b] |= 1 << a
            bit += 1
    return tuple(rows)


def transitive(n):
    return tuple(sum(1 << j for j in range(i + 1, n)) for i in range(n))


def join(left, right):
    n = len(left)
    right_mask = ((1 << len(right)) - 1) << n
    return tuple(row | right_mask for row in left) + tuple(row << n for row in right)


@lru_cache(None)
def incidence(rows):
    """Literal ordered permutations, independently of the transfer formula."""
    n = len(rows)
    full = (1 << n) - 1
    paths = [[] for _ in range(1 << n)]
    paths[0] = [()]
    for size in range(1, n + 1):
        for path in permutations(range(n), size):
            if all((rows[a] >> b) & 1 for a, b in zip(path, path[1:])):
                paths[sum(1 << v for v in path)].append(path)
    h = len(paths[full])
    p = 0
    first_end = [0] * n
    second_start = [0] * n
    joint = [[0] * n for _ in range(n)]
    ends = [0] * n
    starts = [0] * n
    for path in paths[full]:
        ends[path[-1]] += 1
        starts[path[0]] += 1
    for mask, left_paths in enumerate(paths):
        right_paths = paths[full ^ mask]
        for left in left_paths:
            for right in right_paths:
                p += 1
                if left:
                    first_end[left[-1]] += 1
                if right:
                    second_start[right[0]] += 1
                if left and right:
                    joint[left[-1]][right[0]] += 1
    return h, p, first_end, second_start, joint, ends, starts


def signature(rows, mask):
    h, p, fv, gv, jv, ev, sv = incidence(rows)
    ports = [i for i in range(len(rows)) if (mask >> i) & 1]
    return (h, p, sum(fv[i] for i in ports), sum(gv[i] for i in ports),
            sum(jv[i][j] for i in ports for j in ports),
            sum(ev[i] for i in ports), sum(sv[i] for i in ports))


def compose(x, y):
    h, p, f, g, j, e, s = x
    H, P, F, G, J, E, S = y
    return (h*H, p*P, f*H+p*F, g*P+h*G,
            j*H+h*J+e*S+g*F, h*E, s*H)


@lru_cache(None)
def reach(rows):
    out = [row | (1 << i) for i, row in enumerate(rows)]
    for k in range(len(rows)):
        for i in range(len(rows)):
            if (out[i] >> k) & 1:
                out[i] |= out[k]
    return tuple(out)


def terminal_port(rows, mask):
    r = reach(rows)
    return any((mask >> i) & 1 and all((v >> i) & 1 for v in r)
               for i in range(len(rows)))


def attachment(rows, mask, clock):
    """Names are prefix 0..clock-1, suffix clock.., then the physical root."""
    out = list(join(transitive(clock), rows))
    root = len(out)
    root_row = (1 << clock) - 1
    for v in range(len(rows)):
        if (mask >> v) & 1:
            out[clock + v] |= 1 << root
        else:
            root_row |= 1 << (clock + v)
    out.append(root_row)
    return tuple(out)


def strong(rows):
    return all(v == (1 << len(rows)) - 1 for v in reach(rows))


@lru_cache(None)
def h_dp(rows):
    """Independent subset/end-point DP, not the ordered-cover enumerator."""
    n = len(rows)
    dp = [[0] * n for _ in range(1 << n)]
    for v in range(n):
        dp[1 << v][v] = 1
    for mask in range(1, 1 << n):
        for v, value in enumerate(dp[mask]):
            nxt = rows[v] & ~mask
            while value and nxt:
                bit = nxt & -nxt
                nxt -= bit
                dp[mask | bit][bit.bit_length() - 1] += value
    return sum(dp[-1])


def full_state_control(rows, mask, sig):
    h, p, f, g, j, e, s = sig
    for clock in (1, 2):
        actual = attachment(rows, mask, clock)
        check(h_dp(actual) == h-j+(1 << clock)*f, 'actual attachment count')
        check(strong(actual) == terminal_port(rows, mask), 'terminal SCC guard')


def main():
    all_states = []
    small = []
    guarded = []
    tournaments = Counter()
    for n in range(1, 6):
        for code in range(1 << (n*(n-1)//2)):
            rows = tournament(n, code)
            h, p, fv, gv, jv, ev, sv = incidence(rows)
            tournaments[n] += 1
            check(h == h_dp(rows) and h > 0 and h % 2 == 1, 'Hamilton parity')
            check(p % 2 == 0 and p >= 2*h, 'two-slot parity and empty slots')
            check((p == 2) == (n == 1), 'two-slot singleton rigidity')
            for mask in range(1 << n):
                sig = signature(rows, mask)
                h, p, f, g, j, e, s = sig
                check(0 <= f <= p-h, 'first-slot-empty exclusion')
                check(f >= mask.bit_count()*h, 'marked cut injection')
                state = rows, mask, sig
                all_states.append(state)
                if n <= 3:
                    small.append(state)
                if terminal_port(rows, mask):
                    guarded.append(state)
                if n <= 4:
                    full_state_control(*state)

    transfers = 0
    for left, lm, ls in small:
        for right, rm, rs in small:
            if len(left) + len(right) > 5:
                continue
            actual = join(left, right)
            mask = lm | (rm << len(left))
            check(signature(actual, mask) == compose(ls, rs), 'seven-coordinate law')
            transfers += 1

    # Entire declared raw extension universe; actual masks and rows are retained.
    # The equation f'=3f is the slope condition, not an extrapolated finite clock.
    first = []
    modes = Counter()
    digits = {'append': Counter(), 'prepend': Counter()}
    tests = 0
    for rows, mask, sig in guarded:
        h, p, f, g, j, e, s = sig
        for block, bm, bs in small:
            tests += 2
            app = compose(sig, bs)
            if terminal_port(block, bm) and app[2] == 3*f:
                check(len(block) == 1 and bm == 1 and p == 2*f,
                      'append singleton classification')
                check(app == (h,2*p,3*f,2*g+h,j+e+g,h,s), 'append full state')
                delta = app[0]-app[4]-3*(h-j)
                check(delta == 2*j-2*h-e-g, 'append intercept')
                first.append(('append', rows, mask, sig, app))
                modes['append'] += 1
                digits['append'][delta] += 1
            pre = compose(bs, sig)
            if pre[2] == 3*f:
                check(len(block) == 1 and bm == 1 and f == h,
                      'prepend singleton classification')
                check(pre == (h,2*p,3*f,p+g,j+s+f,e,h), 'prepend full state')
                check(mask.bit_count() == 1 and j == 0, 'prepend zero joint defect')
                delta = pre[0]-pre[4]-3*(h-j)
                check(delta == -3*h-s < 0, 'prepend negative intercept')
                first.append(('prepend', rows, mask, sig, pre))
                modes['prepend'] += 1
                digits['prepend'][delta] += 1

    continuations = 0
    for mode, rows, mask, old, mid in first:
        if mode == 'append':
            mr = join(rows, (0,)); mm = mask | (1 << len(rows))
            check(tuple(v & ((1 << len(rows))-1) for v in mr[:-1]) == rows
                  and mm & ((1 << len(rows))-1) == mask, 'append old native state')
        else:
            mr = join((0,), rows); mm = 1 | (mask << 1)
            check(tuple(v >> 1 for v in mr[1:]) == rows and mm >> 1 == mask,
                  'prepend old native state under the displayed name embedding')
        check(signature(mr, mm) == mid, 'actual first-stage mask and adjacency')
        full_state_control(mr, mm, mid)
        for block, bm, bs in small:
            if terminal_port(block, bm):
                check(compose(mid, bs)[2] != 3*mid[2], 'no next append')
                continuations += 1
            check(compose(bs, mid)[2] != 3*mid[2], 'no next prepend')
            continuations += 1

    # Sharpness: one genuine step exists, but its digit is -4.
    one = (0,)
    one_sig = signature(one, 1)
    two_sig = signature(transitive(2), 3)
    check(one_sig == (1,2,1,1,0,1,1), 'literal singleton')
    check(two_sig == (1,4,3,3,2,1,1), 'literal first continuation')
    for clock in range(1, 6):
        check(h_dp(attachment(transitive(2), 3, clock)) ==
              3*h_dp(attachment(one, 1, clock))-4, 'nonvacuous affine tripling')

    # Inherited exact +2/+4 hostiles to an overbroad no-positive-growth claim.
    base = transitive(3); base_mask = 7
    plus4 = transitive(5); plus4_mask = 21
    plus2 = join(transitive(2), (2,4,1)); plus2_mask = 7
    check(signature(base, base_mask) == (1,8,7,7,6,1,1), 'inherited entry')
    check(signature(plus4, plus4_mask) == (1,32,21,21,12,1,1), '+4 actual ports')
    check(signature(plus2, plus2_mask) == (3,48,21,39,16,1,3), '+2 actual ports')
    check((plus4_mask & 7) != base_mask, '+4 exits by changing an old port')
    check((plus2_mask & 7) == base_mask and bool(plus2[4] & (1 << 2)),
          '+2 preserves ports but has a new-to-old suffix arc')
    for out, om, digit in ((plus2, plus2_mask, 2), (plus4, plus4_mask, 4)):
        for clock in (1,2,3):
            check(strong(attachment(out, om, clock)), 'positive hostile strongness')
            check(h_dp(attachment(out, om, clock)) ==
                  3*h_dp(attachment(base, base_mask, clock))+digit,
                  'positive hostile actual count')

    # Formal parity hostile: dropping realizability admits prepend then append.
    fake = (1,3,1,1,0,1,1)
    middle = compose(one_sig, fake)
    end = compose(middle, one_sig)
    check(fake[1] % 2 == 1 and middle[2] == 3*fake[2] and end[2] == 3*middle[2],
          'parity is a load-bearing realization sidecar')

    print('Ordered-end old-port-preserving closure: exact finite controls PASS')
    print('Labeled tournaments by order:', dict(sorted(tournaments.items())))
    print('Actual suffix/port states:', len(all_states), '; strong-attachment states:', len(guarded))
    print('Direct seven-coordinate transfer checks:', transfers)
    print('Raw append/prepend slope tests:', tests)
    print('One-step actual states:', dict(sorted(modes.items())))
    for mode in ('append','prepend'):
        print(mode, 'finite offset distribution:', dict(sorted(digits[mode].items())))
    print('Second-step exclusions in declared finite block bank:', continuations)
    print('Positive +2/+4 controls: actual and strong; outside ordered-end/old-port scope')
    print('Formal odd-p two-step control: rejected by actual two-slot parity')
    print('Unbounded consequence is analytic; finite negative-offset census is not generalized')
    print('Always-active gates:', GATES)


if __name__ == '__main__':
    main()
