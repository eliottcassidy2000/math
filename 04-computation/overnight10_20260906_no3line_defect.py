#!/usr/bin/env python3
"""Exact controls for the Lipschitz diagonal-defect probability theorem.

Standard library only. Complete labelled boards at n=2..5, all row/column
transpositions, rook-polynomial means, small conditional reveal ranges,
and rational exponential enclosures. No producer or repository imports.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, permutations
from math import comb, factorial
import sys

sys.stdout.reconfigure(newline="\n")
GATES = 0


def need(ok, label):
    global GATES
    GATES += 1
    if not ok:
        raise RuntimeError(label)


def falling(n, k):
    if k > n:
        return 0
    return factorial(n) // factorial(n-k)


def partitions(n, minimum=2):
    if n == 0:
        yield ()
    for a in range(minimum, n+1):
        for tail in partitions(n-a, a):
            yield (a,) + tail


def poly_mul(a, b):
    c = [0] * (len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            c[i+j] += x*y
    return c


def rook_profile(parts):
    out = [1]
    for a in parts:
        q = 2*a
        p = [1]
        for k in range(1, a+1):
            p.append(comb(q-k, k) + comb(q-k-1, k-1))
        out = poly_mul(out, p)
    return out


def exact_mean(parts, limit=None):
    n = sum(parts)
    m = rook_profile(parts)
    end = n if limit is None else min(n, limit)
    return sum((Q((-1)**(k+1)*(k-2)*(2*n-k+1), k+1)
                * Q(m[k], falling(n, k))) for k in range(3, end+1))


def factorial_moment(parts, L, k):
    n = sum(parts)
    if k > min(n, L):
        return Q(0)
    m = rook_profile(parts)
    return Q(falling(L,k)*factorial(k)*m[k], falling(n,k)**2)


def board_generator(n):
    pairs = list(combinations(range(n), 2))
    counts = [0]*n
    rows = []
    def visit():
        if len(rows) == n:
            if counts == [2]*n:
                yield tuple(rows)
            return
        remain = n-len(rows)-1
        for a,b in pairs:
            if counts[a] == 2 or counts[b] == 2:
                continue
            counts[a] += 1
            counts[b] += 1
            if all(x+remain >= 2 for x in counts):
                rows.append((a,b))
                yield from visit()
                rows.pop()
            counts[a] -= 1
            counts[b] -= 1
    yield from visit()


def component_type(board):
    n = len(board)
    adjacency = [set() for _ in range(2*n)]
    for r, pair in enumerate(board):
        for c in pair:
            adjacency[r].add(n+c)
            adjacency[n+c].add(r)
    seen = set()
    sizes = []
    for root in range(2*n):
        if root in seen:
            continue
        seen.add(root)
        todo = [root]
        for x in todo:
            for y in adjacency[x]:
                if y not in seen:
                    seen.add(y)
                    todo.append(y)
        need(len(todo)%2 == 0, "bipartite component size")
        sizes.append(len(todo)//2)
    return tuple(sorted(sizes))


def occupancy(board):
    return Counter(c-r for r,pair in enumerate(board) for c in pair)


def defect(board):
    return sum(max(y-2,0) for y in occupancy(board).values())


def no_three(board):
    points = [(r,c) for r,pair in enumerate(board) for c in pair]
    return all((x[0]-z[0])*(y[1]-z[1]) != (y[0]-z[0])*(x[1]-z[1])
               for x,y,z in combinations(points,3))


def class_volume(parts):
    den = 1
    for a,m in Counter(parts).items():
        den *= (2*a)**m * factorial(m)
    value, rem = divmod(factorial(sum(parts))**2, den)
    need(rem == 0, "shore-preserving orbit volume")
    return value


def B4(n, c4):
    return (Q((2*n-5)*(n*n+5*n-11),15*(n-1)*(n-2))
            - Q(2*(2*n-3)*c4,5*n*(n-1)*(n-2)*(n-3)))


def label_board(board, rho, sigma):
    out = [None]*len(board)
    for r,pair in enumerate(board):
        out[rho[r]] = tuple(sorted(sigma[c] for c in pair))
    return tuple(out)


def conditional_range_control(board):
    n = len(board)
    perms = list(permutations(range(n)))
    sums = defaultdict(lambda: [0,0])
    for rho in perms:
        for sigma in perms:
            key = rho[:n-1] + sigma[:n-1]
            value = defect(label_board(board,rho,sigma))
            for depth in range(2*(n-1)+1):
                sums[key[:depth]][0] += value
                sums[key[:depth]][1] += 1
    children = defaultdict(list)
    for key,(total,count) in sums.items():
        if key:
            children[key[:-1]].append(Q(total,count))
    maximum = Q(0)
    for values in children.values():
        width = max(values)-min(values)
        maximum = max(maximum,width)
        need(width <= 4, "actual conditional Doob range")
    need(Q(*sums[()]) == exact_mean(component_type(board)),
         "permutation reveal normalization")
    return len(children),maximum


def main():
    print("DIAGONAL POSITIVE-PART DEFECT: EXACT CONTROLS")
    for y in range(81):
        f = max(y-2,0)
        for K in range(3,17):
            partial = sum((-1)**(k+1)*(k-2)*comb(y,k)
                          for k in range(3,min(y,K)+1))
            need((partial-f)*(-1)**(K+1) >= 0,
                 "pointwise alternating truncation direction")
    print("positive-part alternating truncations: y0..80,K3..16 PASS")

    board_total = 0
    swaps = 0
    global_swap = 0
    representatives = {}
    for n in range(2,6):
        data = defaultdict(lambda: [0,0,0,0])
        row_moments = defaultdict(lambda: Q(0))
        max_swap = 0
        for board in board_generator(n):
            parts = component_type(board)
            representatives.setdefault(parts,board)
            f = defect(board)
            good = no_three(board)
            need(not good or f == 0, "full success implies slope-one success")
            data[parts][0] += 1
            data[parts][1] += f
            data[parts][2] += int(f == 0)
            data[parts][3] += int(good)
            ys = occupancy(board)
            for d in range(0,n):
                L = n-d
                for k in range(0,L+1):
                    row_moments[(parts,L,k)] += falling(ys[d],k)
            for a,b in combinations(range(n),2):
                swapped = list(board)
                swapped[a],swapped[b] = swapped[b],swapped[a]
                dr = abs(defect(swapped)-f)
                col_swapped = tuple(tuple(sorted(b if x==a else a if x==b else x
                                                for x in pair)) for pair in board)
                dc = abs(defect(col_swapped)-f)
                need(dr <= 4 and dc <= 4, "literal row/column transposition range")
                max_swap = max(max_swap,dr,dc)
                swaps += 2
        for parts,(count,total,zero,safe) in sorted(data.items()):
            need(count == class_volume(parts), "complete labelled-board class universe")
            need(Q(total,count) == exact_mean(parts), "literal/rook exact defect mean")
            for L in range(1,n+1):
                for k in range(L+1):
                    need(row_moments[(parts,L,k)]/count == factorial_moment(parts,L,k),
                         "literal single-diagonal factorial moments")
            if n>=4:
                need(exact_mean(parts) >= B4(n,parts.count(2)) >= Q(2*n,15),
                     "finite fourth-order mean lower bound")
            print(f"n={n} type={parts}: boards={count} EF={Q(total,count)} "
                  f"F0={zero} X0={safe}")
            board_total += count
        print(f"  n={n} largest literal swap change={max_swap}")
        global_swap = max(global_swap,max_swap)
    need(global_swap == 4, "range four is attained")
    print(f"complete labelled universe: {board_total} boards; {swaps} swaps")

    sharp = tuple(tuple(sorted((i,(i+1)%5))) for i in range(5))
    changed = list(sharp)
    changed[0],changed[2] = changed[2],changed[0]
    need(defect(sharp) == 5 and defect(changed) == 1, "explicit sharp swap")
    print("sharp n5 identity+shift control: swapping rows0,2 changes F5 to F1")
    for parts in ((4,),(2,2)):
        count,maximum = conditional_range_control(representatives[parts])
        print(f"n4 type={parts}: {count} Doob nodes, max conditional range={maximum}")

    profile_count = 0
    for n in range(2,19):
        for parts in partitions(n):
            profile_count += 1
            m = rook_profile(parts)
            if n>=4:
                need(m[3] == Q(2*n*(n-2)*(2*n-5),3), "universal m3")
                need(m[4] == Q(n*(2*n-5)*(2*n-6)*(2*n-7),12)+parts.count(2),
                     "universal m4 with C4 correction")
                need(exact_mean(parts,4) == B4(n,parts.count(2)), "closed B4 formula")
                need(B4(n,parts.count(2)) >= Q(2*n,15), "all-profile finite lower bound")
            for L in (0,1,n//2,n):
                theta = Q(L,n)
                for k in range(2,n+4):
                    lam = 2*theta
                    actual = factorial_moment(parts,L,k)
                    delta = lam**k-actual
                    bound = 2**k*k*(k-1)*(theta**(k-1)/Q(2*n)
                                                 +theta**k/Q(4*(n-1)))
                    need(0 <= delta <= bound, "uniform all-k factorial deficit")
    print(f"all {profile_count} cycle profiles with total shore size2..18 PASS")

    for n in range(4,101):
        residual = (B4(n,Q(n,2))-Q(2*n,15))
        t = n-4
        need(residual == Q(11*t**3+48*t*t+58*t+12,
                            15*(n-1)*(n-2)*(n-3)), "positive shifted finite numerator")

    # A rational Taylor enclosure proves the two harmless numerical constants.
    K = 12
    e2_low = sum(Q(2**k,factorial(k)) for k in range(K+1))
    first_tail = Q(2**(K+1),factorial(K+1))
    e2_high = e2_low + first_tail/(1-Q(2,K+2))
    need(e2_low > 7 and e2_high < Q(37,5), "rational exp(2) enclosure")
    lower_C4_upper = (13*e2_high-63/e2_high+12)/6
    upper_C4_plus_g2 = (13*e2_high+87/e2_high-12)/6
    need(lower_C4_upper < 17, "uniform lower mean error constant17")
    need(upper_C4_plus_g2 < 16, "uniform upper mean error constant16")
    need(Q(2*2,15)**2 / Q(16) == Q(1,225), "range normalization control")
    need(Q(4,225)/16 == Q(1,900), "finite concentration exponent constant")
    print("rational Taylor bounds: n*alpha-17 <= EF <= n*alpha+16 PASS")
    print("concentration denominator16(n-1); all-n exponent n^2/[900(n-1)]")
    print(f"PASS {GATES} always-active gates")


if __name__ == "__main__":
    main()
