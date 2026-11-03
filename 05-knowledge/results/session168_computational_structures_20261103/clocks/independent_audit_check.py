#!/usr/bin/env python3
"""Independent controls for signed-XOR exponent and literal guard clocks."""
import collections
import itertools
import json
import math


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def guard_clock(x):
    result = []
    for i,(s,u,d) in enumerate(x):
        ns,nu,nd = x[(i+1)%len(x)]
        result.append(((u+ns)%3,(s+u+d+nu)%3,(s+nd)%3))
    return tuple(result)


def main():
    signed_xor_counts = {}
    for level in range(4):
        size = 1 << level
        count = 0
        for shift in range(size):
            for function in itertools.product(range(2), repeat=size):
                p = tuple((s ^ function[x])*size + (x ^ shift) for s in range(2) for x in range(size))
                need(len(set(p)) == 2*size, 'signed-XOR permutation')
                need(all(p[p[p[p[u]]]] == u for u in range(2*size)), 'universal exponent four')
                unseen = set(range(2*size))
                cycles = collections.Counter()
                while unseen:
                    root = min(unseen)
                    u = root
                    length = 0
                    while True:
                        unseen.remove(u)
                        length += 1
                        u = p[u]
                        if u == root:
                            break
                    cycles[length] += 1
                if shift == 0:
                    expected = {1:2*function.count(0),2:function.count(1)}
                else:
                    differences = [function[x]^function[x^shift] for x in range(size)]
                    expected = {2:differences.count(0),4:differences.count(1)//2}
                expected = {length:count for length,count in expected.items() if count}
                need(dict(cycles) == expected, 'cycle multiplicities')
                count += 1
        signed_xor_counts[level] = count
    clocks = []
    nonzero = set(itertools.product(range(3),repeat=3)) - {(0,0,0)}
    for width in (1,3,9,27,81):
        initial = ((0,0,1),) + ((0,0,0),)*(width-1)
        x = initial
        occupancy = collections.Counter()
        for time in range(26*width):
            occupancy[x[0]] += 1
            x = guard_clock(x)
            if time+1 < 26*width:
                need(x != initial, 'guard delta has no earlier return')
        need(x == initial, 'guard delta return')
        need(occupancy == collections.Counter({v:width for v in nonzero}), 'uniform local occupation')
        x, unit_phase = initial, 0
        for time in range(104*width):
            if x[0] == (0,0,1):
                unit_phase = (unit_phase+1)%4
            x = guard_clock(x)
            if time+1 < 104*width:
                need((x,unit_phase) != (initial,0), 'four-lap first return')
        need((x,unit_phase) == (initial,0), 'four-lap return')
        clocks.append({'width':width,'delta_period':26*width,'four_lap_period':104*width,'triggers_per_lap':width})
    valuation_controls = []
    for width in (3,9,27):
        for valuation in sorted({0,1,width-1,width-2,width-3,width//2}):
            # (z-1)^v has valuation v; reversing spatial orientation preserves it.
            initial = tuple((0,0,math.comb(valuation,i)*(-1)**(valuation-i)%3) if i<=valuation else (0,0,0) for i in range(width))
            q = 1
            while q < width-valuation:
                q *= 3
            expected = 26*q
            x = initial
            for time in range(1,expected+1):
                x = guard_clock(x)
                if x == initial:
                    break
            need(x == initial and time == expected, 'independent valuation first return')
            valuation_controls.append({'width':width,'valuation':valuation,'period':time})
    print(json.dumps({'all_signed_xor_maps':signed_xor_counts,'guard_clock_controls':clocks,
                      'valuation_controls':valuation_controls},sort_keys=True,indent=2))


if __name__ == '__main__':
    main()
