#!/usr/bin/env python3
"""One complete regular11 -> strong hub -> transitive certificate.

This is an illustrative route, not the full 1223-class coverage check.
The upper-triangle bit order is (0,1),(0,2),...,(9,10).
"""
from collections import Counter

BITS = '1111100000111100010111110001100111111101111110011100101'
PATHS = [[9,0,5,6,3,10,8,2,7,1,4],[4,9,0,10,6,3,1,8,2,7,5],[5,9,10,3,8,6,1,7,0,2,4],[3,8,6,2,1,7,4,0,9,10,5],[10,2,7,4,6,9,5,8,0,3,1],[1,10,0,5,4,3,2,6,8,7,9],[7,6,4,3,0,8,10,5,9,2,1],[3,7,6,1,0,10,2,8,4,5,9],[8,4,9,6,1,0,3,5,10,7,2],[3,5,7,4,9,6,2,1,0,8,10],[9,6,3,0,8,10,2,5,7,4,1],[3,0,4,1,2,6,7,8,10,9,5],[8,0,4,1,3,7,10,9,5,2,6],[2,4,1,3,7,0,6,10,9,5,8],[9,5,8,10,7,0,6,2,4,1,3]]


def require(ok, message):
    if not ok:
        raise ValueError(message)


def decode(bits):
    require(len(bits) == 55 and set(bits) <= {'0','1'}, 'bad bit word')
    a = [[False]*11 for _ in range(11)]
    t = 0
    for i in range(11):
        for j in range(i+1,11):
            a[i][j] = bits[t] == '1'; a[j][i] = not a[i][j]; t += 1
    return a


def scores(a):
    return list(map(sum,a))


def tariff(a):
    return sum(abs(v-i) for i,v in enumerate(sorted(scores(a))))//2


def move(a,p):
    require(sorted(p) == list(range(11)), 'path repeats or omits a vertex')
    require(all(a[x][y] for x,y in zip(p,p[1:])), 'path arc is backwards')
    keep = {frozenset((x,y)) for x,y in zip(p,p[1:])}
    return [[False if i == j else (a[i][j] if frozenset((i,j)) in keep else not a[i][j]) for j in range(11)] for i in range(11)]


def blue(a,p):
    reverse = {p[i]:p[10-i] for i in range(11)}
    return all(a[i][j] == a[reverse[j]][reverse[i]] for i in range(11) for j in range(11))


def strong(a):
    for root in range(11):
        seen = {root}; todo = [root]
        while todo:
            i = todo.pop()
            for j in range(11):
                if a[i][j] and j not in seen:
                    seen.add(j); todo.append(j)
        if len(seen) != 11:
            return False
    return True


def main():
    a = decode(BITS)
    require(scores(a) == [5]*11 and tariff(a) == 15, 'not regular eleven')
    endpoints = Counter(); increments = [[] for _ in range(11)]
    for t,p in enumerate(PATHS,1):
        require(strong(a), 'route left strong locus before final move')
        require(blue(a,p) == (t == 15), 'incorrect path colour')
        old = scores(a)
        b = move(a,p)
        require(move(b,p) == a, 'native involution failed')
        new = scores(b)
        require(new == [10-old[v]+int(v == p[0])-int(v == p[-1]) for v in range(11)], 'score update')
        require(tariff(b) == 15-t, 'non-tight L change')
        for v in range(11):
            change = (-1)**t*(new[v]-5)-(-1)**(t-1)*(old[v]-5)
            require(change == (-1)**t*(int(v == p[0])-int(v == p[-1])), 'clock identity')
            if change:
                increments[v].append(change)
        endpoints.update((p[0],p[-1])); a = b
    require(sorted(scores(a)) == list(range(11)), 'target not transitive')
    for v,j in enumerate(scores(a)):
        require(endpoints[v] == abs(j-5), 'endpoint multiplicity')
        require(not increments[v] or len(set(increments[v])) == 1, 'coordinate cancellation')
    require(sorted(endpoints[v] for v in range(11)) == [0,1,1,2,2,3,3,4,4,5,5], 'wrong inventory')
    try:
        move(decode(BITS),[0]*11)
    except ValueError:
        pass
    else:
        raise ValueError('malformed path not rejected')
    print('PASS: one complete 14-black + 1-blue route, strong corridor, involutions and endpoint clocks; malformed path rejected')


if __name__ == '__main__':
    main()
