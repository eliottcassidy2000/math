"""Three equal unit atoms, exact nonadjacent fibres, and fixed-source role loss."""
from itertools import product
from math import isqrt

AUTHORITATIVE = 'RKBRKRKBRKBRKRKBRKRKBRKBRKRKBRKBRKB'


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


G = [1, 1, 1]
for _ in range(24):
    G.append(G[-1] + G[-2])


def zeck(n):
    support = []
    for i in range(len(G)-1, 1, -1):
        if G[i] <= n:
            support.append(i)
            n -= G[i]
    need(n == 0, 'Zeckendorf table too small')
    return tuple(reversed(support))


def Q(n):
    m = n + 1
    b = (3*m-isqrt(5*m*m)-1)//2
    return n-2*b, b


def charge(indices):
    A = B = 0
    for i in indices:
        v = (-1,1) if i == 0 else (1,0) if i == 1 else Q(G[i])
        A += v[0]
        B += v[1]
    return A, B


def predicted_fibre(n):
    ordinary = zeck(n)
    lower = zeck(n-1)
    result = {ordinary, (0,) + lower}
    if 2 not in lower:
        result.add((1,) + lower)
    return result


def all_nonadjacent(max_index):
    result = {}
    def visit(i, blocked, n, indices):
        if i > max_index:
            result.setdefault(n, set()).add(indices)
            return
        visit(i+1, False, n, indices)
        if not blocked:
            visit(i+1, True, n+G[i], indices+(i,))
    visit(0, False, 0, ())
    return result


def preferred_words(top):
    words = ['K','B','R']
    for i in range(3, top+1):
        indices = range(i-1, -1, -2)
        words.append(''.join(words[j] for j in indices))
    return words


def row_word(n, words):
    return ''.join(words[i] for i in reversed(zeck(n)))


def v2(n):
    return (n & -n).bit_length()-1


def main():
    need(len(AUTHORITATIVE) == 35, 'data length')
    words = preferred_words(26)
    need(words[26][:34] == AUTHORITATIVE[:34], 'first34')
    need(words[26][34] == 'R' and AUTHORITATIVE[34] == 'B', '35 hostile retained')
    print('authoritative word:', AUTHORITATIVE)
    print('three-unit preferred grammar matches first34; predicted R versus supplied B at35')

    census = all_nonadjacent(20)
    counts = {2:0, 3:0}
    marked = {'B':0, 'K':0}
    for n in range(1, G[20]):
        actual = census[n]
        expected = predicted_fibre(n)
        need(actual == expected, ('exact fibre', n))
        counts[len(actual)] += 1
        good = {rep for rep in actual if charge(rep) == Q(n)}
        need(len(good) == 2, ('two charge sheets', n, good))
        lower = zeck(n-1)
        kind = 'K' if 2 in lower else 'B'
        expected_marker = (0,) if kind == 'K' else (1,)
        need(good == {zeck(n), expected_marker+lower}, ('marked sheet', n))
        marked[kind] += 1
        A,B = Q(n)
        need(A+2*B == n, ('visible value', n))
        prev = Q(n-1)
        need((A-prev[0],B-prev[1]) == ((-1,1) if 2 in lower else (1,0)), ('increment', n))
        need(row_word(n, words) == words[26][:n], ('prefix coherence', n))
    print('all nonadjacent fibres for n1..6764:', counts)
    print('exact two charge-preserving sheets at every n; marked sheet kinds:', marked)
    print('ordinary Zeckendorf expansion gives nested preferred row prefixes through6764')

    # All proper Fibonacci expansion choices, independently composed as languages.
    languages = [{'K'},{'B'},{'R'}]
    for i in range(3, 10):
        reps = predicted_fibre(G[i]) - {(i,)}
        current = set()
        for rep in reps:
            for pieces in product(*(languages[j] for j in reversed(rep))):
                current.add(''.join(pieces))
        expected = {''}
        for symbol in words[i]:
            expected = {prefix+c for prefix in expected for c in ('BK' if symbol == 'B' else symbol)}
        need(current == expected, ('complete word language', i))
        languages.append(current)
    print('complete Fibonacci expansion languages throughF9=34: arbitrary B-to-K replacements only')
    print('language counts for atoms indices3..9:', [len(s) for s in languages[3:]])

    def row_language(n):
        result = set()
        for rep in predicted_fibre(n):
            for pieces in product(*(languages[j] for j in reversed(rep))):
                result.add(''.join(pieces))
        return result

    rows35 = row_language(35)
    rows36 = row_language(36)
    need(AUTHORITATIVE in rows35, 'isolated35 legal')
    need(not any(word.startswith(AUTHORITATIVE) for word in rows36), 'no nested36')
    need(all(word[34] == 'R' for word in rows36), 'forced red35 in row36')
    need(charge((1,9)) == Q(35) == charge((2,9)), 'red blue full-charge collision')
    print('isolated row35 authoritative word is legal and charge-correct; every row36 has R at35')
    print('row36 supports:', sorted(predicted_fibre(36)))

    first = None
    branches = {}
    triple_total = 0
    restricted_total = 0
    for y in range(3, 4096, 2):
        t = v2(3*y+1)
        seen = set()
        count = 0
        restricted = 0
        for x in range(2, y, 2):
            r = v2(x)
            q = x >> r
            u = y-x
            a = v2(3*u+1)
            kind = 'consume' if a < r else 'swap' if a > r else 'collision'
            expected = 'consume' if r > t else 'swap' if r == t else 'collision'
            need(kind == expected, ('branch from fixed source', y,q,r,u))
            seen.add(kind)
            count += 1
            restricted += (q % 3 == 0) != (u % 3 == 0)
        need(count == (y-1)//2, ('triple fibre count', y))
        triple_total += count
        need(restricted == (y//3 if y % 3 else 0), ('one-divisible3 fibre', y))
        restricted_total += restricted
        if len(seen)==3 and len(predicted_fibre(y))==2 and first is None:
            first = y
        if y == 41:
            for q,r in ((1,1),(3,2),(1,3)):
                u=y-(q << r)
                a=v2(3*u+1)
                need((q % 3 == 0) != (u % 3 == 0), ('role-compatible witness',q,r,u))
                branches[r]=(q,r,u,a)
            need(restricted == 13, '41 restricted states')
    need(first == 41, ('minimal two-row three-branch witness', first))
    print('fixed-source decomposition controls:', triple_total, 'triples over odd y3..4095')
    print('exactly one q/u divisible3:', restricted_total, 'states; floor(y/3) for3-units, zero otherwise')
    print('minimum y with two row representations but all3 precision branch types:', first)
    print('41 witnesses (q,r,u,a):', branches)
    print('41 row supports:', sorted(predicted_fibre(41)))
    print('PASS: explicit checks active under -O')


if __name__ == '__main__':
    main()
