"""Exact spatial carry walls, division-budget coding, and temporal hostiles."""
from itertools import product


def need(ok, message):
    if not ok:
        raise RuntimeError(message)


def v2(n):
    return (n & -n).bit_length() - 1


def U(n):
    z = 3*n + 1
    return z >> v2(z)


def carry(c, bit):
    return (3*bit + c) // 2


def first_wall(n):
    previous = n & 1
    for s in range(1, n.bit_length()+2):
        bit = (n >> s) & 1
        if bit == previous:
            return s
        previous = bit
    raise RuntimeError(('no finite wall', n))


def main():
    word_count = 0
    print('Spatial carry scans: all binary words of lengths 1..16')
    for length in range(1, 17):
        nonsync = 0
        for word in product((0, 1), repeat=length):
            states = {0, 1, 2}
            repeated = False
            for i, bit in enumerate(word):
                old_size = len(states)
                states = {carry(c, bit) for c in states}
                need(len(states) <= old_size, ('rank increases', word, i))
                repeated |= i > 0 and bit == word[i-1]
                need((len(states) == 1) == repeated, ('wall equivalence', word, i))
                if not repeated:
                    need(states == ({0, 1} if bit == 0 else {1, 2}), ('alternating bank', word, states))
            nonsync += len(states) != 1
            word_count += 1
        need(nonsync == 2, ('nonsync census', length, nonsync))
    print('words checked:', word_count, '; exactly two nonsynchronizing words at every length')

    print('Physical odd sources: every 1 <= n < 2^18')
    source_count = 0
    parity_counts = [0, 0]
    alternating_count = 0
    max_a = 0
    for n in range(1, 1 << 18, 2):
        a = first_wall(n)
        need(a == v2(3*n+1), ('division identity', n, a))
        h = n >> (a+1)
        epsilon = 1 if a % 2 == 0 else 5
        residue = ((1 if a % 2 == 0 else 5)*(1 << a)-1)//3
        need(n == residue + (h << (a+1)), ('wall cylinder', n, a, h))
        need(U(n) == 6*h+epsilon, ('tail identity', n, a, h))
        c = 1
        for i in range(a+1):
            z = 3*((n >> i) & 1)+c
            need(z % 2 == (i == a), ('physical emission', n, a, i, z))
            c = z//2
        need(c == (0 if a % 2 == 0 else 2), ('outgoing wall carry', n, a, c))
        bits = [(n >> i) & 1 for i in range(n.bit_length())]
        alternating = all(bits[i] != bits[i-1] for i in range(1, len(bits)))
        need(alternating == (U(n) == 1), ('direct basin', n))
        need((a % 2 == 0) == (U(n) % 6 == 1), ('wall orientation', n))
        need(n == 1 or (U(n) < n) == (a >= 2), ('descent criterion', n))
        parity_counts[a % 2] += 1
        alternating_count += alternating
        max_a = max(max_a, a)
        source_count += 1
    print('sources checked:', source_count, '; even/odd wall-index counts:', parity_counts,
          '; direct-basin alternating sources:', alternating_count, '; maximum a:', max_a)

    print('Gilbreath control: every tail over {0,2,4} of lengths 1..9')
    triangles = 0
    protected = destroyed = no_four = 0
    for length in range(1, 10):
        for tail in product((0, 2, 4), repeat=length):
            row = [1]+list(tail)
            columns = [[] for _ in row]
            leading = []
            while row:
                leading.append(row[0])
                for i, value in enumerate(row):
                    columns[i].append(value)
                row = [abs(row[i]-row[i+1]) for i in range(len(row)-1)]
            for i, column in enumerate(columns[1:], start=1):
                seen_two = False
                for value in column:
                    need(not (seen_two and value == 4), ('phase reversal', tail, i, column))
                    seen_two |= value == 2
            if 4 in tail:
                F = tail.index(4)+1
                wall = 2 in tail[:F-1]
                if wall:
                    need(all(x == 1 for x in leading), ('protected Gilbreath edge', tail, leading))
                    protected += 1
                else:
                    need(leading[F] == 3 and all(x == 1 for x in leading[:F]), ('destroyed edge', tail, leading))
                    destroyed += 1
            else:
                need(all(x == 1 for x in leading), ('no-four edge', tail))
                no_four += 1
            triangles += 1
    print('triangles checked:', triangles, '; protected:', protected, '; destroyed:', destroyed,
          '; no four:', no_four)
    hostile = [1, 2, 6]
    row1 = [abs(hostile[i]-hostile[i+1]) for i in range(2)]
    need(row1 == [1, 4] and abs(row1[0]-row1[1]) == 3, 'alphabet-six hostile')
    print('Minimal alphabet-six breach:', hostile, '->', row1, '-> [3]')

    # A spatially synchronized canonical input can have an unsynchronized
    # canonical output. Thus synchronization is not a temporal phase.
    need(first_wall(3) == 1 and U(3) == 5 and first_wall(5) == 4, '3 to 5 control')
    print('Minimal odd canonical-word phase reversal: 3 (11) -> 5 (101); first walls 1 -> 4')

    print('Abundant-wall hostile: n_k=4*8^k-5, k=1..128')
    for k in range(1, 129):
        n = 4*8**k-5
        bits = [(n >> i) & 1 for i in range(n.bit_length())]
        equal_pairs = sum(bits[i] == bits[i-1] for i in range(1, len(bits)))
        need(equal_pairs == 3*k-1, ('pair count', k, equal_pairs))
        need(first_wall(n) == 1, ('fixed earliest wall', k))
        x = n
        for j in range(k):
            need(x == 4*9**j*8**(k-j)-5, ('even shadow', k, j))
            x = U(x)
            need(x == 6*9**j*8**(k-j)-7 and x > n, ('odd shadow', k, j))
            x = U(x)
            need(x == 4*9**(j+1)*8**(k-j-1)-5 and x > n, ('next even shadow', k, j))
    print('all controls passed; 3k-1 internal walls and first wall at1 coexist with tau(n_k)>2k')

    print('Exact suffix universality: a=1..64, h=0..255')
    for a in range(1, 65):
        epsilon = 1 if a % 2 == 0 else 5
        residue = (epsilon*(1 << a)-1)//3
        for h in range(256):
            n = residue + (h << (a+1))
            need(first_wall(n) == a and U(n) == 6*h+epsilon, ('unbounded-tail cylinder', a, h))
    print('16384 branch/tail controls passed; no tail-size inference follows from a fixed wall')


if __name__ == '__main__':
    main()
