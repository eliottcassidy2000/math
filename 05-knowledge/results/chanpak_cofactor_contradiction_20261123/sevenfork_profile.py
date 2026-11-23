def f(address):
    total = 0
    while address:
        high, bit = divmod(address, 2)
        if bit == 0:
            total += high + (high & 1)
        address = high
    return total


def flipped(rows, i, j):
    out = list(rows)
    out[i] ^= 1 << j
    out[j] ^= 1 << i
    return tuple(out)


def star_rows(core, tags, x, y):
    """Ordinal left flank, marked core, right flank, and last vertex as apex."""
    left = bin(x)[2:] if x else ''
    right = bin(y)[2:][::-1].translate(str.maketrans('01', '10')) if y else ''
    tape = left + ''.join(str((tags >> i) & 1) for i in range(len(core))) + right
    n = len(tape)
    mask = (1 << n) - 1
    rows = [mask ^ ((1 << (i + 1)) - 1) for i in range(n)]
    width = len(core)
    core_mask = ((1 << width) - 1) << len(left)
    for i, row in enumerate(core):
        index = len(left) + i
        rows[index] = (rows[index] & ~core_mask) | (row << len(left))
    apex = 0
    for i, bit in enumerate(tape):
        if bit == '1':
            apex |= 1 << i
        else:
            rows[i] |= 1 << n
    return tuple(rows + [apex])


def subpaths(rows, starts=None):
    n = len(rows)
    size = 1 << n
    end = [[0] * n for _ in range(size)]
    for i in range(n):
        if starts is None or starts >> i & 1:
            end[1 << i][i] = 1
    for subset in range(1, size):
        for i, count in enumerate(end[subset]):
            if count:
                choices = rows[i] & ~subset
                while choices:
                    bit = choices & -choices
                    end[subset | bit][bit.bit_length() - 1] += count
                    choices ^= bit
    return [sum(row) for row in end]


def profile(core, tags):
    n = len(core)
    full = (1 << n) - 1
    paths = subpaths(core)
    paths[0] = 1
    h = paths[-1]
    q = sum(paths[s] * paths[full ^ s] for s in range(full + 1))

    def eo(graph, marked):
        rooted = subpaths(graph, marked)
        return tuple(sum(rooted[s] * paths[full ^ s] for s in range(1, full + 1)
                         if s.bit_count() % 2 == parity) for parity in (0, 1))

    reverse = tuple(row ^ full ^ (1 << i) for i, row in enumerate(core))
    k = subpaths(star_rows(core, tags, 0, 0))[-1]
    return h, q, eo(core, tags), eo(reverse, full ^ tags), k


def count_from_profile(data, x, y):
    h, q, (e, o), (end_e, end_o), k = data
    return (k + h * (f(x) + f(y)) + q * x * y + (end_e + end_o) * x +
            (e + o) * y - (end_e - end_o) * (x & 1) - (e - o) * (y & 1))



CENTER=(16,9,67,21,38,15,59)
TAG=68
def toggle(A,e):return flipped(A,*e)
def tagged(A,m):
 h,q,(a,b),(c,d),k=profile(A,m)
 return h,q,a+b,c+d,b-a,d-c,k
