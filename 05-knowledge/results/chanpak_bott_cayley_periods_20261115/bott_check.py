"""Direct checks of the r=1 boundary and r=2 equal-order hostile."""
from itertools import combinations


def incoming(n, code):
    adj = [0] * n
    for i, (a, b) in enumerate(combinations(range(n), 2)):
        if code >> i & 1:
            adj[a] |= 1 << b
        else:
            adj[b] |= 1 << a
    return adj


def subset_path_counts(adj):
    n = len(adj)
    H = [0] * (1 << n)
    end = [[0] * n for _ in H]
    H[0] = 1
    for S in range(1, len(H)):
        V = S
        while V:
            b = V & -V
            V -= b
            j = b.bit_length() - 1
            if b == S:
                end[S][j] = 1
            else:
                T = S ^ b
                A = adj[j] & T
                while A:
                    a = A & -A
                    A -= a
                    end[S][j] += end[T][a.bit_length() - 1]
            H[S] += end[S][j]
    return H


def submasks(S):
    T = S
    while T:
        yield T
        T = (T - 1) & S
    yield 0


def profiles(adj):
    H = subset_path_counts(adj)
    s = len(H) - 1
    q2 = [sum(H[A]*H[S ^ A] for A in submasks(S)) for S in range(len(H))]
    p2 = q2[s] - 2*H[s]
    q3 = sum(H[A] * q2[s ^ A] for A in range(len(H)))
    return H[s], p2, q3 - 3*p2 - 3*H[s]


def cyclic_join(A, B, C):
    offsets = [0, len(A), len(A) + len(B)]
    result = []
    Data = [A, B, C]
    for i in range(3):
        previous = (i - 1) % 3
        pre = ((1 << len(Data[previous])) - 1) << offsets[previous]
        result.extend((x << offsets[i]) | pre for x in Data[i])
    return result


def ck(got, want):
    if got != want:
        raise AssertionError((got, want))


if __name__ == '__main__':
    ac, bc = 38858915863, 15014846659
    A, B = incoming(9, ac), incoming(9, bc)
    T1, T2 = incoming(1, 0), incoming(2, 0)
    ck(profiles(A), (375, 10242, 95562))
    ck(profiles(B), (375, 10242, 98202))
    fa = subset_path_counts(cyclic_join(A, T2, T2))[-1]
    fb = subset_path_counts(cyclic_join(B, T2, T2))[-1]
    ck(fb - fa, 10560)
    ck((
        subset_path_counts(cyclic_join(T1, T1, T1))[-1],
        subset_path_counts(cyclic_join(T2, T1, T1))[-1],
    ), (3, 5))
    print('order-nine P1/P2 tie, P3 gap 2640; literal probe counts:', fa, fb)
    print('singleton H-match crosses by +2')
