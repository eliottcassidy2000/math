from itertools import combinations
from random import Random


def slots(r, mx):
    n = len(r)
    f = (1 << n) - 1
    a = [[0] * n for _ in range(f + 1)]
    H = [1] + [0] * f
    for s in range(1, f + 1):
        for i in range(n):
            if s >> i & 1:
                y = s ^ (1 << i)
                a[s][i] = sum(a[y][j] for j in range(n) if (r[i] & y) >> j & 1) if y else 1
        H[s] = sum(a[s])
    out = [[1] + [0] * f, H]
    for m in range(2, mx + 1):
        u = [0] * (f + 1)
        for s in range(f + 1):
            t = s
            while True:
                u[s] += H[t] * out[-1][s ^ t]
                if not t:
                    break
                t = (t - 1) & s
        out.append(u)
    return out


def fields(r, M, mx):
    n = len(r)
    f = (1 << n) - 1
    Q = slots(r, mx)
    a = [[[0] * n for _ in r] for _ in range(f + 1)]
    for s in range(1, f + 1):
        for i in range(n):
            if not (s & M) >> i & 1:
                continue
            y = s ^ (1 << i)
            # Build with the last vertex, so the initial mark is retained.
            if not y:
                a[s][i][i] = 1
            else:
                for k in range(n):
                    if not y >> k & 1:
                        continue
                    a[s][i][k] = sum(a[s ^ (1 << k)][i][j] for j in range(n) if r[j] >> k & 1)
    result = []
    for m in range(1, mx + 1):
        U = V = J = 0
        for s in range(1, f + 1):
            tot = sum(sum(x) for x in a[s]) * Q[m][f ^ s]
            U += tot
            V += (-1) ** s.bit_count() * tot
            if s.bit_count() % 2 == 0:
                J += 2 * Q[m][f ^ s] * sum(a[s][i][j] for i in range(n) for j in range(n) if M >> j & 1)
        result.append((Q[m][f], U, V, J))
    return result


def turnback(r, M, j, L):
    n = len(r)
    all_tail = ((1 << L) - 1) << n
    v = n + j - 1
    a = [x | all_tail for x in r]
    a += [((1 << (n + L)) - (1 << (n + i + 1))) for i in range(L)]
    for s in range(n):
        if M >> s & 1:
            a[s] ^= 1 << v
    a[v] |= M
    return a


def check(r, M, j, L, mx):
    F = fields(r, M, mx)
    for m, (base,U,V,J) in enumerate(F, 1):
        expected = m**L * base + m**(L-j+1) * ((m+1)**(j-1)*U + (m-1)**(j-1)*V) - m**L*J
        actual = slots(turnback(r, M, j, L), m)[m][-1]
        if expected != actual:
            raise AssertionError((r,M,j,L,m,F,expected,actual))


def main():
    r = (16,21,17,7,8)
    expect = {9:(25,4,122,-34,20), 20:(25,4,111,-31,20)}
    for M in (9,20):
        F = fields(r, M, 2)
        G = (F[0][1], F[0][3], *F[1][1:])
        if G != expect[M]:
            raise AssertionError((M,G))
        print('literal core kernel', M, G)
    for j in range(2,8):
        h1,q1 = [x[-1] for x in slots(turnback(r,9,j,j),2)[1:]]
        h2,q2 = [x[-1] for x in slots(turnback(r,20,j,j),2)[1:]]
        if (h1,h2,q1-q2)!=(5+25*2**(j-1),5+25*2**(j-1),22*3**(j-1)-6):
            raise AssertionError((j,h1,q1,h2,q2))
    print('literal same-H separation verified for j=2..7')
    G = Random(8121)
    ntests=0
    for trial in range(36):
        n=G.randrange(1,6)
        a=[0]*n
        for i,j in combinations(range(n),2):
            if G.randrange(2): a[i]|=1<<j
            else: a[j]|=1<<i
        M=G.randrange(1,1<<n)
        for L in (1,2,3):
            for j in range(1,L+1):
                check(a,M,j,L,3)
                ntests+=3
    print('literal returned-slot checks (all Q1,Q2,Q3):',ntests)

if __name__=='__main__':
    main()
