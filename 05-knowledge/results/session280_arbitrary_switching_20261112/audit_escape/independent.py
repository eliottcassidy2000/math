from functools import lru_cache
import json

def require(ok, message='independent audit gate'):
    if not ok:
        raise RuntimeError(message)

def count_profile(rows):
    n = len(rows)
    full = (1 << n) - 1
    require(all((not rows[i] >> i & 1 for i in range(n))))
    require(all(((rows[i] >> j & 1) + (rows[j] >> i & 1) == 1 for i in range(n) for j in range(i + 1, n))))
    reverse = tuple((sum((1 << j for j in range(n) if rows[j] >> i & 1)) for i in range(n)))

    @lru_cache(None)
    def starts(S, v, rev=False):
        if not S >> v & 1:
            return 0
        if S == 1 << v:
            return 1
        rest = S ^ 1 << v
        neigh = (reverse if rev else rows)[v] & rest
        ans = 0
        while neigh:
            bit = neigh & -neigh
            neigh ^= bit
            ans += starts(rest, bit.bit_length() - 1, rev)
        return ans

    @lru_cache(None)
    def H(S):
        return sum((starts(S, v) for v in range(n) if S >> v & 1)) if S else 1
    U = V = d = e = 0
    for S in range(1, full + 1, 2):
        weight = H(full ^ S)
        sign = 1 if S.bit_count() % 2 else -1
        u = starts(S, 0) * weight
        v = starts(S, 0, True) * weight
        U += u
        V += v
        d += sign * u
        e += sign * v
    R = full ^ 1
    Q = sum((H(S) * H(R ^ S) for S in range(full + 1) if S & 1 == 0))
    return dict(H=H(full), h=H(R), Q=Q, U=U, V=V, d=d, e=e)

def context(n, mask):
    rows = [0] * (n + 1)
    for i in range(n):
        if mask >> i & 1:
            rows[0] |= 1 << i + 1
        else:
            rows[i + 1] |= 1
        for j in range(i + 1, n):
            rows[i + 1] |= 1 << j + 1
    return tuple(rows)

def fuse(A, B):
    n = len(A)
    m = len(B)
    r = [0] * (n + m - 1)
    for X, mp in ((A, list(range(n))), (B, [0] + list(range(n, n + m - 1)))):
        for i in range(len(X)):
            for j in range(len(X)):
                if X[i] >> j & 1:
                    r[mp[i]] |= 1 << mp[j]
    for i in range(1, n):
        for j in range(n, n + m - 1):
            r[i] |= 1 << j
    return tuple(r)

def main():
    old = ((56, 13, 25, 48, 34, 6), (52, 13, 24, 49, 34, 6), (52, 13, 24, 33, 42, 6))
    repair = (386, 508, 489, 497, 485, 193, 257, 320, 32)
    outputs = [fuse(context(3, 1), x) for x in old]
    before = [count_profile(x) for x in outputs]
    after = count_profile(repair)
    require([p['H'] for p in before] == [269, 271, 273])
    require((after['H'], after['h'], after['V'], after['e']) == (273, 15, 441, 15))
    require(all(((p['h'], p['V'], p['e']) == (15, 441, 15) for p in before[:2])))
    diff = lambda a, b: [(i, j) for i in range(len(a)) for j in range(i + 1, len(a)) if (a[i] ^ b[i]) >> j & 1]
    repairs = diff(outputs[2], repair)
    gate = diff(outputs[1], repair)
    require(len(repairs) == 4 and len(gate) == 5)
    frame = (0, 1, 2, 3, 4, 6, 8, 5, 7)
    require(all((X[u] >> v & 1 for X in (outputs[1], repair) for u, v in zip(frame, frame[1:]))))
    return dict(before=before, repair=after, repair_edges=repairs, gate_edges=gate, frame=frame)
def H_only(rows):
    n=len(rows);N=1<<n;dp=[[0]*n for _ in range(N)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,N):
        for v in range(n):
            if not dp[S][v]:continue
            rem=rows[v]&~S
            while rem:
                bit=rem&-rem;rem-=bit
                dp[S|bit][bit.bit_length()-1]+=dp[S][v]
    return sum(dp[-1])

def carry_controls():
    old=((56,13,25,48,34,6),(52,13,24,49,34,6))
    bank=[fuse(context(3,1),x) for x in old]+[(386,508,489,497,485,193,257,320,32)]
    results=[]
    for P in(56,57):
        word=bin(P)[2:];mask=sum(1<<i for i,b in enumerate(word) if b=='1')
        L=context(len(word),mask)
        results.extend(H_only(fuse(L,T)) for T in bank)
    require(results==list(range(24875,24886,2)),"six literal counts")
    # Independently evaluate the arithmetic recurrence to extreme scales.
    @lru_cache(None)
    def g(P):
        if P==0:return 1
        if P%2:return g(P//2)
        z=P//2;return g(z)+z+z%2
    sample=[]
    for k in range(0,102,2):
        V=436*4**k+5
        require((V-21)%30==0,"even z exact divisibility")
        z=(V-21)//15;P=2*z
        def shift(p):return 15*(g(p)-1)+(V-15)*p
        require(shift(P+1)-shift(P)==6,"carry increment")
        h0=269+872*(4**k-1)//3+10*k
        if k in(0,2,10,100):sample.append({'k':k,'N':str(h0+shift(P))})
    return {'literal_six_counts':results,'all_even_k_0_through100_checked':51,'scale_samples':sample}

if __name__ == '__main__':
    print(json.dumps({"seed_controls":main(),"carry_controls":carry_controls()},indent=2))
