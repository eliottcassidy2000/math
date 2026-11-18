"""Literal skeletons, explicit masks, independent subset DP controls."""
from collections import Counter
from functools import cache
from random import Random

B = [(12,17,2,6,13),(8,17,3,6,13)]
D = (2,4,9,3)
RES = [(28,53,56,50,32,1),(14,60,56,48,33,1)]


def chk(ok, *diagnostic):
    if not ok:
        raise RuntimeError(diagnostic)


def skeleton(a, s, repeat):
    fixed = 31 ^ sum(1 << i for i in repeat)
    lim = fixed.bit_count() + s

    @cache
    def go(mask, last, slots, counts):
        if not slots:
            return Counter({counts: 1}) if mask == fixed and all(counts) else Counter()
        out = go(mask, -1, slots - 1, counts).copy()
        for i in range(5):
            if last != -1 and not (a[last] >> i & 1):
                continue
            if i in repeat:
                k = repeat.index(i)
                if counts[k] < lim:
                    new = counts[:k] + (counts[k]+1,) + counts[k+1:]
                    out.update(go(mask, i, slots, new))
            elif not (mask >> i & 1):
                out.update(go(mask | 1 << i, i, slots, counts))
        return out
    return go(0, -1, s, (0,)*len(repeat))


def blow(a, parts):
    starts = [0]
    for p in parts:
        starts.append(starts[-1]+len(p))
    return [v << starts[i] | sum(((1 << len(p))-1) << starts[j]
            for j,p in enumerate(parts) if a[i] >> j & 1)
            for i,b in enumerate(parts) for v in b]


def count(a):
    size = 1 << len(a)
    dp = [[0]*len(a) for _ in range(size)]
    for i in range(len(a)):
        dp[1 << i][i] = 1
    for S,row in enumerate(dp):
        for i,v in enumerate(row):
            if not v: continue
            A = a[i] & (size-1-S)
            while A:
                b = A & -A
                A -= b
                dp[S|b][b.bit_length()-1] += v
    return [1]+[sum(row) for row in dp[1:]]


def pair(a):
    h = count(a)
    return h[-1], sum(x*y for x,y in zip(h,reversed(h)))


def tt(n):
    return [sum(1 << j for j in range(i+1,n)) for i in range(n)]


def sample(n,r):
    a = [0]*n
    for i in range(n):
        for j in range(i+1,n):
            x,y = (i,j) if r.randrange(2) else (j,i)
            a[x] |= 1 << y
    return a


def main():
    for a,d in zip(B,(11,13)):
        chk(skeleton(a,1,(0,4)) ==
            {(1,1):d,(1,2):6,(2,1):6,(2,2):3},a)
    a,b = [skeleton(t,2,(4,)) for t in B]
    chk({k:b[k]-a[k] for k in a.keys()|b.keys() if a[k]!=b[k]}
        == {(1,):4,(2,):2,(3,):-2})
    chk(pair(D)==(5,32))
    chk([pair(a) for a in RES]==[(23,272),(23,296)])
    rng = Random(44)
    trials = [(sample(rng.randint(1,5),rng),sample(rng.randint(1,5),rng))
              for _ in range(60)]
    trials += [(tt(1),tt(k)) for k in range(1,10)]
    trials += [(tt(1),r) for r in RES]
    for x,y in trials:
        hx,px=pair(x)
        hy,py=pair(y)
        for t,sgn in zip(B,(-1,1)):
            v=blow(t,[x,[0],[0],[0],y])
            chk(count(v)[-1]==3*px*py+sgn*hx*hy,x,y,t)
    chk([pair(blow(t,[[0]]*4+[D])) for t in B]==[(187,3736),(197,3728)])
    print("PASS: skeleton coefficients, D and both center-switch profiles,",
          2*len(trials),"literal Hamiltonian counts and signed two-slot change")


if __name__=='__main__':
    main()
