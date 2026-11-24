from itertools import permutations,combinations
from functools import cache
from math import comb
from collections import defaultdict
from counts import flipped

@cache
def slots(T,v):
    a=[0]*len(T)
    def has(i,j): return bool(T[i]>>j&1)
    for p in permutations(i for i in range(len(T)) if i!=v):
        opt=has(v,p[0])+has(p[-1],v); need=0
        for i,j in zip(p,p[1:]):
            h=has(i,j)
            via=has(i,v) and has(v,j)
            if not h and not via:break
            if via:
                if h:opt+=1
                else:need+=1
        else:
            for z in range(opt+1):
                a[need+z]+=comb(opt,z)
    return tuple(a)

def is_strong(T):
    for rev in range(2):
        R={0}
        for _ in T:
            R|={j for i in R for j in range(len(T)) if
                (T[i]>>j&1 if not rev else T[j]>>i&1)}
        if len(R)!=len(T):return False
    return True

def scan():
    canonical=[(106,100,25,82,35,76,20),(114,100,89,19,34,76,24)]
    rails=defaultdict(list)
    for C in canonical:
        for q in [None,*combinations(range(7),2)]:
            T=flipped(C,*q) if q else C
            if not is_strong(T): continue
            for v in range(7):
                a=slots(T,v)
                if not any(a[2:]):continue
                plus,minus=[],[]
                for p in combinations((i for i in range(7) if i!=v),2):
                    U=flipped(T,*p)
                    b=slots(U,v)
                    if a[2:]!=b[2:] or abs(a[1]-b[1])!=2 or not is_strong(U):continue
                    tag=(q,v,p)
                    rails[(min(a[1],b[1]),a[2:])].append(tag)
                    (plus if b[1]>a[1] else minus).append(p)
                if plus and minus:print('two-sided fork:',q,v,a,plus,minus)
    print('rails from both charts',dict(rails))
if __name__=='__main__':scan()
