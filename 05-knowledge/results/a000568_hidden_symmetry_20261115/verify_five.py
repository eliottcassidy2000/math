from itertools import product
from collections import Counter, defaultdict
from verify import hp

P=5
steps=[frozenset(s) for s in product((1,4),(2,3))]
masks=[frozenset(d for d in range(P) if (b>>d)&1) for b in range(1,31)]

def graph(a,b,d,extra=None):
    n=10+(extra is not None)
    out=[0]*n
    for u in range(n):
        for v in range(u+1,n):
            if v==10: forward=(u//5==extra)
            elif u//5==v//5: forward=(v-u)%5 in (a,b)[u//5]
            else: forward=(v-u)%5 in d
            if forward: out[u]|=1<<v
            else: out[v]|=1<<u
    return out

def hp_reverse(out):
    n=len(out)
    inc=[sum(1<<u for u in range(n) if out[u]>>v&1) for v in range(n)]
    # D[S][v] counts paths beginning at v, using the original outgoing arcs.
    # Equivalently grow on the left through each predecessor of v.
    D=[[0]*n for _ in range(1<<n)]
    for v in range(n): D[1<<v][v]=1
    for S in range(1,1<<n):
        for v,value in enumerate(D[S]):
            if not value: continue
            R=inc[v]&~S
            while R:
                b=R&-R
                R-=b
                D[S|b][b.bit_length()-1]+=value
    return sum(D[-1])

def canonical(a,b,d):
    def bits(s): return sum(1<<x for x in s)
    options=[]
    for A,B,D in ((a,b,d),(b,a,frozenset(range(5))-frozenset((-x)%5 for x in d))):
        for u in range(1,5):
            for t in range(5):
                options.append((bits((u*x)%5 for x in A),bits((u*x)%5 for x in B),
                                bits((u*x+t)%5 for x in D)))
    return min(options)

def main():
    hist=Counter()
    classes=defaultdict(set)
    minimum=10**9
    for a,b,d in product(steps,steps,masks):
        out=graph(a,b,d)
        h=hp(out)
        assert h==hp_reverse(out)
        hist[h]+=1
        if h%3: classes[h].add(canonical(a,b,d))
        for extra in (0,1):
            out=graph(a,b,d,extra)
            aug=hp(out)
            assert aug==hp_reverse(out)
            minimum=min(minimum,aug)
    assert hist=={8175:80,8245:40,8265:40,15395:80,15465:120,15505:40,15565:40,15745:40}
    assert {h:len(cc) for h,cc in classes.items()}=={8245:1,15395:2,15505:1,15565:1,15745:1}
    assert minimum==30035
    print("480 mixed C5 tables, 960 distinguished-vertex controls, exact normalizer orbit census passed.")

if __name__=="__main__": main()
