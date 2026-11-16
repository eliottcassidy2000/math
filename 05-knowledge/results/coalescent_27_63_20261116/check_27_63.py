import csv, sys
from collections import Counter, defaultdict
from itertools import permutations, product

def restriction(a,S):
    return tuple(sum(((a[i]>>j)&1)<<k for k,j in enumerate(S)) for i in S)

def profile(a):
    n=len(a); F=(1<<n)-1
    dp=[[0]*n for _ in range(F+1)]
    for i in range(n):dp[1<<i][i]=1
    for S in range(1,F+1):
        for i,x in enumerate(dp[S]):
            if x:
                d=a[i] & (F^S)
                while d:
                    b=d&-d;d^=b
                    dp[S|b][b.bit_length()-1]+=x
    return [1]+[sum(v) for v in dp[1:]]

def key(a):
    d=defaultdict(list)
    for i,x in enumerate(a):d[x.bit_count()].append(i)
    return min(restriction(a,[v for z in blocks for v in z])
               for blocks in product(*(permutations(d[k]) for k in sorted(d))))

def scc_sizes(a):
    rs=[]
    for i in range(len(a)):
        x=1<<i
        while True:
            p=x
            for j in range(len(a)):
                if p>>j&1:x|=a[j]
            if x==p:break
        rs.append(x)
    return tuple(len(v) for _,v in sorted(((k.bit_count(),[r for r in rs if r==k]) for k in set(rs)),reverse=True))

def main(path):
    seen={}; all27={}; pair={}
    for row in csv.DictReader(open(path)):
        h=int(row['H'])
        if h not in (27,63):continue
        a=tuple((int(row['code'])>>(8*i))&255 for i in range(8))
        tab=profile(a)
        assert tab[-1]==h
        for S in range(1,256):
            if tab[S]!=h and not (h==63 and tab[S]==27):continue
            z=restriction(a,[i for i in range(8) if S>>i&1]); k=key(z)
            if tab[S]==63:seen[k]=z
            elif h==27:all27[k]=z
    counts=Counter()
    for k,z in sorted(seen.items(),key=lambda x:(len(x[0]),x[0])):
        sizes=scc_sizes(z)
        if any(i>=2 for i in sizes) and sizes!=(len(z),):continue
        tab=profile(z)
        kids=defaultdict(list)
        for S,h in enumerate(tab):
            if h!=27:continue
            zz=restriction(z,[i for i in range(len(z)) if S>>i&1]); kk=key(zz)
            kids[kk].append(S)
            assert kk in all27
        counts[(len(z),bool(kids))]+=1
        if kids:
            print("host",z,"children",dict(kids))
            for kk,ss in kids.items():pair.setdefault(kk,[]).append((k,ss))
    print('all_target',counts,'all_source',Counter((len(k),scc_sizes(k)) for k in all27))
    print('distinct_classes_matched',Counter((len(k),scc_sizes(k)) for k in pair))
    print('distinct_pairs',Counter((len(k),len(z)) for k,v in pair.items() for z,_ in v))
if __name__=='__main__':main(sys.argv[1])
