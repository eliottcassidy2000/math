from collections import defaultdict
from itertools import combinations

def rows(n,code):
    a=[0]*n
    for k,(i,j) in enumerate(combinations(range(n),2)):
        if code>>k&1:a[i]|=1<<j
        else:a[j]|=1<<i
    return a

def opp(a):
    n=len(a);f=(1<<n)-1
    return [f^row^(1<<i) for i,row in enumerate(a)]

def dp(a):
    n=len(a); f=(1<<n)-1
    h=[1]+[0]*f
    ends=[[0]*n for _ in range(f+1)]
    for i in range(n): ends[1<<i][i]=1
    for s in range(1,f+1):
        h[s]=sum(ends[s])
        for i,x in enumerate(ends[s]):
            more=a[i]&~s
            while x and more:
                z=more&-more; more-=z
                ends[s|z][z.bit_length()-1]+=x
    return h,ends

def packets(a):
    n=len(a); f=(1<<n)-1
    h,end=dp(a); _,start=dp(opp(a))
    q=sum(h[s]*h[f^s] for s in range(f+1))
    bs=[sum(h[f^s]*start[s][i] for s in range(1,f+1)) for i in range(n)]
    bt=[sum(h[f^s]*end[s][i] for s in range(1,f+1)) for i in range(n)]
    pp=[1 if s.bit_count()%2 else -1 for s in range(f+1)]
    es=[sum(pp[s]*h[f^s]*start[s][i] for s in range(1,f+1)) for i in range(n)]
    et=[sum(pp[s]*h[f^s]*end[s][i] for s in range(1,f+1)) for i in range(n)]
    pair=[[sum(end[s][i]*start[f^s][j] for s in range(1,f)) for j in range(n)] for i in range(n)]
    for o in range(f+1):
        b=[i for i in range(n) if o>>i&1]
        R=(sum(bs[i] for i in b),sum(bt[i] for i in range(n) if i not in b))
        Y=(sum(es[i] for i in b),sum(et[i] for i in range(n) if i not in b))
        v=sum(bs[i]-end[f][i] for i in b)-sum(pair[i][j] for i in b for j in b)+h[f]
        yield (h[f],q,R,Y),v,o

def main(path):
    table=defaultdict(lambda:defaultdict(dict))
    for line in open(path):
        n,c=map(int,line.split())
        for packet,H,mask in packets(rows(n,c)):table[n][packet][H]=(c,mask)
    for n,buck in sorted(table.items()):
        W=[]
        for p,d in buck.items():
            for s in d:
                if s-2 in d:continue
                t=s
                while t+2 in d:t+=2
                W.append(((t-s)//2+1,s,t,p,[d[x] for x in range(s,t+1,2)]))
        top=max(W)
        print(n,top)
if __name__=='__main__':
    import sys;main(sys.argv[1])
