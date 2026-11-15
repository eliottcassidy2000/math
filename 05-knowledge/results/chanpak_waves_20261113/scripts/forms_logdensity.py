# Sum of 2/q over distinct carrier forms (a,q,s,e,t,f,n) of strong marked cores, from the wave-26 palette census
# files pal_m.txt (regenerate with: gentourng m | ./palette run m > pal_m.txt).
import sys
def parse(code,n):
    out=[0]*n; k=0
    for i in range(n):
        for j in range(i+1,n):
            if code[k]=='1': out[i]|=1<<j
            else: out[j]|=1<<i
            k+=1
    return out
def strong(out,n):
    def closure(adj):
        reach=1; fr=1
        while fr:
            nf=0
            for v in range(n):
                if (fr>>v)&1: nf|=adj[v]
            nf&=~reach; reach|=nf; fr=nf
        return reach==(1<<n)-1
    inn=[0]*n
    for v in range(n):
        for u in range(n):
            if (out[v]>>u)&1: inn[u]|=1<<v
    return closure(out) and closure(inn)
d=sys.argv[1] if len(sys.argv)>1 else '.'
forms={}
for m in (3,4,5,6):
    S=set(); cnt=0
    for line in open(f'{d}/pal_{m}.txt'):
        p=line.split(); a,q,s,e,t,f,n=map(int,p[:7]); code=p[7]
        if not strong(parse(code,m),m): continue
        cnt+=1; S.add((a,q,s,e,t,f,n))
    forms[m]=S
    print(f"order {m}: strong marked cores {cnt}, distinct forms {len(S)}, sum 2/q = {sum(2/x[1] for x in S):.3f}")
for name,ms in (("{3,4}",(3,4)),("{3,4,5}",(3,4,5)),("{6}",(6,))):
    U=set().union(*[forms[m] for m in ms]); print(f"core set {name}: distinct forms {len(U)}, sum 2/q = {sum(2/x[1] for x in U):.3f}")
