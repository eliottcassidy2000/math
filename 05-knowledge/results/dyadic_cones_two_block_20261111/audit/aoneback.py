# Auditor's brute-force check of README section 4 (one back arc) by enumerating path covers explicitly.
from itertools import permutations
from functools import lru_cache
def tables(m,adj):
    F=1<<m
    end=[[0]*m for _ in range(F)]; st=[[0]*m for _ in range(F)]
    for a in range(m): end[1<<a][a]=1; st[1<<a][a]=1
    for S in range(1,F):
        for a in range(m):
            if not S>>a&1 or S==1<<a: continue
            P=S&~(1<<a)
            end[S][a]=sum(end[P][b] for b in range(m) if P>>b&1 and adj[b][a])
            st[S][a]=sum(st[P][b] for b in range(m) if P>>b&1 and adj[a][b])
    hs=[1]+[sum(end[S]) for S in range(1,F)]
    return end,st,hs
bad=0;cases=0
for m in range(3,10):
    for i in range(m):
        for j in range(i+2,m):
            k=j-i-1; s=m-1-j
            adj=[[a<b for b in range(m)] for a in range(m)]
            adj[i][j]=False; adj[j][i]=True
            end,st,hs=tables(m,adj); F=1<<m; full=F-1
            H=hs[full]
            Q=lambda e,f: sum(end[S][e]*st[full^S][f] for S in range(1,full) if S>>e&1 and not S>>f&1)
            U=[sum(hs[full^C]*st[C][y] for C in range(1,F) if C>>y&1) for y in range(m)]
            W=Q(i,j)+Q(j,i)
            ok = H==1+2**k and W==2**k
            for y in range(m):
                if y<i: f=2**(i-y-1+s+1)*(3**k+2**(k+1)-1)
                elif y==i: f=2**s*(2**(k+1)-1)
                elif y<j:
                    a=y-i-1;b=j-y-1; f=2**s*(2**(b+1)+2**(a+1)*3**b-(a==0))
                elif y==j: f=2**s*(2**k+1)
                else: f=2**(m-1-y)*(2**k+1)
                ok = ok and U[y]==f
            cases+=1; bad+= (not ok)
            # classification of power-of-two responses for |O|=1
            for y in range(m):
                u=U[y]
                if u&(u-1)==0:
                    print('pow2 response: m',m,'i',i,'j',j,'k',k,'y',y,'offset',H,'U',u)
print('cases',cases,'bad',bad)
# number theory brute force
print([k for k in range(1,200) if (lambda x: x&(x-1)==0)(3**k+2**(k+1)-1)])
print([(a,b) for a in range(0,80) for b in range(0,80) if (lambda x: x>0 and x&(x-1)==0)(2**(b+1)+2**(a+1)*3**b-(a==0))])
