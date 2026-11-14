# L(A_x, k) for one-centre A_x (spine 0..m-1, hub u; u->p iff x[p]==1), root = spine vertex k.
from functools import lru_cache
def L_root(x,k):
    m=len(x); pos=[p for p in range(m) if p!=k]
    # case A: u in P. state (Pa_last in {-1,0,1}, Pb_started, Q_started)
    cnt={(-1,0,0):1}
    for p in pos:
        nc={}
        for (pa,pb,q),c in cnt.items():
            # to Pa
            s=(x[p],pb,q); nc[s]=nc.get(s,0)+c
            # to Pb
            if pb or x[p]==1:
                s=(pa,1,q); nc[s]=nc.get(s,0)+c
            # to Q
            if q or p>k:
                s=(pa,pb,1); nc[s]=nc.get(s,0)+c
        cnt=nc
    A=sum(c for (pa,pb,q),c in cnt.items() if pa in (-1,0))
    # case B: u in Q. state (Qa_last, Qb_started); P free
    cnt={(-1,0):1}
    for p in pos:
        nc={}
        for (qa,qb),c in cnt.items():
            s=(qa,qb); nc[s]=nc.get(s,0)+c            # to P
            if qa!=-1 or p>k:                          # to Qa (first element must be > k)
                s=(x[p],qb); nc[s]=nc.get(s,0)+c
            if qb or x[p]==1:                          # to Qb
                s=(qa,1); nc[s]=nc.get(s,0)+c
        cnt=nc
    B=sum(c for (qa,qb),c in cnt.items() if qa==0 or (qa==-1 and x[k]==0))
    return A+B
def W(x):
    m=len(x); s=0
    for i in range(m):
        if x[i]!=1: continue
        for j in range(i+1,m):
            if x[j]==0: s+= 1 if j-i==1 else 2**(j-i-2)
    return s
if __name__=="__main__":
    import itertools
    bad=0; tot=0
    for line in open('oc10.txt'):
        p=line.split(); bits=p[0]; H=int(p[1])
        n=0
        while n*(n-1)//2<len(bits): n+=1
        m=n-1
        # recover x: bit for pair (i, m) is '1' iff i->hub iff x_i=0
        idx=0; x=[None]*m
        for i in range(n):
            for j in range(i+1,n):
                if j==m: x[i]=0 if bits[idx]=='1' else 1
                idx+=1
        assert H==1+2*W(x)
        for it in p[2:]:
            v,r=it.split(':'); v=int(v); h,L,E,S=map(int,r.split(','))
            if v<m:
                tot+=1
                if L_root(x,v)!=L: bad+=1
    print("checked",tot,"spine roots; mismatches:",bad)
