import numpy as np
V={n:np.loadtxt('pairs_%d.txt'%n,dtype=np.int64).reshape(-1,2) for n in (8,9,10)}
for (al,be,name) in [(3,1,'K1'),(5,3,'TT2')]:
    for n in (8,9,10):
        s=al*V[n][:,0]+be*V[n][:,1]
        vals,cnt=np.unique(s,return_counts=True)
        # bulk: central 80% of pairs by statistic
        lo,hi=np.percentile(s,[10,90]).astype(np.int64)
        res=sorted(set((vals%6).tolist()))
        N=np.arange(lo,hi+1); N=N[np.isin(N%6,res)]
        r=np.zeros(len(N),dtype=np.int64); pos=np.searchsorted(vals,N); ok=(pos<len(vals))
        ok[ok]&=(vals[pos[ok]]==N[ok]); r[ok]=cnt[pos[ok]]
        pts=len(V[n]); lines=len(N)
        print('%s order %2d: central-80%% window [%d, %d]: %d admissible lines, %d pairs in all -> mean pairs/line %.2f; '
              'line multiplicity min %d, 1%% %d, median %d; uncovered lines %d'
              %(name,n,lo,hi,lines,pts,((s>=lo)&(s<=hi)).sum()/lines,r.min(),np.percentile(r,1),np.median(r),(r==0).sum()))
