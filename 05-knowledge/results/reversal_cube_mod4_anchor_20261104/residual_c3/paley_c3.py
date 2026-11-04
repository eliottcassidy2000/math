import itertools, json, time
from rw import F_runword
from lib import conv
def paley7():
    return [[1 if i!=j and (j-i)%7 in (1,2,4) else 0 for j in range(7)] for i in range(7)]
C3=[[0,1,0],[0,0,1],[1,0,0]]
def v2(x):
    if x==0: return 99
    k=0
    while x%2==0: x//=2; k+=1
    return k
Q=paley7(); j=0; m=7
blocks={i:C3 for i in range(m) if i!=j}
nv=6+3  # variables: 6 blocks then 3 named
D=[()]+[(i,) for i in range(nv)]+list(itertools.combinations(range(nv),2))
def run(Qx, named):
    jobs=[]
    for S in D:
        for R in D:
            a=[1 if i in S else 0 for i in range(nv)]; c=[1 if i in R else 0 for i in range(nv)]
            jobs.append((Qx,j,named,blocks,a,c))
    t=time.time(); vals=F_runword(jobs); print('ran',len(jobs),'queries in',round(time.time()-t,1),'s',flush=True)
    F={(S,R):v for (S,R),v in zip(itertools.product(D,D),vals)}
    def sub(S):
        out=[]
        for r in range(len(S)+1):
            for T in itertools.combinations(S,r): out.append((T,(-1)**(len(S)-r)))
        return out
    return {(S,R):sum(sa*sc*F[(Sa,Rc)] for Sa,sa in sub(S) for Rc,sc in sub(R)) for S in D for R in D}
CT=run(Q,C3); CU=run(conv(Q),C3); CR=run(conv(Q),conv(C3))
for name,CX in (('unchanged C3',CU),('reversed C3',CR)):
    diff={k:CT[k]-CX[k] for k in CT if CT[k]!=CX[k]}
    vals=sorted({v2(x) for x in diff.values()})
    print(name,': nonzero exact coefficient differences',len(diff),'2-adic valuations present',vals[:6])
    ex=sorted(diff.items(),key=lambda kv:(v2(kv[1]),len(kv[0][0])+len(kv[0][1])))[:6]
    print('  examples (S=a-vars,R=c-vars; vars 0-5 blocks 1..6, 6-8 named):',ex)
print('H(T)=',CT[((),())],'H(U)=',CU[((),())])
