import numpy as np, glob
V10=np.loadtxt('pairs_10.txt',dtype=np.int64)
S11=np.unique(np.concatenate([np.loadtxt(f,dtype=np.int64,usecols=(0,1)) for f in sorted(glob.glob('s11_*.txt'))]),axis=0)
print('order-11 distinct sampled pairs',len(S11))
def run(vals_lists, start, cls):
    U=np.unique(np.concatenate(vals_lists))
    N=np.arange(start,U[-1]+1); N=N[cls(N)]
    hit=np.isin(N,U); miss=N[~hit]
    return (int(miss[0]) if len(miss) else None), int(U[-1])
K1=lambda V:np.unique(3*V[:,0]+V[:,1]); TT2=lambda V:np.unique(5*V[:,0]+3*V[:,1])
tests=[('J1  K1@10 + K1@11, N = 0 mod 3',[K1(V10),K1(S11)],36835,lambda N:(N%2==1)&(N%3==0)),
       ('J2  TT2@10 + TT2@11, all odd N',[TT2(V10),TT2(S11)],188269,lambda N:N%2==1),
       ('J3  K1@10 + TT2@11, all odd N',[K1(V10),TT2(S11)],36835,lambda N:N%2==1),
       ('J3b K1@10 + K1@11 + TT2@11, all odd N',[K1(V10),K1(S11),TT2(S11)],36835,lambda N:N%2==1)]
for name,vl,start,cls in tests:
    first,top=run(vl,start,cls)
    print('%-42s every admissible N in [%d, %s] covered; first gap %s'%(name,start,(first-2) if first else top,first))
# order-11 sampled windows per kernel (longest covered run in the classes it realizes)
for nm,f in (('K1',K1),('TT2',TT2)):
    v=f(S11); res=sorted(set((v%6).tolist()))
    N=np.arange(v[0],v[-1]+1); N=N[np.isin(N%6,res)]; h=np.isin(N,v)
    idx=np.nonzero(~h)[0]; e=np.concatenate(([-1],idx,[len(N)])); g=np.diff(e)-1; k=int(np.argmax(g))
    print('order 11 sampled, %s: residues mod 6 %s, longest fully covered run [%d, %d]'%(nm,res,int(N[e[k]+1]),int(N[e[k+1]-1])))
