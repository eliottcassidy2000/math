import sys
def bm_profile(s):
    N=len(s); S=0
    for i,b in enumerate(s):
        if b: S|=1<<(N-1-i)
    C=1; B=1; L=0; m=1; prof=[]
    for n in range(N):
        d=bin((S>>(N-1-n)) & C).count('1')&1
        if d:
            if 2*L<=n:
                T=C; C^=B<<m; L=n+1-L; B=T; m=1
            else:
                C^=B<<m; m+=1
        else: m+=1
        prof.append(L)
    return prof
s=list(open(sys.argv[1],'rb').read()[:int(sys.argv[2])])
p=bm_profile(s); print("ref final L",p[-1],"n/2",len(s)//2)
open(sys.argv[3],'w').write(' '.join(map(str,p)))
