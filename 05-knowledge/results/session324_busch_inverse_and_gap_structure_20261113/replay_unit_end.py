"""Exhaust rooted labeled tournaments and verify the unit-end grammar through five."""
from itertools import combinations, permutations
from replay_busch_inverse import data, strong, need

def canon(T):
    return min(tuple(T[i][j] for a,i in enumerate(p) for j in p[a+1:]) for q in permutations(range(1,len(T))) for p in [(0,)+q])

def extend(C,b,w):
    c=len(C); n=1+b+c
    T=[[0]*n for _ in range(n)]
    for i in range(1,b+1):
        for j in range(i+1,n): T[i][j]=1
        T[0][i]=int(i==1 or ((w>>(i-2))&1)); T[i][0]=1-T[0][i]
    for i in range(c):
        for j in range(c): T[b+1+i][b+1+j]=C[i][j]
        T[b+1+i][0]=int(i==0); T[0][b+1+i]=int(i!=0)
    return T

def main():
    rows={1:[[[0]]],2:[]}
    for n in range(3,6):
        rows[n]=[]
        for c in range(1,n):
            b=n-1-c
            for C in rows[c]:
                for w in range(1<<max(0,b-1)):
                    T=extend(C,b,w); need(strong(T) and data(T,0)[2]==1)
                    H,_,_,_=data(C,0)
                    N,_,_,_=data(T,0)
                    covers=[(p,cut) for p in permutations(range(c)) for cut in range(c+1)
                            if all(C[p[i]][p[i+1]] for i in range(c-1) if i+1!=cut)]
                    W=sum(cut<c and p[cut]==0 for p,cut in covers)
                    Wdual=sum(cut>0 and p[cut-1]==0 for p,cut in covers)
                    A=[[0]+[int(k==0) for k in range(len(T))]]+[[int(k!=0)]+r for k,r in enumerate(T)]
                    HA,_,_,LA=data(A,0)
                    HP,_,EP,LP=data([r[:b+1] for r in T[:b+1]],0)
                    SC=data(tuple(zip(*C)),0)[2]
                    need(N==LP*Wdual+HP*H-EP*SC and HA==N-1+LP*len(covers)-W and LA-HA==1)
                    rows[n].append(T)
        decoded=set(map(canon,rows[n]))
        need(len(decoded)==len(rows[n]))
        truth=set()
        pairs=list(combinations(range(n),2))
        for m in range(1<<len(pairs)):
            M=[[0]*n for _ in range(n)]
            for i,(a,b) in enumerate(pairs):
                M[a][b]=m>>i&1; M[b][a]=1-M[a][b]
            if strong(M) and data(M,0)[2]==1: truth.add(canon(M))
        need(truth==decoded)
    need([len(rows[i]) for i in range(1,6)]==[1,0,1,3,8])
    print("unit-end grammar: exact labeled/canonical comparison through order five")
if __name__=="__main__": main()
