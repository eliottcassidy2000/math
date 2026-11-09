"""THM-5056: explicit odd realisers of the offset values 2^(L-1) + c (c in {1, 5, 37}) and of the cone-chain values
2^(L-1) + 2^(L-3) + c for 10 <= L <= 20, built by the induction of oddbuild.py.  Output in wcheck format (set O).
usage: python3 oddspecial.py out.txt ; ./wcheck < out.txt"""
import sys
sys.argv=[sys.argv[0],'3','5','0']+sys.argv[1:]
src=open('oddbuild.py').read().split("lo,hi=int(sys.argv[1])")[0]
exec(src)
out=open(sys.argv[4],'w'); cnt=0
for Lk in range(10,21):
    for c in (1,5,37):
        for kk in (2**(Lk-1)+c, 2**(Lk-1)+2**(Lk-3)+c if Lk>=13 else None):
            if kk is None: continue
            G=build(kk); n=len(G); assert n==L(kk)
            edges=' '.join(f'{u}-{v}' for u in range(n) for v in range(u+1,n) if (G[u]>>v)&1)
            out.write(f'odd{kk} O {n} {kk} | {edges}\n'); cnt+=1
out.close(); print('special witnesses written:',cnt)
