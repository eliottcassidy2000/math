exec(open(__import__('os').path.join(__import__('os').path.dirname(__import__('os').path.abspath(__file__)),'w32_unipotent.py')).read().split("I=np.identity(dim,dtype=np.int64)")[0])
import sys
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation/experiments/ribbon_spectrum_20261104')
from ribbon import tournament_from_word
I=np.identity(dim,dtype=np.int64); p=PR[0]
def bulk_strong(ws):
    l=len(ws); k=max(6,30//l+4); word=[]
    for q in range(k*l):
        m=ws[q%l]; word.append(frozenset(s for s in (1,2,3) if (m>>(s-1))&1 and s<=q))
    adj=tournament_from_word(3,word); n=len(adj)
    def reach(adjm,src):
        seen=1<<src; fr=1<<src
        while fr:
            v=(fr&-fr).bit_length()-1; fr&=fr-1
            new=adjm[v]&~seen; seen|=new; fr|=new
        return seen
    radj=[0]*n
    for a in range(n):
        o=adj[a]
        while o:
            b=(o&-o).bit_length()-1; o&=o-1; radj[b]|=1<<a
    mid=n//2; comp=reach(adj,mid)&reach(radj,mid)
    return all((comp>>t)&1 for t in range(n//3,2*n//3))
ev1=0; ev1_strong=[]
for L in range(1,5):
    for w in itertools.product(masks,repeat=L):
        M=I.copy()
        for li in w: M=(MI[li]@M)%p
        if rank_mod((M-I)%p,p)==dim: continue
        ev1+=1
        if bulk_strong(list(w)): ev1_strong.append(w)
print('eigenvalue-1 blocks:',ev1,' of which bulk-strong:',len(ev1_strong), ev1_strong[:10])
