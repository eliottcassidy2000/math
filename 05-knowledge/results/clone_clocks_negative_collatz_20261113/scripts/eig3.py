import sys, itertools
sys.path.insert(0, 'scripts')
from funnel import word_polys
import sympy as sp
def mats(w):
    t=len(w)
    L1,D1,B1=word_polys(w,1); L2,D2,B2=word_polys(w,2); L3,D3,B3=word_polys(w,3)
    z=1 if t==0 else 0
    return sp.Matrix([[B1,-D1,L2],[L1,-z,0],[D2,0,2**t]]), sp.Matrix([[2*B2,-2*D2,2*L3],[L2,-1,0],[D3,0,3**t]])
words=['']
for L in range(1,6):
    for bits in itertools.product('01',repeat=L):
        w=''.join(bits)
        if w[-1]=='1': words.append(w)
cache={w:mats(w) for w in words}
hits=[]
for m in (1,2,3):
    for seq in itertools.product(words, repeat=m):
        if m>=2 and sum(len(w)+1 for w in seq)>9: continue
        A1=sp.eye(3); A2=sp.eye(3)
        for w in seq:
            M1,M2=cache[w]; A1=M1*A1; A2=M2*A2
        for lam in [2]+[3**a for a in range(1,6)]:
            d1=(A1-lam*sp.eye(3)).det(); d2=(A2-lam*sp.eye(3)).det()
            if d1==0 or d2==0:
                hits.append((m,seq,lam,d1==0,d2==0))
        nblocks = globals().get('nblocks', 0) + 1; globals()['nblocks'] = nblocks
both=[h for h in hits if h[3] and h[4]]
print('blocks examined:', globals().get('nblocks'), ' eigenvalue-2 hits:', sum(1 for h in hits if h[2]==2))
print('blocks with 3^a eigenvalue in block1:',sum(1 for h in hits if h[3]),' block2:',sum(1 for h in hits if h[4]),' both:',len(both))
for h in both[:40]: print(h)
for h in [h for h in hits if h[4]][:40]: print('blk2', h)
