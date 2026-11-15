import pickle
from collections import Counter
D=pickle.load(open('memneed_1.765.pkl','rb')); need=D['need']
code=[(i,a) for i in range(3) for a in range(1,8)]; cidx={c:k for k,c in enumerate(code)}
syl=['0.2','1.2','0.123','2.13','0.1.3','0.3.1','0.13.2','0.1.123','0.3.12','0.12.3','1.2.3','1.12.3','0.1.23']
def mask(t): return sum(1<<(int(ch)-1) for ch in t if ch in '123')
def tocode(s):
    letters=[mask(t) for t in s.split('.')]; out=[]; z=0
    for l in letters:
        if l==0: z+=1
        else: out.append(cidx[(z,l)]); z=0
    assert z==0; return tuple(out)
SC=[tocode(s) for s in syl]
def cyc_parses(c):
    n=len(c)
    for start in range(n):
        w=c[start:]+c[:start]
        ok=[False]*(n+1); ok[0]=True
        for i in range(n):
            if ok[i]:
                for s in SC:
                    if tuple(w[i:i+len(s)])==s: ok[i+len(s)]=True
        if ok[n]: return True
    return False
tab=Counter()
for c,L in need.items(): tab[(min(L,13), cyc_parses(c))]+=1
print('cheap cycles by bottleneck length and whether they parse cyclically into floor syllables:')
for L in sorted(set(k[0] for k in tab)): print('  L=%2d: syllable-composed %6d, other %6d'%(L,tab[(L,True)],tab[(L,False)]))
print('syllable code forms:',dict(zip(syl,SC)))
