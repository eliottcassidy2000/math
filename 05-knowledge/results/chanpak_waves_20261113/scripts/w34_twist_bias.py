import sys, random, numpy as np
sys.path.insert(0,'/Users/e/Documents/GitHub/math/04-computation/experiments/ribbon_spectrum_20261104')
from ribbon import tournament_from_word, is_strong
R='/Users/e/Documents/GitHub/math/04-computation/experiments/ribbon_spectrum_20261104/span2_intrep.txt'
tok=open(R).read().split(); pos=0
def nxt():
    global pos; v=tok[pos]; pos+=1; return int(v)
dim=nxt(); nl=nxt()
init=np.array([nxt() for _ in range(dim)],dtype=np.int64); fin=np.array([nxt() for _ in range(dim)],dtype=np.int64)
MI={}
for _ in range(nl):
    li=nxt(); MI[li]=np.array([[nxt() for _ in range(dim)] for _ in range(dim)],dtype=np.int64)
MOD=225
L=[frozenset(),frozenset([1]),frozenset([2]),frozenset([1,2])]
def Hmod(word):
    y=init%MOD
    for q,a in enumerate(word):
        a=frozenset(s for s in a if s<=q); y=(MI[L.index(a)]@y)%MOD
    return int(fin@y)%MOD
runs=['1.2'*1]
def s(word): return '.'.join(''.join(map(str,sorted(a))) or '0' for a in word)
FLOOR=['1.2.1.2.1.2','2.1.2.1.2.1','0.2.0.2.0.2','12.0.12.0.12.0']
rng=random.Random(5); n=int(sys.argv[1]); trials=int(sys.argv[2])
stats={'all':[0,0,0],'strong':[0,0,0],'strong, no floor run':[0,0,0]}
for t in range(trials):
    w=[L[rng.randrange(4)] for _ in range(n)]
    h=Hmod(w); st=is_strong(tournament_from_word(2,w))
    keys=['all']
    if st:
        keys.append('strong')
        ws='.'+s(w)+'.'
        if not any(('.'+f+'.') in ws for f in FLOOR): keys.append('strong, no floor run')
    for k in keys:
        stats[k][0]+=1; stats[k][1]+= (h%3==0); stats[k][2]+= (h%5==0)
for k,(c,z3,z5) in stats.items():
    print(f'{k:22s} n={n}: samples {c:7d}  P(3|H) = {z3/max(c,1):.4f} (uniform 0.3333)   P(5|H) = {z5/max(c,1):.4f} (uniform 0.2000)')
