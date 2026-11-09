from lemmacheck import *
import random
random.seed(7)
def R(adj):
    f=flags(adj); p=f['pnu'][0]; return ph(adj,p,1)
def PH32(adj):
    f=flags(adj); w=f['r3'][0]; return ph(adj,w,2)
fails=0; tests=0; viol=0
# seeds: case A  K1 v (X u E2) u K1 ; case B  E2 v (Z u K1) u K1  for random small even X, Z
for trial in range(40):
    X=rand_even(random.randint(2,4))
    if X is None: continue
    A=pad(cone(pad(X,2)),1)      # Seed2
    B=pad(e2join(pad(X,1)),1)    # Seed3
    for name,S,op0 in (('A',A,R),('B',B,PH32)):
        fS=flags(S)
        need = (len(fS['pnu'])>=2) if name=='A' else (len(fS['r3'])>=1)
        if not need or not is_even(S): viol+=1; print('seed flag/even fail',name); continue
        n=len(S); off=dcount(S)-2**(n-1)
        G=op0(S); seq=[op0.__name__]
        for step in range(3):
            tests+=1
            fG=flags(G); n=len(G); d=dcount(G)
            if not(fG['good'] and is_even(G) and d-2**(n-1)==off): fails+=1; print('chain fail',name,seq,fG['good'],d-2**(n-1),off)
            if len(G)>=13: break
            op=random.choice([R,PH32]); G=op(G); seq.append(op.__name__)
print('chain tests',tests,'fails',fails,'seed violations',viol)
