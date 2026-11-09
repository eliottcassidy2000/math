"""THM-5056 seeds: the four odd seeds (offsets 5, 37, 37, 1), their first PH(3,2) step and two further offset moves,
checked by brute force (d by enumeration; oddness by a persistent adjacent-twin pair and by Aut(G))."""
import sys
sys.path.insert(0,'.')
from lemmacheck import dcount, automorphisms, is_even, pad, cone, ph   # guarded import: no tests run
from collections import defaultdict
def L(k): return k.bit_length()
def union(a,b): na=len(a); return a+[x<<na for x in b]
def K(n): return [((1<<n)-1)^(1<<v) for v in range(n)]
def E(n): return [0]*n
def join(a,b):
    na,nb=len(a),len(b); fa=(1<<na)-1; fb=((1<<nb)-1)<<na
    return [x|fb for x in a]+[(y<<na)|fa for y in b]
def threshold(k):   # threshold graph realising odd k on L(k) vertices (Lemma B1): bit i = 1 -> vertex i dominates 0..i-1
    adj=[0]
    for i in range(1,L(k)):
        if (k>>i)&1: adj=cone(adj)
        else: adj=pad(adj)
    return adj
def odd_by_twins(adj,x,y):
    return (adj[x]|(1<<x))==(adj[y]|(1<<y)) and (adj[x]>>y)&1
def flags(adj):
    n=len(adj); deg=[bin(a).count('1') for a in adj]
    return deg
# ---- (2) seeds
P3=[0b110,0b001,0b001]   # centre 0
K4e=[0b1110,0b1101,0b0011,0b0011]  # 0,1 adjacent twins of degree 3; 2,3 degree 2
C4=[0b1010,0b0101,0b1010,0b0101]
seeds={}
seeds[5]=(K4e,(0,1),2)                                     # order 4, d = 13 = 2^3 + 5; w = 2
s1=join(K(2),union(P3,E(2)))                               # K_2 v (P_3 u E_2): order 7, d = 101 = 2^6 + 37
seeds['37a']=(s1,(0,1),2+3)                                # w = an E_2 vertex (index 2 + 3 + 0 = 5)
g55=threshold(55)
s2=union(K(2),g55)                                         # K_2 u G_55: order 8, d = 165 = 2^7 + 37
w55=[v for v in range(len(g55)) if bin(g55[v]).count('1')==2][0]+2
seeds['37b']=(s2,(0,1),w55)
s3=union(K(2),C4)                                          # K_2 u C_4: order 6, d = 33 = 2^5 + 1
seeds[1]=(s3,(0,1),2)
ok=True
for name,(G,(x,y),w) in seeds.items():
    n=len(G); d=dcount(G); deg=flags(G)
    assert deg[w]==2 and w not in (x,y), name
    assert odd_by_twins(G,x,y), name
    H=ph(G,w,2); dH=dcount(H)
    c=d-2**(n-1)
    good = odd_by_twins(H,x,y) and dH==2**(len(H)-1)+c and not is_even(H)
    print('seed',name,'order',n,'d',d,'offset',c,'-> PH(3,2): order',len(H),'d',dH,'odd',not is_even(H),'twins kept',odd_by_twins(H,x,y))
    ok &= good
    # two more moves: R excluding a hub leaf, then PH(3,2) excluding a fresh r3 vertex outside {x,y}
    cur=H
    for step in range(2):
        degc=[bin(a).count('1') for a in cur]; n2=len(cur)
        if step==0:
            p=[v for v in range(n2) if degc[v]==1 and v not in (x,y)][0]; nxt=ph(cur,p,1)
        else:
            w2=[v for v in range(n2) if degc[v]==2 and v not in (x,y)][0]; nxt=ph(cur,w2,2)
        dn=dcount(nxt)
        good2 = odd_by_twins(nxt,x,y) and dn==2**(len(nxt)-1)+c
        print('   step',step,'order',len(nxt),'d',dn,'offset',dn-2**(len(nxt)-1),'twins kept',odd_by_twins(nxt,x,y))
        ok&=good2; cur=nxt
print('seeds and chains OK' if ok else 'SEED FAILURE')
