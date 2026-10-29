from itertools import permutations
import json

def req(x,msg):
 if not x:raise RuntimeError(msg)
def edges(n):return {(a,b) for a in range(n) for b in range(n) if (b-a)%n in range(1,(n+1)//2)}
def partner(E,P):
 seams=set(zip(P,P[1:]));req(seams<=E,'guard')
 return {(b,a) if (a,b) not in seams else (a,b) for a,b in E}
def h(E,n):return sum(all(e in E for e in zip(p,p[1:])) for p in permutations(range(n)))
E=edges(5);P=(0,1,2,3,4);Q=(0,2,3,4,1)
a,b=partner(E,P),partner(E,Q)
req((h(a,5),h(b,5))==(11,15),'brute H')
# Padding gives a universal source in both partners, independently of all other arcs.
for n in range(5,13):
 x=n;old=E
 E=E|{(x,P[0])}|{(v,x) for v in range(n) if v!=P[0]}
 PP=(x,)+P;QQ=(x,)+Q
 aa,bb=partner(E,PP),partner(E,QQ)
 req(aa==a|{(x,v) for v in range(n)},'source sum P')
 req(bb==b|{(x,v) for v in range(n)},'source sum Q')
 a,b,P,Q=aa,bb,PP,QQ
print(json.dumps({'brute_permutations':240,'base_H':[11,15],'padding_n_through':13,'method':'edge sets, permutation count; exact source-join equality'}))
