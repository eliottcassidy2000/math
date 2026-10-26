from itertools import permutations,combinations
import json

def s(n):return sum(d for d in range(1,n) if n%d==0)
def need(x,m):
 if not x:raise RuntimeError(m)
need([s(n) for n in [4,9,57,65,8,18,6,25]]==[3,4,23,19,7,21,6,6],'aliquot witnesses')
need(2%7==23%7 and 5%7==19%7,'prime target residue collisions')
need((2*2)%7==4 and (2*5)%7==3,'multiplier predicate hostile')
need(2**2%7==5**2%7==4,'square quotient collapse')
R={1,2,4};rows=[sum(1<<j for j in range(7) if (j-i)%7 in R) for i in range(7)]
H=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(6)) for p in permutations(range(7)))
Aut=sum(all((rows[i]>>j&1)==(rows[p[i]]>>p[j]&1) for i,j in combinations(range(7),2)) for p in permutations(range(7)))
need((H,Aut)==(189,21),'Paley direct control')
need(H%7==0 and H//Aut==9,'raw H and group order differ')
print(json.dumps({'status':'PASS','aliquot_witnesses':{str(n):s(n) for n in [4,9,57,65,8,18,6,25]},'Paley7_H':H,'Paley7_Aut':Aut,'Paley7_tc':H//Aut,'scope':'specific residue/action descent maps refuted; not all possible bridges'},sort_keys=True))
