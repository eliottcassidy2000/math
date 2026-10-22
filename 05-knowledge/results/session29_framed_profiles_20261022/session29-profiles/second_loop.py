import audit as a
import itertools as it,json,math
from pathlib import Path
# Independent brute permutation Hamiltonian witness verification.
checks={}
for m in (40,41):
 paths=[p for p in it.permutations(range(5)) if all(a.arc(m,5,p[i],p[i+1]) for i in range(4))]
 aut=[p for p in it.permutations(range(5)) if a.induce(m,5,p)==m]
 a.require(len(paths)==a.hp(m,5),'DP versus brute Hamiltonian')
 checks[m]={'arcs':[[i,j] if a.arc(m,5,i,j) else [j,i] for i,j in a.edges(5)],'paths':paths,'automorphisms':aut}
# explicit deletion isomorphisms, independently permuting full card arrays
cards={m:[a.induce(m,5,[w for w in range(5) if w!=v]) for v in range(5)] for m in (40,41)}
for perm in it.permutations(range(5)):
 maps=[]
 for v,w in enumerate(perm):
  good=next((p for p in it.permutations(range(4)) if a.induce(cards[41][w],4,p)==cards[40][v]),None)
  if good is None:break
  maps.append(good)
 else:break
a.require(len(maps)==5,'deck matching')
checks['matching']={'deleted_vertex_map':perm,'card_local_isomorphisms':maps}
def join(ms):
 n=5*len(ms)
 return sum((a.arc(ms[i//5],5,i%5,j%5) if i//5==j//5 else 1)<<b for b,(i,j) in enumerate(a.edges(n)))
words={}
for word in it.product((40,41),repeat=2):
 m=join(word); p=a.profile(m,10,4);s=a.score(m,10)
 words[str(word)]={'D4':p,'scores':s,'H':a.hp(m,10)}
 a.require(a.hp(m,10)==a.hp(word[0],5)*a.hp(word[1],5),'join product')
a.require(len(set(str(x['D4']) for x in words.values()))==1,'join D4 common')
a.require(len(set(str(x['scores']) for x in words.values()))==1,'join scores common')
checks['join_words']=words
Path(__file__).with_name('second-results.json').write_text(json.dumps(checks,indent=2)+'\n')
print('SECOND LOOP PASS; direct paths, explicit deck matching, all 4 two-block joins')
# Correct automorphism and converse sidecars, plus Burnside word tests.
for m,autn in ((40,3),(41,1)):
 a.require(len(checks[m]['automorphisms'])==autn,'core automorphism order')
 a.require(a.L[5][m^1023]==m,'self converse core')
for r in range(1,9):
 ws=list(it.product((0,1),repeat=r)); orbits={min(w,w[::-1]) for w in ws}
 a.require(len(orbits)==(2**r+2**((r+1)//2))//2,'merged word count')
 for weights in ((15,13),(5,13)):
  mass=lambda w:math.prod(weights[i] for i in w)
  total=sum(mass(w) for w in orbits)
  p=(weights[0]**2+weights[1]**2)**(r//2)*(sum(weights) if r%2 else 1)
  a.require(2*total==sum(weights)**r+p,'weighted Burnside')
print('AUTOMORPHISM AND BURNSIDE CHECKS PASS r=1..8')
