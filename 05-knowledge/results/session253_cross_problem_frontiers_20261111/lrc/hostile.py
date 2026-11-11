import sys,json
from fractions import Fraction
from pathlib import Path
_REL='05-knowledge/results/session251_creative_completion_20261111/ternary'
_ROOTS=list(Path(__file__).resolve().parents)
_COMPILER=next((p/_REL for p in _ROOTS if (p/_REL/'compiler.py').exists()),None)
if _COMPILER is None: raise RuntimeError('Run inside the repository with the Session251 compiler available')
sys.path.insert(0,str(_COMPILER))
from compiler import signature,generate,hp_subsets,construct
CORES=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
N=52477
out=[]
def failure(n):
 w=bin(n)[2:] if n else ''
 for length in range(1,len(w)+1):
  for i in range(len(w)-length+1):
   v=sum(int(b)*Fraction(2,3)**(j+1) for j,b in enumerate(w[i:i+length]))
   if v>=1:return {'substring':w[i:i+length],'start':i,'sum':str(v)}
 return None
for name,c in zip(('K1','C3','S4','S5a','S5b'),CORES):
 for o in range(1<<len(c)):
  for h,a,b in generate(signature(c,o),N):
   if h==N:out.append({'core':name,'cut':o,'A':a,'B':b,'A_word':bin(a)[2:],'B_word':bin(b)[2:],'A_failure':failure(a),'B_failure':failure(b),'order':len(c)+1+a.bit_length()+b.bit_length()})
out.sort(key=lambda x:(x['order'],max(x['A'],x['B']),x['core'],x['cut']))
# Independent physical DP for a smallest-order certificate only.
w=out[0];c=CORES[('K1','C3','S4','S5a','S5b').index(w['core'])]
actual=hp_subsets(construct(c,w['cut'],w['A'],w['B']))[-1]
if actual!=N:raise RuntimeError(actual)
print(json.dumps({'target':N,'canonical_certificates':len(out),'minimum_order':out[0]['order'],'independent_smallest_graph_count':actual,'certificates':out},indent=2))
