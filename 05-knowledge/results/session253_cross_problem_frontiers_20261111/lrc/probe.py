import sys,json,functools
from fractions import Fraction
from pathlib import Path
_REL='05-knowledge/results/session251_creative_completion_20261111/ternary'
_ROOTS=list(Path(__file__).resolve().parents)
_COMPILER=next((p/_REL for p in _ROOTS if (p/_REL/'compiler.py').exists()),None)
if _COMPILER is None: raise RuntimeError('Run inside the repository with the Session251 compiler available')
sys.path.insert(0,str(_COMPILER))
from compiler import signature,g,count,generate
CORES=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
sigs=set(signature(c,o) for c in CORES for o in range(1<<len(c)))
def no11(n):return not (n&(n>>1))
@functools.lru_cache(None)
def safe(n):
 m=Fraction(1)
 for b in bin(n)[2:] if n else '':
  m=min(Fraction(3,2)*m-int(b),Fraction(1))
  if m<=0:return False
 return True
# NEW restrictions only; no inherited spectrum verification.
N=65535
for n in range(N+1):
 if no11(n):
  if 3*g(n)<n:raise RuntimeError(('coercivity',n))
  if not n%2 and 3*g(n)<2*n+1:raise RuntimeError(('even',n))
 if safe(n) and not no11(n):raise RuntimeError(('safe adjacency',n))
res={'coercivity_universe':N,'signature_count':len(sigs),'bounds':{}}
for bound in (255,1023,4095):
 sets={k:set() for k in ('full','one_no11','both_no11','one_safe','both_safe')}
 # Complete inherited canonical inverse enumerator. Axis reductions delete
 # trailing ones, preserving prefix-closed no11 and strict-safe restrictions.
 for s in sigs:
  for h,a,b in generate(s,bound):
   sets['full'].add(h)
   na,nb=no11(a),no11(b)
   if na or nb:sets['one_no11'].add(h)
   if na and nb:sets['both_no11'].add(h)
   sa,sb=safe(a),safe(b)
   if sa or sb:sets['one_safe'].add(h)
   if sa and sb:sets['both_safe'].add(h)
 res['bounds'][bound]={k:{'count':len(v),'first_missing':[x for x in range(1,bound+1,2) if x not in v][:16]} for k,v in sets.items()}
print(json.dumps(res,indent=2,sort_keys=True))
