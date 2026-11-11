import json,sys,functools,time
from fractions import Fraction
from pathlib import Path
_REL='05-knowledge/results/session251_creative_completion_20261111/ternary'
_ROOTS=list(Path(__file__).resolve().parents)
_COMPILER=next((p/_REL for p in _ROOTS if (p/_REL/'compiler.py').exists()),None)
if _COMPILER is None: raise RuntimeError('Run inside the repository with the Session251 compiler available')
sys.path.insert(0,str(_COMPILER))
from compiler import signature,g,generate
CORES=[(0,),(2,4,1),(2,4,9,3),(2,4,9,19,7),(2,4,25,3,11)]
sigs=set(signature(c,o) for c in CORES for o in range(1<<len(c)))
@functools.lru_cache(None)
def safe(n):
 m=Fraction(1)
 for b in bin(n)[2:] if n else '':
  m=min(Fraction(3,2)*m-int(b),Fraction(1))
  if m<=0:return False
 return True
X=65535
# Minimum safe flank value, infinity if no safe-flank certificate.
best=[X*10]*(X+1)
full=bytearray(X+1)
none11=bytearray(X+1)
t=time.time()
for s in sigs:
 for h,a,b in generate(s,X):
  full[h]=1
  if not(a&(a>>1)) or not(b&(b>>1)):none11[h]=1
  if safe(a) and a<best[h]:best[h]=a
  if safe(b) and b<best[h]:best[h]=b
out={'bound':X,'full_count':sum(full),'one_safe_count':sum(full[h] and best[h]<X*10 for h in range(1,X+1,2)),'one_no11_extra_holes':[h for h in range(1,X+1,2) if full[h] and not none11[h]],'one_safe_extra_holes':[h for h in range(1,X+1,2) if full[h] and best[h]==X*10][:40], 'bounded_safe_flanks':{K:{'covered':sum(full[h] and best[h]<=K for h in range(1,X+1,2)), 'extra_holes':[h for h in range(1,X+1,2) if full[h] and best[h]>K][:16]} for K in (0,1,2,4,8,16,32,64)},'largest_minimum_safe_flank':max(v for h,v in enumerate(best) if full[h] and v<X*10)}
print(json.dumps(out,indent=2,sort_keys=True))
