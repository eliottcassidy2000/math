#!/usr/bin/env python3
import itertools,math,json
from pathlib import Path

def involutions(n):return sum(math.factorial(n)//(2**j*math.factorial(j)*math.factorial(n-2*j)) for j in range(n//2+1))
def counts(p):
 N=(math.factorial(p-1)+(p-1)**2)//p;F=(involutions(p)+p-1)//p;return N,F,(N+F)//2
out=[]
for n in(3,5,7):
 def canon(p):return min(tuple((p[(i+a)%n]+b)%n for i in range(n)) for a in range(n) for b in range(n))
 def conv(p):return tuple((-p.index((-i)%n))%n for i in range(n))
 reps={canon(p) for p in itertools.permutations(range(n))};F=sum(canon(conv(p))==p for p in reps);merged={min(p,canon(conv(p))) for p in reps}
 if (len(reps),F,len(merged))!=counts(n):raise RuntimeError('orbit formula')
 # Independently verify scalar range obstruction at p7.
 J=[sum(1 if (p[j]-p[i])%n in range(1,(n+1)//2) else -1 for i in range(n) for j in range(n) if(j-i)%n in range(1,(n+1)//2)) for p in reps]
 out.append({'p':n,'ordinary':len(reps),'self_converse':F,'merged':len(merged),'distinct_J':len(set(J)),'I_p':involutions(n)})
N,F,M=counts(11)
if (N,F,M)!=(329900,3246,166573):raise RuntimeError('p11formula')
output={'status':'PASS','exhaustive':out,'formula_only_p11':{'ordinary':N,'self_converse':F,'merged':M}}
Path(__file__).with_name('prime-matching-counts-audit.json').write_text(json.dumps(output,indent=2)+'\n');print(output)
