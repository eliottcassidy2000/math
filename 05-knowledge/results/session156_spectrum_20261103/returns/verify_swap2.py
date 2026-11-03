def require(x):
 if not x: raise AssertionError("swap audit failed")
import sys,itertools
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parent))
import actual_swap_check as s
c=s.c
D=c.W['L2']+('b',)+c.W['L2']
for word in [('q1',),('q5','a'),('q5','a')+D,('q1','b')+D,('b','a'),('b','a','b')+D,('q5','a')+D+('b',),('q1','b')+D+('b',),('b','a','b')+D+('b',)]:
 rows=c.graph(word); n=len(rows)
 require(all(not rows[i]>>i&1 and rows[i] < 1<<n for i in range(n)))
 require(all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i,j in itertools.combinations(range(n),2)))
 # Independent graph queries for H,H_a,H_aa, not matrix recurrence.
 a=c.add_a(rows);aa=c.add_a(a)
 h=s.count(rows)[0];ha=s.count(a)[0];haa=s.count(aa)[0]
 require((ha-h)%2==0 and (haa-h)%2==0)
 state=(h,(ha-h)//2,(haa-ha)//2)
 print('word',word,'state',state,'port',n-1)
