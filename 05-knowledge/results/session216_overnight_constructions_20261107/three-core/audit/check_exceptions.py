from itertools import permutations
from pathlib import Path
import json,sys
source=Path(sys.argv[1])
out=[]
for value,record in json.loads(source.read_text()).items():
 rows=record['rows'];n=len(rows)
 assert all(not(rows[i]>>i&1) for i in range(n))
 assert all(((rows[i]>>j&1)+(rows[j]>>i&1)==1) for i in range(n) for j in range(i))
 count=sum(all(rows[p[i]]>>p[i+1]&1 for i in range(n-1)) for p in permutations(range(n)))
 assert count==int(value)
 cycle=record['cycle'];assert sorted(cycle)==list(range(n))
 assert all(rows[cycle[i]]>>cycle[(i+1)%n]&1 for i in range(n))
 out.append(dict(value=int(value),order=n,permutation_count=count,cycle_valid=True))
print(json.dumps(out,indent=2))
