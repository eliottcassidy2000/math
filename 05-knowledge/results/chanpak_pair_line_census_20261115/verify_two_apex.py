import itertools
import random
import subprocess
import sys
from reproduce import stats

exe=sys.argv[1]
r=random.Random(115203)
tests=0
for m in range(2,8):
    masks=list(itertools.product(range(1<<m),repeat=2)) if m<=3 else [
        (r.randrange(1<<m),r.randrange(1<<m)) for _ in range(35)
    ]
    for a,b in masks:
        c=sum(1<<j for j,(u,v) in enumerate(itertools.combinations(range(m+2),2))
                if u<m and ((v==m and a>>u&1) or (v==m+1 and b>>u&1)))
        h,q,_=stats(m+2,c)
        out=subprocess.check_output([exe,str(m),str(a),str(b)]).decode().strip()
        assert out==f"{h},{q}", (m,a,b,h,q,out)
        tests+=1
print("Direct subset Hamilton/cover controls passed:",tests)
