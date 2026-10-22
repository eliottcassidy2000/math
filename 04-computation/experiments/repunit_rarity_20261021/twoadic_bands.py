# Per fixed order j: residual of sigma(j,d) vs local mean, split by distance above the floor (d/f(j)) and by v2(d+1).
import json, math, random
from collections import defaultdict
tab={int(k):{int(a):b for a,b in v.items()} for k,v in json.load(open('strong_table_lt375.json')).items()}
F={5:9,6:15,7:25,8:45,9:75,10:125,11:225}
def v2(x):
    c=0
    while x%2==0: x//=2; c+=1
    return c
out=defaultdict(lambda: defaultdict(list))
for j in [6,7,8,9,10]:
    ds=[d for d in range(F[j],375,2)]
    cnt={d:tab.get(d,{}).get(j,0) for d in ds}
    top=max(d for d in ds if cnt[d]>0)
    ds=[d for d in ds if d<=top]
    for d in ds:
        w=[cnt[e] for e in ds if abs(e-d)<=max(6,0.12*d) and e!=d]
        if len(w)<4: continue
        r=math.log((cnt[d]+0.5)/(sum(w)/len(w)+0.5))
        ratio=d/F[j]; band='near (<2f)' if ratio<2 else ('mid (2-4f)' if ratio<4 else 'far (>4f)')
        out[band][min(v2(d+1),4)].append(r)
for band in ['near (<2f)','mid (2-4f)','far (>4f)']:
    g=out[band]
    print(f'{band:12s}: mean residual by v2(d+1): '+' '.join(f'{k}:{sum(v)/len(v):+.2f}(n={len(v)})' for k,v in sorted(g.items())))
