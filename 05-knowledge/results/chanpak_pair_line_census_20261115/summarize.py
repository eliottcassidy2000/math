"""Tight projected coverage percentages/gaps for the selected bands."""
import sys,json
from score import band
rows=[list(map(int,l.split())) for l in open(sys.argv[1])]
subsets={'all masks':(1,),'masks<=63 + anti':(2,3),'all masks + anti':(1,3),'H+Q':(4,),'all anti':(3,)}
out={}
for lo in (69,32769):
 y=[a for a in rows if a[0]>=lo]
 row={}
 for name,ks in subsets.items():
  A={a[0] for a in y if any(a[k] for k in ks)}
  row[name]={'any order':band(A,lo,59999),
       'through 10':band({a[0] for a in y if any(a[k]&2047 for k in ks)},lo,59999),
       'per order':{n:sum(any(a[k] & 1<<n for k in ks) for a in y) for n in range(1,16)}}
 out[lo]=row
print(json.dumps(out,indent=2))
