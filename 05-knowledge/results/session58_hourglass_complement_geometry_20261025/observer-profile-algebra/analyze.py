import itertools,json,sys
from pathlib import Path

def req(v,m):
 if not v:raise RuntimeError(m)
def analyze(n):
 pairs=list(itertools.combinations(range(n),2));idx={p:i for i,p in enumerate(pairs)};maps=[]
 for p in itertools.permutations(range(n)):
  dst=[idx[tuple(sorted((p[i],p[j])))] for i,j in pairs]
  inv=sum(1<<k for k,(i,j) in enumerate(pairs) if p[i]>p[j])
  maps.append((dst,inv))
 seen=set();rows=[]
 for x in range(1<<len(pairs)):
  if x in seen:continue
  bits=[i for i in range(len(pairs)) if x>>i&1];orb=set();signed=0;even=True;aut=0
  for dst,inv in maps:
   y=sum(1<<dst[k] for k in bits);orb.add(y);s=(-1)**((x&inv).bit_count()%2);signed+=s
   if y==x:
    aut+=1
    if s<0:even=False
  seen.update(orb)
  if even and len(bits)%2==0:
   rows.append({'mask':x,'edges':[pairs[k] for k in bits],'aut':aut,'signed_order_sum':signed,'orbit_size':len(orb)})
 out={'n':n,'dimension_converse_invariants':len(rows),'rank':sum(r['signed_order_sum']!=0 for r in rows),'rows':rows}
 return out
if __name__=='__main__':
 out=[analyze(n) for n in range(1,int(sys.argv[1] if len(sys.argv)>1 else 6)+1)]
 expected=[(1,1),(1,1),(2,2),(3,3),(10,10),(34,32)]
 for d,e in zip(out,expected):req((d['dimension_converse_invariants'],d['rank'])==e,'expected dimensions/ranks')
 Path(__file__).with_name('walsh_certificate.json').write_text(json.dumps(out,indent=2)+'\n')
 for d in out:
  print(d['n'],d['dimension_converse_invariants'],d['rank'],[r for r in d['rows'] if not r['signed_order_sum']])
