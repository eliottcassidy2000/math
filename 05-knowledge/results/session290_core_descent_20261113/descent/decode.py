"""Exact target recognition and physical construction in the inherited 70 charts.
This is NOT a universal odd-count selector. See PROOF.md.
"""
from pathlib import Path
import gzip,json,argparse
HERE=Path(__file__).resolve().parent

def need(ok,why):
 if not ok:raise RuntimeError(why)
def load():return json.loads(gzip.decompress((HERE/'profiles.json.gz').read_bytes()))
def profile(r,kind):return r['U' if kind=='even' else 'child_U']
def val(E,m):return sum(c*b**m for b,c in E)
def expand(rows,selected,ports,m):
 labels=[v for v in range(len(rows)) for _ in range(m if v in ports else 1)]
 out=[sum(1<<j for j,w in enumerate(labels) if (v==w and i<j) or (v!=w and rows[v]>>w&1)) for i,v in enumerate(labels)]
 mask=sum(1<<i for i,v in enumerate(labels) if selected>>v&1)
 return out,mask

def tail(rows,mask,t):
 need(t>=1,'positive physical tail length')
 n=len(rows);root=n+t;out=rows.copy()+[0]*(t+1)
 for v in range(n):
  for w in range(n,root):out[v]|=1<<w
  if mask>>v&1:out[root]|=1<<v
  else:out[v]|=1<<root
 for v in range(n,root):
  for w in range(v+1,root+1):out[v]|=1<<w
 return out

def construct(data,kind,index,m,t):
 need(m>=1 and t>=1,'positive clocks')
 r=data['even_records' if kind=='even' else 'records'][index]
 rows,mask=expand(r['rows'],r['selected'],r['ports'],m)
 u=val(profile(r,kind),m);rr=tail(rows,mask,t)
 return {'kind':kind,'family':index,'m':m,'t':t,'U':u,'H':1+(1<<t)*u,'rows':rr,'source_rows':rows,'selected':mask,'order':len(rr)}

def recognize(N,data=None):
 data=load() if data is None else data
 if N<=1 or N%2==0:return []
 v=((N-1)&-(N-1)).bit_length()-1;hits=[]
 for t in range(1,v+1):
  u=(N-1)>>t
  for kind,rs in [('even',data['even_records']),('odd',data['records'])]:
   for i,r in enumerate(rs):
    E=profile(r,kind);lo=1;hi=u.bit_length()-1
    while lo<=hi:
     m=(lo+hi)//2;x=val(E,m)
     if x==u:hits.append({'kind':kind,'family':i,'m':m,'t':t});break
     if x<u:lo=m+1
     else:hi=m-1
 return hits

def descent(hit,data=None):
 data=load() if data is None else data
 kind=hit['kind'];i=hit['family'];m=hit['m'];t=hit['t'];trace=[]
 while True:
  r=data['even_records' if kind=='even' else 'records'][i];u=val(profile(r,kind),m)
  trace.append({'kind':kind,'family':i,'m':m,'t':t,'U':u,'H':1+(1<<t)*u})
  if kind=='odd':i=r['parent_synchrony_family'];kind='even';m+=1
  elif m>1:m-=1
  elif t>1:t-=1
  else:break
 need(all(a['H']>b['H'] for a,b in zip(trace,trace[1:])),'strict actual target descent')
 return trace

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('target',type=int);p.add_argument('--rows',action='store_true');a=p.parse_args();d=load();hits=recognize(a.target,d)
 out={'target':a.target,'scope':'membership only in the70 inherited charts with the free port fixed at1; not arbitrary spectrum membership','representations':hits}
 if hits:
  out['descent']=descent(hits[0],d)
  if a.rows:out['construction']=construct(d,**{'kind':hits[0]['kind'],'index':hits[0]['family'],'m':hits[0]['m'],'t':hits[0]['t']})
 print(json.dumps(out,indent=2))
