"""Exact equinumerosity/trace census and marked two-port Fourier bridge at n=5."""
import argparse
from collections import Counter
import gzip
from hashlib import sha256
import importlib.util
from itertools import combinations,permutations
import json
from pathlib import Path
import subprocess
import sys
import tempfile

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'homogeneous_pair_response_20261103'
sys.path.insert(0,str(OLD))
import check
import response as poly
need=check.need
QUARTIC=HERE.parent/'session188_structural_bridges_20261104/response/verify.py'
spec=importlib.util.spec_from_file_location('cover_form',QUARTIC)
cover=importlib.util.module_from_spec(spec);spec.loader.exec_module(cover)


def action(n):
 pairs=list(combinations(range(n),2));index={e:i for i,e in enumerate(pairs)};out=[]
 for p in permutations(range(n)):
  targets=[index[tuple(sorted((p[u],p[v])))] for u,v in pairs]
  inv=sum(1<<i for i,(u,v) in enumerate(pairs) if p[u]>p[v])
  def image(x):return sum(1<<targets[i] for i in range(len(pairs)) if x>>i&1)
  out.append((inv,image(inv),[image(x) for x in range(1<<len(pairs))]))
 return out


def census(n):
 acts=action(n);size=1<<(n*(n-1)//2);tc={};gc={};tours=[];graphs=[]
 for x in range(size):
  if x not in tc:
   orb={images[x]^target for inv,target,images in acts}
   for y in orb:tc[y]=x
   tours.append(dict(code=x,orbit=len(orb),aut=len(acts)//len(orb),H=check.h(check.tournament(n,x))))
  if x not in gc:
   orb={images[x] for inv,target,images in acts}
   for y in orb:gc[y]=x
   even=all((x&inv).bit_count()%2==0 for inv,target,images in acts if images[x]==x)
   graphs.append(dict(code=x,orbit=len(orb),aut=len(acts)//len(orb),sign_even=even,edges=x.bit_count()))
 for t in tours:t['converse']=tc[(size-1)^t['code']]
 good=[g for g in graphs if g['sign_even']]
 sc=sum(t['converse']==t['code'] for t in tours);ep=sum(g['edges']%2==0 for g in good);op=len(good)-ep
 need(len(good)==len(tours),'same-order type equinumerosity')
 need(sc==ep-op,'converse signed trace')
 return dict(n=n,tournament_types=tours,graph_types=graphs,self_converse=sc,
             sign_even_even_edges=ep,sign_even_odd_edges=op,labelled_sign_even=sum(g['orbit'] for g in good))


def hadamard(rows):
 out=[r[:] for r in rows];step=1
 while step<len(out):
  for b in range(0,len(out),2*step):
   for j in range(step):
    x,y=out[b+j],out[b+j+step]
    out[b+j]=[u+v for u,v in zip(x,y)];out[b+j+step]=[u-v for u,v in zip(x,y)]
  step*=2
 return out


def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path,default=HERE/'receipt.json')
 p.add_argument('--tables',type=Path,default=HERE/'full_response.json.gz');args=p.parse_args()
 small=[census(n) for n in range(6)]
 need([(x['self_converse'],x['sign_even_even_edges'],x['sign_even_odd_edges']) for x in small[3:]]==[(2,2,0),(2,3,1),(8,10,2)],'small trace controls')
 n=5;ms=poly.masks(n);battery=[(a,c) for a in ms for c in ms];high=[(31,31),(31,0),(0,31)];qs=battery+high
 rows=[check.tournament(n,x) for x in range(1024)]
 inp=['1024']
 for a in rows:inp += [f'5 {len(qs)}',' '.join(map(str,a))]+[f'{x} {y}' for x,y in qs]
 with tempfile.TemporaryDirectory(prefix='even-response-') as folder:
  exe=Path(folder)/'counter.exe'
  subprocess.run(['g++','-O3','-std=c++17',str(OLD/'response.cpp'),'-o',str(exe)],check=True)
  run=subprocess.run([str(exe)],input='\n'.join(inp)+'\n',text=True,capture_output=True,check=True)
 vals=list(map(int,run.stdout.split()));need(len(vals)==1024*len(qs),'actual query stream')
 values=[vals[i*len(qs):(i+1)*len(qs)] for i in range(1024)]
 co=[poly.coefficients(v[:256],ms) for v in values]
 quartics=0
 for x in range(1024):
  form=cover.quartic_normal_form(rows[x])
  for j,a in enumerate(ms):
   for k,c in enumerate(ms):
    if a.bit_count()==c.bit_count()==2:
     need(co[x][j*len(ms)+k]==form[a|(c<<n)],'independent full named quartic');quartics+=1
  for q,v in zip(high,values[x][256:]):need(poly.evaluate(co[x],ms,*q)==v,'dense full-response reconstruction')
 # Raw transform uses (-1)^popcount(code & support); correct to the upward +1 gauge.
 numerator=hadamard(co)
 numerator=[[-z if support.bit_count()%2 else z for z in row] for support,row in enumerate(numerator)]
 inverses=hadamard([[-z if support.bit_count()%2 else z for z in row] for support,row in enumerate(numerator)])
 need(all(inverses[i]==[1024*z for z in co[i]] for i in range(1024)),'full old-edge Fourier inverse')
 nonzero=Counter();quartic_odd=0
 for support,vec in enumerate(numerator):
  for k,z in enumerate(vec):
   degree=ms[k//len(ms)].bit_count()+ms[k%len(ms)].bit_count()
   if z:nonzero[degree,support.bit_count()]+=1
   if degree==4 and support.bit_count()%2:need(z==0,'quartics are converse-even');quartic_odd+=1
 for x in range(1024):
  need(all(co[x][k]==co[1023^x][k] for k in range(256) if ms[k//16].bit_count()+ms[k%16].bit_count()==4),'whole-converse cubic')
 # A fixed named query does not descend to unmarked tournament types.
 t=[6,4,0];u=check.opposite(t);query=(1,0)
 hostile=[check.literal_h(check.extend(a,*query)) for a in (t,u)]
 need(hostile==[3,9],'isomorphic named-query hostile')
 # Literal edge-support expansion at n=3: a globally sign-odd one-edge support survives marking.
 f=[check.h(check.extend(check.tournament(3,x),*query)) for x in range(8)]
 hat=hadamard([[z] for z in f]);hat=[((-1)**s.bit_count())*v[0] for s,v in enumerate(hat)]
 marked=next((s,z) for s,z in enumerate(hat) if s.bit_count()==1 and z)
 # Independent full-word polynomial at twelve actual n=5 class representatives.
 literal=0
 for t in small[-1]['tournament_types']:
  x=t['code'];direct=cover.polynomial(rows[x])
  need(all(co[x][i*16+j]==direct.get(a|(c<<5),0) for i,a in enumerate(ms) for j,c in enumerate(ms)),'literal full named polynomial')
  literal+=1
 data=dict(n=5,old_code='bit1 is upward arc; pairs lexicographic',masks=ms,
           attachment_basis='row-major (a,c), full bidegree (2,2), zeros retained',
           values=values,queries=qs,coefficients=co,old_edge_walsh_numerators=numerator,denominator=1024)
 raw=(json.dumps(data,sort_keys=True,separators=(',',':'))+'\n').encode();packed=gzip.compress(raw,mtime=0)
 args.tables.write_bytes(packed)
 sources=[Path(__file__).resolve(),OLD/'check.py',OLD/'response.py',OLD/'response.cpp',QUARTIC]
 result=dict(status='PASS complete n<=5 type census; complete n=5 named two-port Fourier transform',
             census=small,actual_queries=len(vals),full_named_responses=1024,coefficients_per_response=256,
             quartic_normal_form_checks=quartics,odd_old_support_quartic_zeros=quartic_odd,
             literal_full_polynomials=literal,nonzero_bidegree_old_support={f'{a},{b}':v for (a,b),v in sorted(nonzero.items())},
             named_isomorphic_hostile=dict(rows=[[6,4,0],[0,1,3]],a=1,c=0,counts=hostile),
             relative_even_one_edge=dict(support=marked[0],walsh_numerator=marked[1],denominator=8),
             sources={str(s.relative_to(HERE.parent)).replace('\\','/'):sha256(s.read_bytes()).hexdigest() for s in sources},
             tables_sha256=sha256(packed).hexdigest(),uncompressed_sha256=sha256(raw).hexdigest())
 args.output.write_text(json.dumps(result,sort_keys=True,indent=2)+'\n',encoding='utf8',newline='\n')
 print(json.dumps({k:v for k,v in result.items() if k not in ('census','sources','nonzero_bidegree_old_support')},indent=2))


if __name__=='__main__':main()
