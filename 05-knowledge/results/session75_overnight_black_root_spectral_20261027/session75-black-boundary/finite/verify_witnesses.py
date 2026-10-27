#!/usr/bin/env python3
"""Portable independent exhaustive certificate replay; no producer imports."""
from arcset_checks import arcs,opposite,connected,core,blue,replay,require,cut_parameters
import json,csv,gzip,pathlib,hashlib,collections
H=pathlib.Path(__file__).resolve().parent.parent

def source_rows(n,code):
 a=[0]*n;bit=n*(n-1)//2-1
 for i in range(n):
  for j in range(i+1,n):
   u,v=(i,j)if code>>bit&1 else(j,i);a[u]|=1<<v;bit-=1
 return tuple(a)
def hp(e,n):
 def go(p,left):
  if not left:yield tuple(p);return
  for v in sorted(left):
   if not p or(p[-1],v)in e:yield from go(p+[v],left-{v})
 yield from go([],set(range(n)))
def main():
 stats={};words=0;digest=hashlib.sha256()
 with gzip.open(H/'finite/witnesses.jsonl.gz','rt')as certs:
  for n in range(4,10):
   p=H/'inputs'/f'classes_n{n}.tsv'if n<9 else H/'inputs/n9_classes.tsv.gz';counts=collections.Counter();params=collections.Counter()
   with(gzip.open if n==9 else open)(p,'rt')as f:
    for ix,r in enumerate(csv.DictReader(f,delimiter='\t')):
     line=next(certs);digest.update(line.encode());z=json.loads(line);code=int(r['canon'if n<9 else'code']);require((z['n'],z['i'],z['c'])==(n,ix,code),'source index or code mismatch');e=arcs(source_rows(n,code));st=connected(e,range(n));counts['ordinary_classes']+=1
     if z['t']=='N':require(not st,'wrong nonstrong filter');counts['non_strong']+=1;continue
     require(st,'retained nonstrong source')
     if z['t']=='P':
      pp=list(hp(e,n));require(pp and all(blue(e,q)for q in pp),'false pureblue exclusion');require(n==5 and len(pp)==15,'unexpected pureblue exception');counts['pure_blue']+=1;continue
     require(z['t']=='W','unknown tag');delta,lam,_=cut_parameters(e,n);require(delta==lam,'lower-bound mismatch');initial=opposite(e)if z['d']else e;require(n-1-sum(u==z['v']for u,w in initial)==delta,'port not extremal')
     count=replay(e,n,{'route':z['p'],'dual':z['d'],'port':z['v']},'B',delta);words+=count;counts['strong_black_active']+=1;counts['literal_R_and_K_edges']+=count;params[str(delta)]+=1
   stats[str(n)]={**dict(counts),'exact_distance_counts':dict(sorted(params.items()))}
  require(next(certs,None)is None,'extra unbound certificate')
 # The only strong order3 class, explicitly excluded from the black-active domain.
 e=arcs((2,4,1));pp=list(hp(e,3));require(connected(e,range(3))and len(pp)==3 and all(blue(e,p)for p in pp),'C3 boundary failed')
 require([stats[str(n)]['ordinary_classes']for n in range(4,10)]==[4,12,56,456,6880,191536],'class counts changed')
 print(json.dumps({'status':'PASS','scope':'all ordinary strong black-active classes at orders4 through9','source_coverage':'every source index exactly once; completeness independently proved by authenticate_universe.cpp','pure_blue_exclusions':['C3 (3 Hamilton words)','R5 (15 Hamilton words)'],'certified_targets':['nonstrong distance=lambda','singleton plus strong core distance=delta'],'certificate_uncompressed_sha256':digest.hexdigest(),'literal_R_edges_checked':words,'literal_K_edges_checked':words,'orders':stats},sort_keys=True,indent=2))
if __name__=='__main__':main()
