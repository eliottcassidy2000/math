#!/usr/bin/env python3
"""Independent set-of-arcs certificate replay. No producer imports. -O safe."""
import json,itertools,hashlib,sys,os,csv,gzip,math,collections
H=os.path.dirname(os.path.abspath(__file__))
def require(q,m):
 if not q:raise RuntimeError(m)
def arcs(rows):return frozenset((u,v)for u,r in enumerate(rows)for v in range(len(rows))if r>>v&1)
def opposite(e):return frozenset((v,u)for u,v in e)
def connected(e,vs):
 vs=set(vs)
 if not vs:return False
 for E in (e,opposite(e)):
  got={min(vs)}
  while True:
   new=got|{v for u,v in E if u in got and v in vs}
   if new==got:break
   got=new
  if got!=vs:return False
 return True
def core(e,n):
 return any(all((v,u)in e for u in range(n)if u!=v)and connected(e,set(range(n))-{v})or all((u,v)in e for u in range(n)if u!=v)and connected(e,set(range(n))-{v})for v in range(n))
def blue(e,p):
 image={p[i]:p[-1-i]for i in range(len(p))}
 return all((image[v],image[u])in e for u,v in e)
def reverse(e,p):
 flips=set(zip(p,p[1:]));return frozenset((e-flips)|{(v,u)for u,v in flips})
def cut_parameters(e,n):
 ds=[sum(u==v for u,w in e)for v in range(n)];delta=min(ds+[n-1-d for d in ds]);full=(1<<n)-1;sums=[0]*(full+1);minimum=n*n;shores=[]
 for s in range(1,full):
  b=s&-s;v=b.bit_length()-1;sums[s]=sums[s^b]+ds[v];k=s.bit_count();c=sums[s]-k*(k-1)//2
  if c<minimum:minimum=c;shores=[s]
  elif c==minimum:shores.append(s)
 return delta,minimum,shores

def replay(e,n,record,lane,depth):
 pps=record['route'];require(len(pps)==depth,'wrong length');raw=e
 if lane=='B'and record['dual']:e=opposite(e)
 initial=e;states=[];failures=[]
 for j,p in enumerate(pps):
  require(sorted(p)==list(range(n)),'not named permutation');require(all((u,v)in e for u,v in zip(p,p[1:])),'not directed');require(not blue(e,p),'blue producer witness')
  if lane=='A':
   s=record['shore'];require(s>>p[0]&1 and not(s>>p[-1]&1),'lost minimum-shore port')
  else:require(p[-1]==record['port'],'lost fixed endpoint')
  e2=reverse(e,p);require(len(e2)==n*(n-1)//2,'arc count changed');require(not blue(e2,p[::-1]),'blackness failed across R')
  # Literal K edge: alternate converse, independent nonpath complement.
  s=e if j%2==0 else opposite(e);q=p if j%2==0 else p[::-1];pathpairs={frozenset((u,v))for u,v in zip(q,q[1:])}
  k=frozenset((u,v)for u in range(n)for v in range(n)if u!=v and (((u,v)in s)==(frozenset((u,v))in pathpairs)))
  expected=opposite(e2)if j%2==0 else e2
  require(k==expected,'literal K mismatch');require(not blue(s,q)and not blue(k,q),'literal K blackness')
  require(connected(e,range(n)),'nonstrong intermediate on geodesic witness')
  states.append(connected(e2,range(n)));e=e2
 require(not connected(e,range(n)),'A terminal still strong')
 if lane=='B':require(core(e,n),'not singleton and strong core')
 return len(pps)

def audit_file(name,expected=None,source=None):
 path=f'{H}/{name}.jsonl';report={'records':0,'strong':0,'non_strong':0,'A_edges':0,'B_edges':0,'parameter_counts':collections.Counter(),'failed_A_cones':0,'failed_B_cones':0};digest=hashlib.sha256();seen=set()
 stream=None
 if source:
  p=f'{H}/inputs/{source}';stream=(gzip.open if p.endswith('.gz')else open)(p,'rt');table=iter(csv.DictReader(stream,delimiter='\t'))
 for line in open(path):
  r=json.loads(line);digest.update(line.encode());a=r['rows'];n=len(a);e=arcs(a);require(len(e)==n*(n-1)//2,'not tournament');report['records']+=1
  if source:
   row=next(table);code=int(row['canon'if n==8 else 'code']);require(code==r['code'],'source table mismatch');bits=0
   for u in range(n):
    for v in range(u+1,n):bits=(bits<<1)+int((u,v)in e)
   require(code==bits,'encoding mismatch');require(code not in seen,'duplicate code');seen.add(code)
  st=connected(e,range(n))
  if not st:require(r['verdict']=='NOT_STRONG','strong filter mismatch');report['non_strong']+=1;continue
  report['strong']+=1;delta,lam,shores=cut_parameters(e,n);require((delta,lam,shores)==(r['delta'],r['lambda'],r['minimum_shores']),'parameter mismatch');report['parameter_counts'][f'{lam},{delta}']+=1
  for lane,dep in [('A',lam),('B',delta)]:
   witness=next((c for c in r[lane]if c['verdict']=='WITNESS'),None);require(witness is not None,'missing global witness')
   report[f'{lane}_edges']+=replay(e,n,witness,lane,dep);report[f'failed_{lane}_cones']+=sum(c['verdict']=='EXCLUDED'for c in r[lane])
 if expected is not None:require(report['records']==expected,'incomplete run')
 if source:
  require(next(table,None)is None,'table not exhausted');stream.close()
 report['receipt_sha256']=digest.hexdigest();report['parameter_counts']=dict(sorted(report['parameter_counts'].items()));return report

def main():
 names=sys.argv[1:]or['probe','selfconverse','classes8','classes9'];result={}
 for name in names:
  cfg={'probe':(1268,None),'selfconverse':(5000,None),'classes8':(6880,'classes_n8.tsv'),'classes9':(191536,'n9_classes.tsv.gz')}[name]
  result[name]=audit_file(name,*cfg)
 print(json.dumps({'status':'PASS','independent_representation':'sets of directed arcs; literal K nonpath complement','source_universe_completeness':'inherited authenticated class tables, exact record correspondence checked','results':result},sort_keys=True,indent=2))
if __name__=='__main__':main()
