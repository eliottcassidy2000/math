#!/usr/bin/env python3
import argparse,hashlib,itertools,json
from pathlib import Path
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser()
parser.add_argument('--source',type=Path,default=HERE.parent/'session195-eliminate-grammar-exceptions')
SRC=parser.parse_args().source.resolve()
def need(ok,msg):
 if not ok:raise RuntimeError(msg)
def graphcheck(a):
 n=len(a)
 for i,x in enumerate(a):
  need(0<=x<1<<n and not(x>>i&1),'row range/loop')
  for j in range(i):need(bool(x>>j&1)!=bool(a[j]>>i&1),'arc')
def data(a):
 graphcheck(a);n=len(a);N=1<<n;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)];H=[0]*N;H[0]=1
 for v in range(n):E[1<<v][v]=S[1<<v][v]=1
 for mask in range(1,N):
  ins=[i for i in range(n) if mask>>i&1];outs=[j for j in range(n) if not(mask>>j&1)]
  H[mask]=sum(E[mask])
  for i in ins:
   for j in outs:
    if a[i]>>j&1:E[mask|1<<j][j]+=E[mask][i]
    else:S[mask|1<<j][j]+=S[mask][i]
 full=N-1;result=[]
 for r in range(n):
  L=sum(E[m][r]*H[full^m] for m in range(N));R=sum(S[m][r]*H[full^m] for m in range(N));de=full^1<<r
  B=sum(H[m]*H[de^m] for m in range(N) if m&de==m)
  result.append((H[full],S[full][r],L,E[full][r],R,B,int(n==1),H[de]))
 return result

def literal(a,r):
 n=len(a);h=e=s=L=R=0
 def path(p):return all(a[i]>>j&1 for i,j in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  k=p.index(r)
  if path(p):h+=1;e+=p[-1]==r;s+=p[0]==r
  if path(p[:k+1]) and path(p[k+1:]):L+=1
  if path(p[:k]) and path(p[k:]):R+=1
 B=d=0
 for p in itertools.permutations(i for i in range(n) if i!=r):
  d+=path(p)
  B+=sum(path(p[:k]) and path(p[k:]) for k in range(n))
 return h,s,L,e,R,B,int(n==1),d

def state(p):return p[0],p[3],p[4],p[7]
def signature(g):return tuple(g['signature'])+(g['H_deleted_root'],)
def seam(a,r,b,t):
 n=len(a);m=len(b);rows=[x|(((1<<m)-1)<<n) for x in a]+[x<<n for x in b]
 rows[r]&=~(1<<(n+t));rows[n+t]|=1<<r
 return rows,n+t

def strong(a):
 R=[x|1<<i for i,x in enumerate(a)]
 for k in range(len(a)):
  for i in range(len(a)):
   if R[i]>>k&1:R[i]|=R[k]
 return all(x==(1<<len(a))-1 for x in R)
def induced(a,V):return [sum(1<<j for j,y in enumerate(V) if a[x]>>y&1) for x in V]
def graphs(n):
 pairs=list(itertools.combinations(range(n),2))
 for v in range(1<<len(pairs)):
  a=[0]*n
  for k,(i,j) in enumerate(pairs):
   if v>>k&1:a[i]|=1<<j
   else:a[j]|=1<<i
  yield a

def main():
 names=['PROOF.md','verify.py','census.cpp','enumerate.cpp','base_inventory.json','gadgets.json','inventory6.json','bank67.txt','witnesses.json','single_addition_audit.json','receipt-full-inventory-1000.json','census-20000.json']
 hashes={name:hashlib.sha256((SRC/name).read_bytes()).hexdigest() for name in names}
 base=json.loads((SRC/'base_inventory.json').read_text());six=json.loads((SRC/'inventory6.json').read_text());new=json.loads((SRC/'gadgets.json').read_text());wits=json.loads((SRC/'witnesses.json').read_text())
 need(len(base)==65 and len(six)==280 and len(new)==2,'inventory sizes')
 need([signature(g) for g in base+new]==[tuple(map(int,line.split())) for line in (SRC/'bank67.txt').read_text().splitlines()],'raw bank mismatch')
 physical=[]
 for g in new:
  a=g['rows'];r=g['root'];actual=data(a)[r]
  need(actual==literal(a,r)==signature(g),'gadget signature')
  physical.append({'rows':a,'root':r,'signature':actual})
 for ix,w in enumerate(wits):
  g=new[ix];a,r=([0],0) if ix==0 else ([2,0],0)
  rows,b=seam(a,r,g['rows'],g['root'])
  need(rows==w['rows'] and b==w['root'],'witness graph reconstruction')
  p=data(rows)[b];need(p==literal(rows,b),'witness literal')
  need(state(p)==tuple(w['state']) and p[0]==w['target'] and strong(rows),'witness state/strong')
 need(data(induced(new[0]['rows'],[2,3,4,5]))[0][0]==1,'A feedback deletion')
 need(data(induced(new[1]['rows'],[1,2,3,4,5]))[0][0]==1,'B feedback deletion')
 checks=0
 for n in range(1,5):
  for a in graphs(n):
   P=data(a)
   for r,p in enumerate(P):
    h,e,R,d=state(p)
    for z in new:
     g,s,L,E,Q,B,delta,dg=signature(z);pred=(g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,dg*h)
     rows,b=seam(a,r,z['rows'],z['root']);need(state(data(rows)[b])==pred,'all-small transfer');checks+=1
 independent_base={tuple(map(int,l.split())) for l in (HERE/'independent-base.txt').read_text().splitlines()}
 independent_six=sorted(tuple(map(int,l.split())) for l in (HERE/'independent-six.txt').read_text().splitlines())
 need(independent_base==set(map(signature,base)),'full base inventory')
 need(set(independent_six)==set(map(signature,six)),'full six inventory')
 independent_single={independent_six[j]:(count,bool(a),bool(b)) for j,count,a,b in [map(int,l.split()) for l in (HERE/'independent-single.txt').read_text().splitlines()]}
 for row in json.loads((SRC/'single_addition_audit.json').read_text()):
  key=signature(six[row['inventory_index']]);count,a,b=independent_single[key]
  need(count==row['retained_states'] and a==(67 in row['attained_targets']) and b==(173 in row['attained_targets']),'single candidate mismatch')
 stdout=(HERE/'independent.stdout').read_text().splitlines()
 need('SINGLE 280 6 5' in stdout,'single aggregate')
 need('EXTENDED 1000 21353 7 21' in stdout,'1000 census')
 need('EXTENDED 20000 2584754 7 21' in stdout,'20000 census')
 need(all(hashlib.sha256((SRC/name).read_bytes()).hexdigest()==sha for name,sha in hashes.items()),'source changed during audit')
 receipt={'status':'PASS','source_sha256':hashes,'physical_gadgets':physical,'physical_target_witnesses':[67,173],'independent_small_input_graph_transfers':checks,'base_full_inventory':65,'order_six_full_inventory':280,'labelled_order_six_graphs':32768,'rooted_order_six_gadgets':196608,'single_candidate_runs':280,'single_candidates_covering_both':0,'single_candidates_covering67':6,'single_candidates_covering173':5,'census1000_states':21353,'census20000_states':2584754,'missing_at_both_bounds':[7,21],'scope':'Exact stated finite grammar and target bounds; no global spectrum or unbounded coverage claim.'}
 (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,sort_keys=True)+'\n');print(json.dumps(receipt,sort_keys=True))
if __name__=='__main__':main()
