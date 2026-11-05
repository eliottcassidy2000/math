#!/usr/bin/env python3
"""Independent exact replay: two physical six-vertex additions to rooted grammar.
No imports from other packets; explicit checks survive python -O.
"""
import argparse, hashlib, itertools, json, time
from collections import deque
from pathlib import Path
HERE=Path(__file__).resolve().parent

def need(ok,msg):
 if not ok:raise RuntimeError(msg)

def profile(a):
 n=len(a);N=1<<n;H=[0]*N;H[0]=1;E=[[0]*n for _ in range(N)];S=[[0]*n for _ in range(N)]
 for i in range(n):
  need(a[i]>>i&1==0,'loop')
  for j in range(i):need(((a[i]>>j)^(a[j]>>i))&1,'invalid tournament')
 for m in range(1,N):
  vertices=[v for v in range(n) if m>>v&1]
  for v in vertices:
   q=m^(1<<v)
   if not q:E[m][v]=S[m][v]=1
   else:
    E[m][v]=sum(E[q][u] for u in vertices if u!=v and a[u]>>v&1)
    S[m][v]=sum(S[q][u] for u in vertices if u!=v and a[v]>>u&1)
  H[m]=sum(E[m]);need(H[m]==sum(S[m]),'two DP directions disagree')
 full=N-1;out=[]
 for b in range(n):
  L=sum(E[u][b]*H[full^u] for u in range(N) if u>>b&1)
  Q=sum(S[u][b]*H[full^u] for u in range(N) if u>>b&1)
  de=full^(1<<b);B=sum(H[u]*H[de^u] for u in range(N) if u&de==u)
  out.append((H[full],S[full][b],L,E[full][b],Q,B,int(n==1),H[de]))
 return out

def literal(a):
 n=len(a);h=0;e=[0]*n;s=[0]*n;L=[0]*n;R=[0]*n
 def path(p):return all(a[x]>>y&1 for x,y in zip(p,p[1:]))
 for p in itertools.permutations(range(n)):
  if path(p):h+=1;s[p[0]]+=1;e[p[-1]]+=1
  for k in range(n+1):
   if path(p[:k]) and path(p[k:]):
    if k:L[p[k-1]]+=1
    if k<n:R[p[k]]+=1
 out=[]
 for b in range(n):
  d=B=0
  for p in itertools.permutations(v for v in range(n) if v!=b):
   if path(p):d+=1
   for k in range(n):B+=path(p[:k]) and path(p[k:])
  out.append((h,s[b],L[b],e[b],R[b],B,int(n==1),d))
 return out

def state(p):return p[0],p[3],p[4],p[7]
def key(g):return tuple(g['signature']+[g['H_deleted_root']])
def transfer(x,p):
 h,e,R,d=x;g,s,L,E,Q,B,delta,dg=p
 return g*h-s*e+L*R,E*h-delta*e,Q*h+B*R,dg*h

def seam(a,r,g,b):
 n=len(a);m=len(g);new=((1<<m)-1)<<n
 out=[v|new for v in a]+[v<<n for v in g]
 out[r]^=1<<(n+b);out[n+b]|=1<<r
 return tuple(out),n+b

def closure(bank,X):
 pbank=[key(g) for g in bank];start=(1,1,1,1);q=deque([start]);seen={start};vals={1};hist={};wit={1:start};terminal={}
 while q:
  x=q.popleft();h,e,R,d=x
  successors=[(('old',None),(h+R-e,d,2*R-d,d)),(('sink',None),(h,0,2*R,d))]
  successors += [(('bank',i),transfer(x,p)) for i,p in enumerate(pbank)]
  for op,z in successors:
   H,E,Q,D=z
   need(H>=h and 0<=E<=H<=Q and 1<=D<=H,'physical bounds failed')
   if H>X:continue
   if H not in vals:wit[H]=z;terminal[z]=(x,op)
   vals.add(H)
   if Q<=X and z not in seen:seen.add(z);q.append(z);hist[z]=(x,op)
 return seen,vals,hist,wit,terminal

def tournaments(n):
 for mask in range(1<<(n*(n-1)//2)):
  a=[0]*n;k=0
  for i in range(n):
   for j in range(i+1,n):
    if mask>>k&1:a[i]|=1<<j
    else:a[j]|=1<<i
    k+=1
  yield tuple(a)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--bound',type=int,default=1000);ap.add_argument('--full-inventory',action='store_true');args=ap.parse_args()
 base=json.loads((HERE/'base_inventory.json').read_text());new=json.loads((HERE/'gadgets.json').read_text());six=json.loads((HERE/'inventory6.json').read_text())
 need(len(base)==65 and len(new)==2 and len(six)==280,'bank size')
 literal_checks=0
 for g in base+six:
  a=tuple(g['rows']);p=profile(a);l=literal(a)
  need(p==l,'independent profile mismatch');need(p[g['root']]==key(g),'inventory representative mismatch');literal_checks+=len(a)
 full_result=None
 if args.full_inventory:
  computed=set()
  for a in tournaments(6):computed.update(profile(a))
  need(computed==set(map(key,six)),'exhaustive six-vertex inventory mismatch')
  full_result={'labelled_graphs':32768,'rooted_gadgets':196608,'distinct_full_signatures':len(computed)}
 graphs=[]
 for target,g in zip([67,173],new):
  a=(0,);r=0
  if target==173:a=(2,0)
  a,r=seam(a,r,tuple(g['rows']),g['root']);p=profile(a)
  need(p==literal(a),'output graph literal/DP mismatch')
  want=(67,4,134,9) if target==67 else (173,8,192,9)
  need(state(p[r])==want,'output state mismatch')
  graphs.append({'target':target,'rows':list(a),'root':r,'state':list(want),'gadget':g})
 # Both all-order laws checked on every <=3 input and every root.
 checks=0
 for n in range(1,4):
  for a in tournaments(n):
   ps=profile(a)
   for r,p in enumerate(ps):
    for g in new:
     out,b=seam(a,r,tuple(g['rows']),g['root'])
     need(state(profile(out)[b])==transfer(state(p),key(g)),'actual composition mismatch');checks+=1
 # Exact arbitrary-depth count-bound closure independently confirms smaller bank omission.
 oldstates,oldvalues,*_=closure(base,173)
 need(67 not in oldvalues and 173 not in oldvalues,'old-bank omission differs')
 single=[]
 for i,g in enumerate(six):
  seen,vals,*_=closure(base+[g],173)
  have=[v for v in [67,173] if v in vals]
  need(len(have)<2,'unexpected single-gadget improvement')
  single.append({'inventory_index':i,'attained_targets':have,'retained_states':len(seen)})
 seen,vals,*_=closure(base+new,args.bound)
 missing=sorted(set(range(1,args.bound+1,2))-vals)
 need(missing==[v for v in [7,21] if v<=args.bound], 'expanded-bank census mismatch')
 receipt={'status':'PASS','bound':args.bound,'base_matrices':65,'additional_matrices':2,'retained_states':len(seen),'missing_odds':missing,'literal_root_profile_checks':literal_checks,'actual_output_literal_checks':2,'all_small_input_composition_checks':checks,'exhaustive_inventory':full_result,'single_six_vertex_addition_audit':{'candidates':len(six),'covering_both_targets':0,'covering_67':sum(67 in t['attained_targets'] for t in single),'covering_173':sum(173 in t['attained_targets'] for t in single)},'scope':'Exact arbitrary-depth closure of this fixed grammar only; no infinite coverage theorem.'}
 (HERE/'witnesses.json').write_text(json.dumps(graphs,indent=2)+'\n');(HERE/'single_addition_audit.json').write_text(json.dumps(single,indent=2)+'\n');(HERE/f'receipt-{args.bound}.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
