#!/usr/bin/env python3
"""Literal strong-tree controls for the separate thickening integration."""
import itertools,json,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(x,*msg):
 if not x:raise RuntimeError(msg)
def tournament(n,mask):
 Q=[0]*n;k=0
 for i in range(n):
  for j in range(i+1,n):
   if mask>>k&1:Q[i]|=1<<j
   else:Q[j]|=1<<i
   k+=1
 return tuple(Q)
def carousel(n):return tuple(sum(1<<((i+d)%n) for d in range(1,(n+1)//2)) for i in range(n))
def tt(n):return tuple(sum(1<<j for j in range(i+1,n)) for i in range(n))
def substitution(Q,children):
 out=[];blocks=[]
 for C in children:blocks.append(list(range(len(out),len(out)+len(C))));out.extend([0]*len(C))
 for i,C in enumerate(children):
  for j,r in enumerate(C):
   out[blocks[i][j]]=sum(1<<blocks[i][k] for k in range(len(C)) if r>>k&1)
   for h,B in enumerate(blocks):
    if Q[i]>>h&1:out[blocks[i][j]]|=sum(1<<v for v in B)
 return tuple(out),blocks

def module_masks(Q):
 n=len(Q);full=(1<<n)-1;out=[]
 for mask in range(1,full+1):
  rem=full^mask;good=True
  while rem:
   b=rem&-rem;rem-=b;v=b.bit_length()-1;r=Q[v]&mask
   if r and r!=mask:good=False;break
  if good:out.append(mask)
 return out

def structure(Q):
 modules=module_masks(Q);full=(1<<len(Q))-1
 strong=[a for a in modules if not any((a&b) and (a&~b) and (b&~a) for b in modules)]
 children={a:[] for a in strong}
 for a in strong:
  if a==full:continue
  parent=min((b for b in strong if a&b==a and b!=a),key=int.bit_count)
  children[parent].append(a)
 prime=[];linear=[]
 for a,cs in children.items():
  cs.sort(key=lambda x:(x&-x).bit_length())
  if not cs:continue
  need(sum(cs)==a,'child partition')
  reps=[(x&-x).bit_length()-1 for x in cs]
  rows=tuple(sum(1<<j for j,v in enumerate(reps) if Q[u]>>v&1) for u in reps)
  islinear=sorted(r.bit_count() for r in rows)==list(range(len(rows)))
  if islinear:linear.append(dict(node=a,children=cs,rows=rows))
  else:
   need(set(module_masks(rows))=={*(1<<i for i in range(len(rows))), (1<<len(rows))-1},'nonprime nonlinear strong quotient')
   prime.append(dict(node=a,children=cs,rows=rows))
 return dict(modules=modules,strong=strong,prime=prime,linear=linear)

def check_case(Q,N,label):
 old=structure(Q);new,blocks=substitution(Q,[tt(N)]*len(Q));tree=structure(new)
 masks=[sum(1<<v for v in b) for b in blocks]
 def lift(mask):return sum(masks[i] for i in range(len(Q)) if mask>>i&1)
 expected={lift(r['node']):(tuple(lift(c) for c in r['children']),tuple(r['rows'])) for r in old['prime']}
 found={r['node']:(tuple(r['children']),tuple(r['rows'])) for r in tree['prime']}
 need(found==expected,'prime quotient/child preservation',Q,N,expected,found)
 return dict(label=label,old_rows=list(Q),inflation=N,physical_rows=list(new),physical_cells=blocks,old_modules=old['modules'],old_strong=old['strong'],old_prime_nodes=old['prime'],inflated_modules=tree['modules'],inflated_strong=tree['strong'],inflated_prime_nodes=tree['prime'],tested_subsets=(1<<len(Q))-1+(1<<len(new))-1)

def main():
 records=[]
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   Q=tournament(n,code)
   for N in (2,3):records.append(check_case(Q,N,f'all_order{n}_code{code}'))
 for name,Q in [('carousel5',carousel(5)),('carousel7',carousel(7)),('prime8',(238,4,240,118,195,146,34,74)),('carousel9',carousel(9))]:records.append(check_case(Q,2,name))
 nested,_=substitution(carousel(3),[tt(2),carousel(3),tt(1)]);records.append(check_case(nested,2,'nested_prime_and_linear'))
 data=json.dumps(records,sort_keys=True,separators=(',',':'))+'\n';(ROOT/'structural_controls.json').write_text(data)
 result=dict(status='PASS',scope='Literal verification of prime-quotient and physical-child preservation under transitive inflation. These small N=2,3 controls do not claim to meet the larger Seidel sufficient exponents or exhaust module-family targets.',literal_inflations=len(records),maximum_physical_order=max(len(r['physical_rows']) for r in records),nonempty_subsets_tested=sum(r['tested_subsets'] for r in records),prime_node_comparisons=sum(len(r['old_prime_nodes']) for r in records),all_small_base_tournaments=75,source_sha256=hashlib.sha256((ROOT/'PROOF.md').read_bytes()).hexdigest(),data_sha256=hashlib.sha256(data.encode()).hexdigest())
 print(json.dumps(result,sort_keys=True,indent=2))
if __name__=='__main__':main()
