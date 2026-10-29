"""Exact weighted fixed sets, common transport kernels, and run-refinement controls."""
from pathlib import Path
from fractions import Fraction as Q
from itertools import product,permutations
from collections import Counter
import json,sys,random
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'moving_selector_atoms_20261028'))
import frames
core=frames.core;oracle=frames.r.oracle
CHECKS=0
def need(x,m):
 global CHECKS
 CHECKS+=1
 if not x:raise RuntimeError(m)
def neg(v,k):
 out=0
 for i in range(k):out+=((-v%3)%3)*3**i;v//=3
 return out
def aut(k):
 if k==0:return [(0,)]
 small=aut(k-1);m=3**(k-1)
 return [tuple(((i+a)%3)*m+gs[i][j] for i in range(3) for j in range(m))
         for a in range(3) for gs in product(small,repeat=3)]
def central(k):
 if k==0:return [(0,)]
 m=3**(k-1);out=[]
 for g in central(k-1):
  for h in aut(k-1):out.append(tuple(g)+tuple(m+h[j] for j in range(m))+tuple(2*m+neg(h[neg(j,k-1)],k-1) for j in range(m)))
 return out
def fixed_poly(k,z,kind):
 if k==0:return Q(z[0])
 m=3**(k-1)
 if kind=='A':return Q(2,3)+fixed_poly(k-1,z[:m],'A')*fixed_poly(k-1,z[m:2*m],'A')*fixed_poly(k-1,z[2*m:],'A')/3
 paired=[z[m+j]*z[2*m+neg(j,k-1)] for j in range(m)]
 return fixed_poly(k-1,z[:m],'C')*fixed_poly(k-1,paired,'A')
def kernel(k,u,v,kind):
 if len(u)!=len(v):return Q(0)
 if k==0:return Q(1)
 m=3**(k-1)
 if kind=='C':
  if any(a//m!=b//m for a,b in zip(u,v)):return Q(0)
  zero=[i for i,x in enumerate(u) if x//m==0];other=[i for i,x in enumerate(u) if x//m!=0]
  fold=lambda x:x%m if x//m==1 else neg(x%m,k-1)
  return kernel(k-1,tuple(u[i]%m for i in zero),tuple(v[i]%m for i in zero),'C')*kernel(k-1,tuple(fold(u[i]) for i in other),tuple(fold(v[i]) for i in other),'A')
 total=Q(0)
 for shift in range(3):
  if any((a//m+shift)%3!=b//m for a,b in zip(u,v)):continue
  value=Q(1)
  for child in range(3):
   ids=[i for i,x in enumerate(u) if x//m==child]
   value*=kernel(k-1,tuple(u[i]%m for i in ids),tuple(v[i]%m for i in ids),'A')
  total+=value/3
 return total
def edge(v,w,k):
 if v==w:return False
 m=3**(k-1)
 return (w//m-v//m)%3==1 if v//m!=w//m else edge(v%m,w%m,k-1)
def tournament(k):
 n=3**k
 return tuple(sum(1<<j for j in range(n) if edge(i,j,k)) for i in range(n))
def permute(T,g):
 out=[0]*len(T)
 for i in range(len(T)):
  for j in range(len(T)):
   if T[i]>>j&1:out[g[i]]|=1<<g[j]
 return tuple(out)
def paths(T):
 n=len(T);full=(1<<n)-1
 def rec(w,mask):
  if mask==full:yield tuple(w);return
  options=(T[w[-1]] if w else full)&~mask
  while options:
   b=options&-options;options-=b;j=b.bit_length()-1
   yield from rec(w+[j],mask|b)
 yield from rec([],0)
def controls():
 report=[];rng=random.Random(11212)
 for kind,k in [('A',1),('A',2),('C',1),('C',2),('C',3)]:
  G=aut(k) if kind=='A' else central(k);n=3**k
  for trial in range(6):
   z=[rng.randrange(-2,4) for _ in range(n)];s=0
   for g in G:
    p=1
    for v in range(n):
     if g[v]==v:p*=z[v]
    s+=p
   need(fixed_poly(k,z,kind)==Q(s,len(G)),'multivariate fixed-set recurrence')
  sources=list(product(range(n),repeat=2)) if k<=2 else [tuple(rng.randrange(n) for _ in range(3)) for _ in range(16)]
  for u in sources:
   counts=Counter(tuple(g[x] for x in u) for g in G)
   targets=list(product(range(n),repeat=2)) if k<=2 else list(counts)+[tuple(rng.randrange(n) for _ in u) for _ in range(16)]
   for v in targets:need(kernel(k,u,v,kind)==Q(counts[v],len(G)),'common transporter kernel')
  fixed=[sum(g[v]==v for v in range(n)) for g in G]
  mean=Q(sum(fixed),len(G));second=Q(sum(v*v for v in fixed),len(G))
  if kind=='C':need(mean==2*k+1 and second==8*k*k+1,'inherited port moments')
  report.append({'group':kind,'depth':k,'order':len(G),'mean_fixed':str(mean),'second_moment':str(second)})
 # The same entire weighted fixed-set polynomial can hide a different joint action.
 plus=[tuple(range(3))+tuple(3+(j+a)%3 for j in range(3))+tuple(6+(j+a)%3 for j in range(3)) for a in range(3)]
 minus=central(2);T=tournament(2);t=tuple(neg(v,2) for v in range(9))
 def masks(G):return Counter(tuple(v for v in range(9) if g[v]==v) for g in G)
 def cycle_types(G):
  out=Counter()
  for g in G:
   seen=set();lengths=[]
   for v in range(len(g)):
    if v in seen:continue
    w=v;length=0
    while w not in seen:seen.add(w);length+=1;w=g[w]
    lengths.append(length)
   out[tuple(sorted(lengths))]+=1
  return out
 need(masks(plus)==masks(minus),'all weighted fixed-set observations coincide')
 need(cycle_types(plus)==cycle_types(minus),'entire ordinary cycle indices also coincide')
 need(all(permute(T,g)==T for g in plus+minus),'both are actual tournament subgroups')
 target=(4,7);source=(3,6)
 need(sum(tuple(g[x] for x in source)==target for g in plus)==1,'parallel transport exists')
 need(sum(tuple(g[x] for x in source)==target for g in minus)==0,'reversal transport forbids parallel matching')
 need(all(tuple(g[t[v]] for v in range(9))==tuple(t[g[v]] for v in range(9)) for g in minus),'centralizer realization guard')
 need(any(tuple(g[t[v]] for v in range(9))!=tuple(t[g[v]] for v in range(9)) for g in plus),'fixed-set law does not supply reversal guard')
 # Actual native complement, not a re-marking loop.
 selected=[]
 for P in paths(T):
  selected.append(P)
  if len(selected)==12:break
 for P in selected:
  B=core.complement(T,P)
  need(core.complement(B,P)==T,'native complement involution')
  for g in aut(2):
   need(core.complement(permute(T,g),tuple(g[v] for v in P))==permute(B,g),'native complement common-source covariance')
   need((permute(B,g)==B)==(g==tuple(range(9))),'aligned full native pair is rigid')
 return {'weighted_controls':report,'fixed_set_hostile':{'source':source,'target':target,'transporter_counts':[1,0],
          'common_fixed_set_distribution':[[list(k),v] for k,v in masks(plus).items()],
          'common_cycle_inventory':[[list(k),v] for k,v in cycle_types(plus).items()]},'native_paths_tested':len(selected)}

def partner(A,runs):
 B=list(core.converse(tuple(A)))
 for run in runs:
  for a,b in zip(run,run[1:]):B[a]^=1<<b;B[b]^=1<<a
 return tuple(B)
def endpoint_iso(B,U,C,V):
 return [g for g in oracle.isomorphisms(B,C,True) if all(g[a[0]]==b[0] and g[a[-1]]==b[-1] and len(a)==len(b) for a,b in zip(U,V))]
def cut(runs,i,j):return runs[:i]+[runs[i][:j],runs[i][j:]]+runs[i+1:]
def endpoint_tournament(n):
 T=[sum(1<<j for j in range(i+1,n)) for i in range(n)]
 T[0]^=1<<(n-1);T[n-1]^=1
 return tuple(T)
def assembly(A,runs,close=True):
 k=len(runs);q=k if close else k-1;Q3=(2,4,1);E=endpoint_tournament(q)
 T,_=core.substitute(Q3,(tuple(A),E,E));n=len(A)
 P=tuple(x for i,run in enumerate(runs) for x in list(run)+([n+i,n+q+i] if i<q else []))
 return T,P,core.complement(T,P)
def cuts_control():
 data=json.loads((HERE.parent/'cursor_cm_phase_20261028/endpoint_kernel.json').read_text())
 row=next(x for x in data['coupled_family'] if x['kind']=='length_aware' and len(x['A'])==8)
 A=tuple(row['A']);U=row['F'];V=row['second_order'];B=partner(A,U)
 need(B==tuple(row['B']) and partner(A,V)==B,'actual inherited aligned pair')
 before=endpoint_iso(B,U,B,V);need(bool(before),'one-stage length-aware responses agree')
 found=[]
 for i,run in enumerate(U):
  for j in range(1,len(run)):
   UU=cut(U,i,j);VV=cut(V,i,j);BU=partner(A,UU);BV=partner(A,VV)
   maps=endpoint_iso(BU,UU,BV,VV)
   if not maps:found.append({'run':i,'cut_after':j,'left_runs':UU,'right_runs':VV,'partners':[BU,BV]})
 need(bool(found),'legal run split defeats endpoint quotient')
 first=found[0];SU,PU,CU=assembly(A,first['left_runs']);SV,PV,CV=assembly(A,first['right_runs'])
 need(SU==SV,'same actual source for separating assemblies')
 need(all(core.edge(SU,a,b) for a,b in zip(PU,PU[1:])) and all(core.edge(SV,a,b) for a,b in zip(PV,PV[1:])),'actual Hamilton paths')
 H=[oracle.hamilton(CU),oracle.hamilton(CV)]
 small=[]
 for R in (first['left_runs'],first['right_runs']):
  S,P,C=assembly(A,R,False)
  need(core.is_strong(A) and core.is_strong(endpoint_tournament(len(R)-1)) and core.is_strong(S),'actual all-strong context')
  small.append({'source':S,'path':P,'H':oracle.hamilton(C)})
 need(small[0]['source']==small[1]['source'] and small[0]['H']!=small[1]['H'],'sixteen-vertex actual complement separation')
 G=oracle.isomorphisms(B,B,True);UU=cut(U,0,1);BB=partner(A,UU);GG=oracle.isomorphisms(BB,BB,True)
 def signature(R,group):return min(tuple((len(run),g[run[0]],g[run[-1]]) for run in R) for g in group)
 before_cells=set();after_cells=set()
 for order in permutations(U):
  R=list(order);before_cells.add(signature(R,G));i=next(i for i,r in enumerate(R) if len(r)>1)
  after_cells.add(signature(cut(R,i,1),GG))
 need((len(before_cells),len(after_cells),len(GG))==(8,24,1),'one intrinsic cut restores the full trit')
 # If equal, E2's intrinsic large decoder still proves class separation; report
 # the cheap context honestly instead of treating a failed scalar as a proof.
 return {'order':8,'source':A,'original_runs':[U,V],'initial_endpoint_transporters':before,
         'separating_cuts':found,'cheap_actual_context':{'order':len(SU),'source':SU,'paths':[PU,PV],
             'partner_H':H,'H_separates':H[0]!=H[1]},'smaller_context':small,
         'all_orderings_before_after':[len(before_cells),len(after_cells)],
         'scope':'one legal run split followed by an actual C3 native-complement context'}
def main(write=True):
 global CHECKS
 CHECKS=0
 inherited=frames.controls()
 out={'status':'FINITE-EXACT with declared menus','weights':controls(),'cuts':cuts_control(),'hostile':inherited['hostile'],
      'normal_basis_orientation_guards':{'F64_shifts':inherited['F64_shifts'],'F256_shifts':inherited['F256_shifts']},
      'nonflat_even_phase':inherited['holonomy']}
 out['checks']=CHECKS
 if write:(HERE/'ports.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
 print(json.dumps({'checks':CHECKS,'cut_positions':[r['cut_after'] for r in out['cuts']['separating_cuts']],
                   'native_partner_H':out['cuts']['cheap_actual_context']['partner_H'],'hostile':out['hostile']}))
 return out
if __name__=='__main__':main()
