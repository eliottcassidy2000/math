"""Exact disjoint-kernel, ordered-hub positivity and pair-free root-family controls."""
from pathlib import Path
from itertools import combinations,permutations,product
from math import comb
import hashlib,json,subprocess,tempfile,sys
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
from kernel import kernel,blow,flip,extend

def need(x,msg):
 if not x:raise RuntimeError(msg)
def rows(n,c):
 R=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if c>>k&1 else(j,i);R[u]|=1<<v
 return tuple(R)
def op(R):return tuple(((1<<len(R))-1)^a^(1<<i) for i,a in enumerate(R))
def make(A,R,out):
 a=len(A);n=len(R);T=list(A)+[v<<a for v in R]
 for i in range(a):
  for j in range(n):
   if out>>j&1:T[i]|=1<<(a+j)
   else:T[a+j]|=1<<i
 return tuple(T)
def profile(T):
 n=len(T);q=[0]*4
 for p in permutations(range(n)):
  for k in range(1,min(n,3)+1):
   for cuts in combinations(range(1,n),k-1):
    bounds=(0,)+cuts+(n,)
    if all(T[p[t]]>>p[t+1]&1 for u,v in zip(bounds,bounds[1:]) for t in range(u,v-1)):q[k]+=1
 return q
def factor(k,a,b):
 m=[0]*k;o=[0]*k
 if a is None:m[b]=1;o[:b]=[1]*b
 elif b is None:m[a]=1;o[a+1:]=[1]*(k-a-1)
 elif a==b:m[a]=1
 else:
  need(a<b,'ordered path');m[a]=m[b]=1;o[a+1:b]=[1]*(b-a-1)
 return m,o
def evaluate(profiles,*terms):
 k=len(profiles);out=1
 for v in range(k):
  m=sum(t[0][v] for t in terms);o=sum(t[1][v] for t in terms)
  out*=sum(comb(o,t)*(profiles[v][m+t] if m+t<4 else 0) for t in range(o+1))
 return out
def hub_formula(profiles):
 k=len(profiles);unit=([0]*k,[0]*k);Es=[factor(k,None,v) for v in range(k)];Ss=[factor(k,u,None) for u in range(k)];G=lambda u,v:factor(k,u,v);alpha=[unit]+[Es[v] for v in range(1,k,2)];beta=[unit]+[Ss[u] for u in range(0,k,2)];tau=[G(u,v) for u in range(0,k,2) for v in range(1,k,2) if u<v];rho0=G(0,0);rholast=[G(u,k-1) for u in range(0,k,2)];elllast=G(k-1,k-1);ell0=[G(0,v) for v in range(1,k,2)];g=G(0,k-1);ev=lambda *ts:evaluate(profiles,*ts)
 t1=2*sum(ev(a,b,g) for a in alpha for b in beta)
 t2=-ev(Es[0],Ss[-1])+ev(Es[-1],Ss[0])-2*sum(ev(t,Es[0],Ss[-1])-ev(t,Es[-1],Ss[0]) for t in tau)
 t3=sum(ev(a,rho0,Ss[-1])-sum(ev(a,r,Ss[0]) for r in rholast) for a in alpha)
 t4=sum(ev(b,Es[0],elllast)-sum(ev(b,Es[-1],e) for e in ell0) for b in beta)
 p1=sum(ev(a,g) for a in alpha)+sum(ev(b,g) for b in beta);p2=ev(Es[-1],Ss[0])+2*sum(ev(t,Es[-1],Ss[0]) for t in tau)
 middle=1
 for v in range(1,k-1):middle*=profiles[v][1]
 remainder=(profiles[0][1]*profiles[-1][2]+profiles[0][2]*profiles[-1][1])*middle
 if k==2:remainder+=profiles[0][1]*profiles[-1][1]
 need(sum((t1,t2,t3,t4))==p1+p2+remainder and min(p1,p2,remainder)>=0 and p1>0,'positive polynomial reduction')
 return {'kernel_terms':[t1,t2,t3,t4],'positive_terms':[p1,p2,remainder],'coefficient':p1+p2+remainder,'strict_lower_bound':2*profiles[0][1]*profiles[-1][1]*middle}
def W(k):
 Q=[0]*(k+1)
 for u in range(k):
  for v in range(u+1,k):Q[u+1]|=1<<(v+1)
  if u%2:Q[u+1]|=1
  else:Q[0]|=1<<(u+1)
 return tuple(Q)
def main():
 C3=(2,4,1);C5=tuple(sum(1<<((i+d)%5) for d in(1,2)) for i in range(5));p5=profile(C5);need(all(v%5==0 for v in p5),'vertex transitive normalization');b5=[v//5 for v in p5]
 families=[];graphs=[];ids={}
 def add(G):
  G=tuple(G)
  if G not in ids:ids[G]=len(graphs);graphs.append(G)
  return ids[G]
 def quartic(T,U,anames,cnames):
  data=[]
  for bits in range(16):
   a=sum(1<<anames[t] for t in range(2) if bits>>t&1);c=sum(1<<cnames[t] for t in range(2) if bits>>(2+t)&1)
   data.append({'a':a,'c':c,'sign':(-1)**(4-bits.bit_count()),'source_id':add(extend(T,a,c)),'target_id':add(extend(U,a,c))})
  return data
 local=[]
 for n in range(2,6):
  codes=range(2) if n==2 else(range(8) if n<=4 else(0,1,81,233,601,1023))
  for code in codes:
   R=rows(n,code)
   for out in(range(4) if n==2 else((code*5+3)&((1<<n)-1),)):
    T=make(C3,R,out);U=flip(T,7);terms=kernel(R,[v for v in range(n) if out>>v&1],[v for v in range(n) if not out>>v&1],0,1)
    local.append({'R':R,'out_mask':out,'T':T,'U':U,'kernel_terms':terms,'coefficient':sum(terms),'queries':quartic(T,U,(0,3),(1,4))})
 zero=[]
 for A in (C3,(6,4,0),(14,12,8,0),C5):
  for code in(0,3,7):
   R=rows(3,code);T=make(A,R,(code+1)%8);U=flip(T,(1<<len(A))-1);zero.append({'A':A,'R':R,'T':T,'U':U,'queries':quartic(T,U,(0,1),(len(A),len(A)+1))})
 root=[]
 for Q in((2,4,1),(14,12,8,0),(6,12,24,17,3)):
  T=blow(Q,[C3]*len(Q));U=blow(op(Q),[C3]*len(Q));R=tuple(row>>3 for row in T[3:]);out=[v for v in range(len(R)) if T[0]>>(v+3)&1];inn=[v for v in range(len(R)) if T[v+3]&1];terms=kernel(R,out,inn,0,3);root.append({'Q':Q,'T':T,'U':U,'kernel_terms':terms,'coefficient':sum(terms),'queries':quartic(T,U,(0,3),(1,6))})
 hub=[]
 specs=[(2,[C3]*2),(4,[C3]*4),(6,[C3]+[(0,)]*4+[C3]),(8,[C3]+[(0,)]*6+[C3]),(2,[C5,C5]),(4,[C5,(0,),(0,),C5])]
 for k,children in specs:
  Q=W(k);T=blow(Q,[C3]+children);U=blow(op(Q),[C3]+children);i=3;j=3+sum(map(len,children[:-1]));profiles=[profile(D) for D in children];profiles[0]=[v//len(children[0]) for v in profiles[0]];profiles[-1]=[v//len(children[-1]) for v in profiles[-1]];f=hub_formula(profiles);hub.append({'chain_children':children,'Q':Q,'T':T,'U':U,'marked_names':[0,i,1,j],'profile_formula':f,'queries':quartic(T,U,(0,i),(1,j))})
 positive=[]
 for r in range(1,21):
  k=2*r;p=[[0,3,6,6] for _ in range(k)];p[0]=p[-1]=[0,1,2,2];positive.append({'r':r,'all_C3':hub_formula(p)})
 for r in range(1,11):
  k=2*r;p=[[0,1,2**(1+t%7)-2,3**(1+t%7)-3*2**(1+t%7)+3] for t in range(k)];p[0]=[0,1,2,2];p[-1]=b5;positive.append({'r':r,'mixed_profiles':p,'formula':hub_formula(p)})
 # Complete root-containing families for two smallest W quotients, all C3.
 for k in(2,4):
  Q=W(k);m=k+1;T=blow(Q,[C3]*m);blocks=[7<<(3*i) for i in range(m)];p=0;compressed=[1<<1,1<<2]+blocks[1:];d=sum(bool(T[p]&B) for B in compressed);anchor=[(B|((1<<p) if T[p]&B else 0),1<<p) for B in compressed]
  if d%2:anchor.append((1<<p,1<<p))
  if(len(compressed)+d)%2:anchor.append((0,1<<p))
  anchor.append((1<<p,0));quart=[(sum(1<<v for t,v in enumerate((0,3)) if bits>>t&1),sum(1<<v for t,v in enumerate((1,3*k)) if bits>>(2+t)&1),(-1)**(4-bits.bit_count())) for bits in range(16)];battery=list(dict.fromkeys([(a,c) for a,c,s in quart]+anchor));need(len(battery)<=m+20,'compressed fixed root-family battery')
  f={'Q':Q,'T':T,'battery':battery,'quartic':quart,'anchor':anchor,'baseline_ids':[add(extend(T,a,c)) for a,c in battery],'records':[]}
  for bits in range(1<<m):
   U=flip(T,(1<<len(T))-1)
   for z in range(m):
    if bits>>z&1:U=flip(U,blocks[z])
   f['records'].append({'selected_children_bits':bits,'target_rows':U,'value_ids':[add(extend(U,a,c)) for a,c in battery]})
  families.append(f)
 with tempfile.TemporaryDirectory() as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True)
  wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(len(values)==len(graphs),'all literal values')
 for category in(local,zero,root,hub):
  for rec in category:
   total=0
   for q in rec['queries']:q['source']=values[q.pop('source_id')];q['target']=values[q.pop('target_id')];total+=q['sign']*(q['source']-q['target'])
   want=0 if category is zero else(rec['profile_formula']['coefficient'] if category is hub else rec['coefficient']);need(total==want,'complete literal disjoint quartic');rec['literal_coefficient']=total
 for f in families:
  b=f['baseline_values']=[values[i] for i in f.pop('baseline_ids')];index={tuple(q):i for i,q in enumerate(f['battery'])}
  for rec in f['records']:
   v=rec['values']=[values[i] for i in rec.pop('value_ids')];need(b!=v,'pair-free root-family response collision')
   if rec['selected_children_bits']&1:
    gap=sum(s*(b[index[a,c]]-v[index[a,c]]) for a,c,s in f['quartic']);need(gap>0,'selected hub quartic');rec['quartic_gap']=gap
   else:
    residues=[(b[index[a,c]]-v[index[a,c]])%4 for a,c in f['anchor']];need(2 in residues,'unselected hub compressed anchor');rec['anchor_residues']=residues
 paths=[ROOT/n for n in('PROOF.md','verify.py','kernel.py','count.cpp')]+[BASE/'arbitrary_third_root_family'/'verify.py',BASE/'response'/'PROOF.md',BASE/'returns_next'/'UNIFORM.md',BASE/'alternating_prime_quotient'/'PROOF.md']
 print(json.dumps({'status':'PASS','literal_augmented_tournaments':len(graphs),'local_kernel_controls':local,'same_port_zero_controls':zero,'balanced_root_controls':root,'ordered_hub_literal_controls':hub,'factorized_positive_controls':positive,'pair_free_root_families':families,'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
