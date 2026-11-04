"""Independent W-chain positivity and all-C3 root-family audit, no producer imports."""
from pathlib import Path
from itertools import combinations,permutations,product
import hashlib,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'cyclic_child_disjoint'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def plus(*ps):
 d={}
 for p in ps:
  for e,z in p.items():d[e]=d.get(e,0)+z
 return {e:z for e,z in d.items() if z}
def scale(p,c):return {e:c*z for e,z in p.items() if c*z}
def mul(p,q):
 d={}
 for a,x in p.items():
  for b,y in q.items():
   e=tuple(u+v for u,v in zip(a,b));d[e]=d.get(e,0)+x*y
 return {e:z for e,z in d.items() if z}
def factor(k,mandatory,optional):
 out={}
 for bits in range(1<<len(optional)):
  e=[0]*k
  for u in mandatory:e[u]+=1
  for t,u in enumerate(optional):e[u]+=bool(bits>>t&1)
  out[tuple(e)]=1
 return out
def symbolic(k):
 one={(0,)*k:1};E=lambda v:factor(k,[v],list(range(v)));S=lambda u:factor(k,[u],list(range(u+1,k)));G=lambda u,v:factor(k,[u] if u==v else[u,v],list(range(u+1,v)))
 alpha=plus(one,*[E(v) for v in range(1,k,2)]);beta=plus(one,*[S(u) for u in range(0,k,2)]);tau=plus(*[G(u,v) for u in range(0,k,2) for v in range(1,k,2) if u<v]);g=G(0,k-1);e=E(k-1);s=S(0);ends=mul(E(0),S(k-1));rho=plus(*[G(u,k-1) for u in range(0,k,2)]);ell=plus(*[G(0,v) for v in range(1,k,2)])
 original=plus(scale(mul(mul(alpha,beta),g),2),scale(mul(plus(one,scale(tau,2)),plus(ends,scale(mul(e,s),-1))),-1),mul(alpha,plus(ends,scale(mul(rho,s),-1))),mul(beta,plus(ends,scale(mul(e,ell),-1))))
 first=mul(plus(alpha,beta),g);second=mul(plus(one,scale(tau,2)),mul(e,s));last=mul(plus(alpha,beta,scale(one,-1),scale(tau,-2)),ends);need(original==plus(first,second,last),'complete polynomial cancellation')
 full=lambda p:{a:c for a,c in p.items() if all(a)};want={}
 for l,r in ((1,2),(2,1)):
  a=[1]*k;a[0]=l;a[-1]=r;want[tuple(a)]=1
 if k==2:want[(1,1)]=1
 need(full(last)==want,'complete-support last-term polynomial, including k2');need(all(c>=0 for p in (first,second,full(last)) for c in p.values()),'coefficientwise positivity after full-support projection');need(full(first).get((1,)*k,0)>=2,'strict all-once lower bound')
 packed=json.dumps([[*a,c] for a,c in sorted(original.items())],separators=(',',':')).encode();return {'chain_order':k,'full_kernel_monomials':len(original),'complete_support_monomials':len(full(original)),'kernel_sha256':hashlib.sha256(packed).hexdigest(),'projected_last_terms':[[*a,c] for a,c in sorted(full(last).items())]}
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def W(k):
 return tuple(sum(1<<v for v in range(k+1) if (u and v and u<v) or(u==0 and v and (v-1)%2==0) or(v==0 and u and (u-1)%2==1)) for u in range(k+1))
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def marked(T):
 n=len(T);q=[0]*(n+1);E=[[0]*(n+1) for _ in T];S=[[0]*(n+1) for _ in T]
 for word in permutations(range(n)):
  bad=sum(1<<j for j,(u,v) in enumerate(zip(word,word[1:])) if not T[u]>>v&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   h=1+cuts.bit_count();last=(cuts&-cuts).bit_length()-1 if cuts else n-1;q[h]+=1;E[word[last]][h]+=1;S[word[0]][h]+=1
 return q,E,S
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def main():
 names=('PROOF.md','verify.py','kernel.py','count.cpp');paths=[SRC/n for n in names]+[BASE/'response'/'PROOF.md',BASE/'returns_next'/'UNIFORM.md',BASE/'alternating_prime_quotient'/'PROOF.md'];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());polys=[symbolic(k) for k in(2,4,6,8)];queries={};hub=[];profiles_checked=0
 def add(T,a,c,z):
  key=(tuple(T),a,c);need(key not in queries or queries[key]==z,'consistent retained query');queries[key]=z
 for rec in data['ordered_hub_literal_controls']:
  children=list(map(tuple,rec['chain_children']));k=len(children);Q=W(k);need(Q==tuple(rec['Q']),'independent W physical quotient');T,off=blow(Q,[(2,4,1)]+children);U,_=blow(op(Q),[(2,4,1)]+children);need(T==tuple(rec['T']) and U==tuple(rec['U']),'physical quotient-converse blowups');b=[];middle=1
  for t,D in enumerate(children):
   q,E,S=marked(D);profiles_checked+=1
   if t in(0,k-1):need(E[0]==S[0] and E[0][1]>0,'marked endpoint balance and positive first coordinate');b.append(E[0])
   else:middle*=q[1]
  last=(b[0][1]*(b[1][2] if len(b[1])>2 else 0)+(b[0][2] if len(b[0])>2 else 0)*b[1][1])*middle+(b[0][1]*b[1][1] if k==2 else 0);lower=2*b[0][1]*b[1][1]*middle;f=rec['profile_formula'];need(f['positive_terms'][-1]==last and f['strict_lower_bound']==lower and f['coefficient']>=lower>0,'exact last functional and strict bound');delta=0
  need(len(rec['queries'])==16,'complete disjoint quartic battery')
  for q in rec['queries']:
   add(T,q['a'],q['c'],q['source']);add(U,q['a'],q['c'],q['target']);delta+=q['sign']*(q['source']-q['target'])
  need(delta==rec['literal_coefficient']==f['coefficient'],'complete positive literal coefficient');hub.append({'chain_order':k,'physical_order':len(T),'coefficient':delta,'lower_bound':lower,'last_term':last})
 family_records=[]
 for rec in data['pair_free_root_families']:
  Q=tuple(rec['Q']);k=len(Q)-1;m=k+1;need(Q==W(k),'family W convention');T,_=blow(Q,[(2,4,1)]*m);need(T==tuple(rec['T']),'family physical source');n=len(T);full=(1<<n)-1;blocks=[7<<(3*i) for i in range(m)];mods=[M for M in range(1,full) if M.bit_count()>1 and module(T,M)];need(mods==sorted(blocks),'exact proper nontrivial original modules');p=0;compressed=[2,4]+blocks[1:];d=sum(bool(T[p]&B) for B in compressed);anchor=[(B|(1 if T[p]&B else 0),1) for B in compressed]
  if d%2:anchor.append((1,1))
  if(len(compressed)+d)%2:anchor.append((0,1))
  anchor.append((1,0));need(anchor==list(map(tuple,rec['anchor'])),'compressed anchor physical masks');battery=list(map(tuple,rec['battery']));need(len(battery)<=m+20 and set(battery)=={(a,c) for a,c,s in rec['quartic']}|set(anchor),'fixed battery completeness/bound');at={q:i for i,q in enumerate(battery)};base=rec['baseline_values']
  for(a,c),v in zip(battery,base):add(T,a,c,v)
  need(len(rec['records'])==1<<m,'all root-containing child subsets retained');seen=set();selected=unselected=0
  for f in rec['records']:
   bits=f['selected_children_bits'];need(bits not in seen,'unique selected family');seen.add(bits);U=list(op(T))
   for i,M in enumerate(blocks):
    if bits>>i&1:
     for u in range(n):
      if M>>u&1:U[u]^=M^(1<<u)
   U=tuple(U);need(U==tuple(f['target_rows']),'exact root/child XOR with physical labels');vals=f['values'];need(vals!=base,'fixed battery detection')
   for(a,c),v in zip(battery,vals):add(U,a,c,v)
   if bits&1:
    selected+=1;gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in rec['quartic']);need(gap==f['quartic_gap'] and gap>0,'selected hub positive disjoint quartic')
   else:
    unselected+=1;need([(base[at[q]]-vals[at[q]])%4 for q in anchor]==f['anchor_residues'] and 2 in f['anchor_residues'],'unselected hub whole-converse anchor')
  family_records.append({'quotient_order':m,'families':len(seen),'hub_selected':selected,'hub_unselected':unselected,'fixed_query_count':len(battery)})
 inputs=[str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in queries]
 with tempfile.TemporaryDirectory(prefix='cyclic-hub-audit-') as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True,capture_output=True);out=subprocess.check_output([str(exe)],input=(str(len(inputs))+'\n'+'\n'.join(inputs)+'\n').encode());need(list(map(int,out.split()))==list(queries.values()),'independent predecessor-DP actual physical counts')
 need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired receipts agree');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','scope':'Independent Sections4–6 positivity and all-C3 root-family/cube audit; Sections1–3 are reviewed separately.','complete_color_polynomial_checks':polys,'explicit_marked_child_profiles':profiles_checked,'ordered_hub_literal_cases':hub,'root_family_controls':family_records,'independent_physical_counts':len(queries),'source_sha256':pins,'audit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in(Path(__file__),HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
