"""Independent XOR/profile-restoration/fixed-battery audit with literal counts."""
from pathlib import Path
from itertools import combinations,permutations,product
from functools import lru_cache
import hashlib,json,subprocess,sys,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'prime_pair_root_family'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def induce(T,vs):return tuple(sum(1<<j for j,v in enumerate(vs) if T[u]>>v&1) for u in vs)
@lru_cache(None)
def profile(T):
 n=len(T);q=[0]*(n+1)
 for word in permutations(range(n)):
  bad=sum(1<<j for j,(u,v) in enumerate(zip(word,word[1:])) if not T[u]>>v&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad==bad:q[1+cuts.bit_count()]+=1
 return q
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,row in enumerate(A):
   T[off[i]+u]=row<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off

def main():
 paths=[SRC/name for name in ('PROOF.md','verify.py','count.cpp')]+[BASE/d/'PROOF.md' for d in ('module_decomposition','marked_pair_quotient','ordered_module_skew','critical_quotient_synthesis')]+[BASE/'returns_next'/'UNIFORM.md',BASE/'arbitrary_third_root_family'/'verify.py'];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths};data=json.loads((SRC/'normal.json').read_text());counts={'families':0,'scalar_strict':0,'equal_pair_selected':0,'equal_pair_not_selected':0};queries={};case_records=[];restorations=0
 def add(G,a,c,z):
  key=(tuple(G),a,c);need(key not in queries or queries[key]==z,'consistent physical count');queries[key]=z
 for rec in data['cases']:
  Q=tuple(rec['quotient']);children=list(map(tuple,rec['children']));T,off=blow(Q,children);need(T==tuple(rec['source_rows']),'physical source blow-up');fixed,_=blow(op(Q),children);n=len(T);full=(1<<n)-1;m=len(Q);z=rec['root_pair_child'];p,q=off[z],off[z]+1;blocks=[((1<<len(A))-1)<<off[i] for i,A in enumerate(children)];need(children[z]==(2,0),'named p to q pair');mods=[M for M in range(1,full) if M.bit_count()>1 and module(T,M)];need(mods==rec['proper_modules'],'complete independent original module list');need(all(any(M&B==M for B in blocks) for M in mods),'no proper cross-child module');need(len(rec['families'])==1<<len(mods),'all proper-module subsets retained');I,J=rec['marked_exterior_child_sets'];X=sum(blocks[i] for i in I);Y=sum(blocks[j] for j in J);need(X and Y and not X&Y and not(X|Y)&blocks[z] and (X|Y).bit_count()>=3,'exterior module size/scope');need(all(T[u]&Y==Y for u in range(n) if X>>u&1),'ordered X to Y');need(all(T[u]&(X|Y) in(0,X|Y) for u in range(n) if not (blocks[z]|X|Y)>>u&1),'exterior module homogeneity')
  quartet=[((bits[0]<<q)|(X if bits[1] else 0),(bits[2]<<p)|(Y if bits[3] else 0),(-1)**sum(bits)) for bits in product((0,1),repeat=4)];need(quartet==list(map(tuple,rec['quartic'])),'four-bit coefficient battery');other=[B for i,B in enumerate(blocks) if i!=z]+[1<<q];d=sum(bool(T[p]&B) for B in other);anchor=[(B|((1<<p) if T[p]&B else 0),1<<p) for B in other]
  if d%2:anchor.append((1<<p,1<<p))
  if (len(other)+d)%2:anchor.append((0,1<<p))
  anchor.append((1<<p,0));need(anchor==list(map(tuple,rec['anchor'])),'compressed block-count anchor formula');qs=list(map(tuple,rec['queries']));need(set(qs)=={(a,c) for a,c,s in quartet}|set(anchor),'complete fixed union battery');need(len(qs)<=m+18,'odd critical quotient sharpened query bound');base=rec['baseline_values'];at={a:i for i,a in enumerate(qs)};h=base[at[(0,0)]]
  for (a,c),v in zip(qs,base):add(T,a,c,v)
  seen=set();branches={k:0 for k in counts if k!='families'}
  for f in rec['families']:
   sel=tuple(f['selected_proper_modules']);need(sel not in seen and all(M in mods for M in sel),'unique original-module family');seen.add(sel);U=list(op(T))
   for M in sel:
    for u in range(n):
     if M>>u&1:U[u]^=M^(1<<u)
   U=tuple(U);need(U==tuple(f['target_rows']),'independent simultaneous XOR including root');need((blocks[z] in sel)==f['pair_selected'],'pair selection parity');need(induce(U,(p,q))==((2,0) if f['pair_selected'] else (0,1)),'pair final orientation');vals=f['values'];g=vals[at[(0,0)]];need(g>=h and vals!=base,'scalar monotonicity and fixed-battery detection');counts['families']+=1
   for (a,c),v in zip(qs,vals):add(U,a,c,v)
   if g>h:branch='scalar_strict'
   else:
    for i,A in enumerate(children):need(profile(A)==profile(induce(U,range(off[i],off[i]+len(A)))),'hereditary child profile equality')
    if f['pair_selected']:
     branch='equal_pair_selected';gap=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in quartet);need(gap==f['quartic_gap'] and gap>0,'positive marked quartic')
     for a,c,s in quartet:add(fixed,a,c,vals[at[(a,c)]]);restorations+=1
    else:
     branch='equal_pair_not_selected';res=[(base[at[(a,c)]]-vals[at[(a,c)]])%4 for a,c in anchor];need(res==f['anchor_residues_mod4'] and 2 in res,'whole-converse anchor residue')
     for a,c in anchor:add(op(T),a,c,vals[at[(a,c)]]);restorations+=1
   need(f['branch']==branch,'recorded branch');counts[branch]+=1;branches[branch]+=1
  case_records.append({'label':rec['label'],'old_order':n,'proper_module_count':len(mods),'fixed_query_count':len(qs),'family_count':len(seen),'branches':branches})
 need(counts==data['totals'],'independent branch totals');inputs=[str(len(G))+' '+' '.join(map(str,G))+f' {a} {c}' for G,a,c in queries]
 with tempfile.TemporaryDirectory(prefix='prime-pair-audit-') as tmp:
  exe=Path(tmp)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True,capture_output=True);out=subprocess.check_output([str(exe)],input=(str(len(inputs))+'\n'+'\n'.join(inputs)+'\n').encode());need(list(map(int,out.split()))==list(queries.values()),'independent predecessor DP all physical counts')
 normal=subprocess.check_output([sys.executable,str(SRC/'verify.py')]);optimized=subprocess.check_output([sys.executable,'-O',str(SRC/'verify.py')]);need(normal==optimized==(SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer modes replay');need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},'source freeze')
 print(json.dumps({'status':'PASS','totals':counts,'independent_cases':case_records,'independent_literal_queries':len(queries),'uniform_profile_restoration_controls':restorations,'producer_modes_equal':True,'source_sha256':pins,'audit_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
