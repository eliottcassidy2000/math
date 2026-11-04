from pathlib import Path
from itertools import permutations,combinations
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'aggregated_cyclic_module'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def profiles(T):
 n=len(T);q=[0]*(n+1);E=[[0]*(n+1) for _ in T];S=[[0]*(n+1) for _ in T]
 for w in permutations(range(n)):
  bad=sum(1<<j for j in range(n-1) if not T[w[j]]>>w[j+1]&1)
  for cut in range(1<<(n-1)):
   if cut&bad!=bad:continue
   k=cut.bit_count()+1;t=(cut&-cut).bit_length()-1 if cut else n-1;q[k]+=1;E[w[t]][k]+=1;S[w[0]][k]+=1
 need(all(sum(e[k] for e in E)==sum(s[k] for s in S)==q[k] for k in range(n+1)),'complete prescribed-slot sum identity');return q,E,S
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,r in enumerate(A):
   T[off[i]+u]=r<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def submasks(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M
def main():
 names=('PROOF.md','verify.py','run_algebra.py','kernel.py','count.cpp','normal.json','optimized.json');pins={n:sha(SRC/n) for n in names};D=json.loads((SRC/'normal.json').read_text());need(D['status']=='PASS' and (SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired receipts');need(all(sha(BASE/p)==h for p,h in D['source_sha256'].items()),'producer source pins');queries={};child_profiles={};unbalanced=0
 def add(T,a,c,v):
  key=(tuple(T),a,c);need(key not in queries or queries[key]==v,'query consistency');queries[key]=v
 def profile(T):
  T=tuple(T)
  if T not in child_profiles:child_profiles[T]=profiles(T)
  return child_profiles[T]
 for rec in D['physical_controls']:
  X,Y=tuple(rec['X']),tuple(rec['Y']);px,py=profile(X),profile(Y);need(px[0]==rec['qX'] and py[0]==rec['qY'],'independent complete child q profiles');unbalanced+=int(px[1]!=px[2] or py[1]!=py[2]);n=len(rec['E']);Q=[2,4,1]+[r<<3 for r in rec['E']]
  for j in range(n):
   if rec['Aout']>>j&1:Q[0]|=1<<(j+3)
   else:Q[j+3]|=1
   if rec['Mout']>>j&1:Q[1]|=1<<(j+3);Q[2]|=1<<(j+3)
   else:Q[j+3]|=6
  T,off=blow(Q,[(2,4,1),X,Y]+[(0,)]*n);UL=flip(T,7);UQ,_=blow(op(Q),[(2,4,1),X,Y]+[(0,)]*n);need(T==tuple(rec['T']) and UL==tuple(rec['local_target']) and UQ==tuple(rec['quotient_target']),'physical local/quotient target construction')
  for key,U in [('local_queries',UL),('quotient_queries',UQ)]:
   need(len(rec[key])==16,'complete actual finite difference');delta=0
   for r in rec[key]:add(T,r['a'],r['c'],r['source']);add(U,r['a'],r['c'],r['target']);delta+=r['sign']*(r['source']-r['target'])
   need(delta==rec['positive_formula']['coefficient']>=rec['positive_formula']['lower_bound']>0,'exact aggregate positivity')
 totals={'families':0,'strict':0,'equal':0};cube=D['whole_child_cube']
 for rec in D['root_families']:
  T,off=blow(rec['Q'],rec['blocks']);need(T==tuple(rec['T']),'root source rows');n=len(T);mods=[M for M in range(1,(1<<n)-1) if M.bit_count()>1 and module(T,M)];need(mods==rec['proper_original_modules'],'complete original modules');need(all(any(M&B==M for B in rec['block_masks']) for M in mods),'proper module containment');battery=list(map(tuple,rec['battery']));at={q:i for i,q in enumerate(battery)};base=rec['baseline_values'];need(len(battery)<=len(rec['Q'])+20,'fixed battery bound');need(len(rec['families'])==1<<len(mods),'complete original root families')
  for (a,c),v in zip(battery,base):add(T,a,c,v)
  seen=set();hub=rec['block_masks'][rec['witness'][0]];hi=mods.index(hub)
  for f in rec['families']:
   bits=f['selection'];seen.add(bits);U=op(T)
   for j,M in enumerate(mods):
    if bits>>j&1:U=flip(U,M)
   need(U==tuple(f['target_rows']),'exact physical original-module XOR');z=f['values'];need(z!=base,'root-family fixed-battery detection');h=z[at[(0,0)]]-base[at[(0,0)]];need(h==f['scalar_increase']>=0,'scalar branch');totals['families']+=1;totals['strict' if h else 'equal']+=1
   for (a,c),v in zip(battery,z):add(U,a,c,v)
   for B,o,qnew in zip(rec['blocks'],off,f['target_child_profiles']):
    C=tuple((U[o+i]>>o)&((1<<len(B))-1) for i in range(len(B)));need(profile(C)[0]==qnew,'independent target child profile')
    if not h:need(qnew==profile(B)[0],'scalar-equality restoration')
   if not h and bits>>hi&1:need(sum(s*(base[at[(a,c)]]-z[at[(a,c)]]) for a,c,s in rec['quartic'])==f['quartic_gap']>0,'selected hub grouped quartic')
   elif not h:need([(base[at[tuple(q)]]-z[at[tuple(q)]])%4 for q in rec['anchor']]==f['anchor_residues'] and 2 in f['anchor_residues'],'unselected hub anchor')
  need(len(seen)==1<<len(mods),'unique complete family bank')
 for rec in cube['states']:
  bits=rec['bits'];Q=op((2,4,1)) if bits&1 else(2,4,1);B=[op(A) if bits>>(j+1)&1 else tuple(A) for j,A in enumerate(cube['blocks'])];T,_=blow(Q,B);need(T==tuple(rec['T']),'cube physical state');need(len(rec['values'])==len(cube['queries']),'cube battery shape')
  for (a,c),v in zip(cube['queries'],rec['values']):add(T,a,c,v)
 need(len(cube['states'])==16 and len({tuple(r['values']) for r in cube['states']})==16,'all complete cube vectors distinct')
 bank=D['complete_named_bank'];states=[bank['source'],bank['baseline'],bank['strict_target']];table={};coeff={}
 for r in bank['coordinates']:
  a,c=r['a'],r['c'];table[a,c]=r['values']
  for T,v in zip(states,r['values']):add(T,a,c,v)
 for r in bank['coordinates']:
  a,c=r['a'],r['c'];want=[sum((-1)**(a.bit_count()+c.bit_count()-s.bit_count()-t.bit_count())*table[s,t][j] for s in submasks(a) for t in submasks(c)) for j in range(3)];need(want==r['coefficients'],'every exact named cubic/quartic Mobius coordinate');coeff[a,c]=want
 for r in bank['full_mask_controls']:
  for j,T in enumerate(states):
   v=sum(z[j] for (a,c),z in coeff.items() if a&~r['a']==0 and c&~r['c']==0);need(v==r['values'][j],'full polynomial arbitrary-mask reconstruction');add(T,r['a'],r['c'],v)
 # Explicit independent alias on the complete unbalanced TT3 physical-name bank.
 alias=[];X=56;Y=64
 for j,T in enumerate(states):
  physical=sum(coeff[1|(1<<i),2|Y][j] for i in(3,4,5));grouped=0
  for bits in range(16):
   a=(1 if bits&1 else 0)|(X if bits&2 else 0);c=(2 if bits&4 else 0)|(Y if bits&8 else 0);v=sum(z[j] for (s,t),z in coeff.items() if s&~a==0 and t&~c==0);grouped+=(-1)**bits.bit_count()*v;add(T,a,c,v)
  need(physical==grouped,'full-degree exact alias under Boolean grouping');alias.append(grouped)
 with tempfile.TemporaryDirectory(prefix='aggregate-cyclic-audit-') as td:
  exe=Path(td)/'literal';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'literal.cpp'),'-o',str(exe)],check=True,capture_output=True);wire=str(len(queries))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in queries)+'\n';out=list(map(int,subprocess.check_output([str(exe)],input=wire.encode()).split()));need(out==list(queries.values()),'independent predecessor DP every physical query')
 need(pins=={n:sha(SRC/n) for n in names},'source freeze');need(totals=={'families':D['root_family_count'],'strict':D['scalar_strict_families'],'equal':D['scalar_equal_families']},'branch totals')
 print(json.dumps({'status':'PASS','physical_positive_controls':len(D['physical_controls']),'controls_with_unbalanced_marks':unbalanced,'independent_complete_q_E_S_profiles':len(child_profiles),'root_families':totals,'cube_states':len(cube['states']),'cube_actual_battery':len(cube['queries']),'complete_named_coordinates':len(coeff),'grouped_coefficients_from_all_named_terms':alias,'independent_physical_queries':len(queries),'source_sha256':pins,'audit_sha256':{p.name:sha(p) for p in(Path(__file__),HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
