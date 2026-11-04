from pathlib import Path
from itertools import permutations,combinations
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'quartic_cover_transfer'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def flip(T,M):return tuple(r^(M^(1<<u) if M>>u&1 else 0) for u,r in enumerate(T))
def module(T,M):return all(T[u]&M in(0,M) for u in range(len(T)) if not M>>u&1)
def induce(T,V):return tuple(sum(1<<j for j,v in enumerate(V) if T[u]>>v&1) for u in V)
def blow(Q,blocks):
 off=[sum(map(len,blocks[:i])) for i in range(len(blocks))];T=[0]*sum(map(len,blocks))
 for i,A in enumerate(blocks):
  for u,r in enumerate(A):
   T[off[i]+u]=r<<off[i]
   for j,B in enumerate(blocks):
    if Q[i]>>j&1:T[off[i]+u]|=((1<<len(B))-1)<<off[j]
 return tuple(T),off
def attach(T,a,c):
 n=len(T);U=list(T)+[1<<(n+1),0]
 for u in range(n):
  if a>>u&1:U[u]|=1<<n
  else:U[n]|=1<<u
  if c>>u&1:U[u]|=1<<(n+1)
  else:U[n+1]|=1<<u
 return tuple(U)
def ham(T):
 n=len(T);d=[[0]*n for _ in range(1<<n)]
 for u in range(n):d[1<<u][u]=1
 for M in range(1,1<<n):
  for u in range(n):
   if not M>>u&1:continue
   for v in range(n):
    if T[u]>>v&1 and not M>>v&1:d[M|1<<v][v]+=d[M][u]
 return sum(d[-1])
def marked(T,u):
 n=len(T);q=[0]*(n+1);E=[0]*(n+1);S=[0]*(n+1)
 for w in permutations(range(n)):
  bad=sum(1<<j for j in range(n-1) if not T[w[j]]>>w[j+1]&1)
  for cuts in range(1<<(n-1)):
   if cuts&bad!=bad:continue
   h=cuts.bit_count()+1;t=(cuts&-cuts).bit_length()-1 if cuts else n-1;q[h]+=1;E[h]+=w[t]==u;S[h]+=w[0]==u
 return q,E,S
def main():
 names=('PROOF.md','verify.py','count.cpp','normal.json','optimized.json');pins={n:sha(SRC/n) for n in names};D=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'producer paired receipts');need(all(sha(BASE/p)==h for p,h in D['source_sha256'].items()),'frozen source dependencies');P={};H={}
 def profile(T,p):
  T=tuple(T);need(T not in P or P[T]==p,'profile consistency');P[T]=p
 def query(T,a,c,v):
  key=(tuple(T),a,c);need(key not in H or H[key]==v,'Hamilton consistency');H[key]=v
 for r in D['cover_quartic_controls']:
  T=tuple(r['T']);U=op(T);A=[i for i in range(len(T)) if r['a_support']>>i&1];C=[i for i in range(len(T)) if r['c_support']>>i&1];co=[[0]*(len(T)+3) for _ in range(2)]
  need(len(A)==len(C)==2 and len(r['all16_profiles'])==16,'quartic support/battery')
  for bits,z in enumerate(r['all16_profiles']):
   a=sum(1<<v for j,v in enumerate(A) if bits>>j&1);c=sum(1<<v for j,v in enumerate(C) if bits>>(j+2)&1);need((a,c)==(z['a'],z['c']),'literal q profile masks');profile(attach(T,a,c),z['source']);profile(attach(U,a,c),z['converse'])
   for k,p in enumerate((z['source'],z['converse'])):
    for h,v in enumerate(p):co[k][h]+=(-1)**bits.bit_count()*v
  need(co[0]==co[1]==r['quartic_profile'] and all(v>=0 for v in co[0]),'all-h whole-converse quartic identity')
 for r in D['ambient_local_quotient_controls']:
  Q=[2,4,1]+[v<<3 for v in r['E']]
  for j in range(len(r['E'])):
   if r['Aout']>>j&1:Q[0]|=1<<(j+3)
   else:Q[j+3]|=1
   if r['Mout']>>j&1:Q[1]|=1<<(j+3);Q[2]|=1<<(j+3)
   else:Q[j+3]|=6
  B=[(2,4,1),r['X'],r['Y']]+[(0,)]*len(r['E']);T,off=blow(Q,B);U,_=blow(op(Q),B);W=flip(T,7);need([T,U,W]==list(map(tuple,[r['core_source'],r['core_quotient'],r['core_local']])),'literal core triple');outer=r['outer'];old=[blow(outer,[V]+[(0,)]*(len(outer)-1))[0] for V in(T,U,W)];need(old==list(map(tuple,r['ambient_rows'])),'exact fixed ambient rows');ext=sum(1<<(len(T)+j-1) for j in range(1,len(outer)) if outer[j]&1)
  px=marked(r['X'],r['marked_X']);py=marked(r['Y'],r['marked_Y']);bx,by=(px[0],py[0]) if r['grouped'] else(px[1],py[1]);need([bx,by]==r['marked_profiles'],'independent complete group/physical profiles')
  if not r['grouped']:need(px[1]==px[2] and py[1]==py[2],'physical balance')
  X=((1<<len(r['X']))-1)<<off[1] if r['grouped'] else 1<<(off[1]+r['marked_X']);Y=((1<<len(r['Y']))-1)<<off[2] if r['grouped'] else 1<<(off[2]+r['marked_Y']);Etotal=[0]+list(range(3+len(r['X'])+len(r['Y']),len(old[0])));lower=2*bx[1]*by[1]*ham(induce(old[0],Etotal));need(lower==r['strict_lower_bound']>0,'strict full ambient Hamilton bound');co=[0]*3
  for bits,z in enumerate(r['queries']):
   a=ext|(1 if bits&1 else 0)|(X if bits&2 else 0);c=ext|(2 if bits&4 else 0)|(Y if bits&8 else 0);need((a,c)==(z['a'],z['c']),'actual lifted physical masks')
   for k,V in enumerate(old):query(V,a,c,z['values'][k]);co[k]+=(-1)**bits.bit_count()*z['values'][k]
  need(co==r['quartics'] and co[0]-co[1]==co[0]-co[2]==r['positive_gap']>=lower,'exact ambient quotient/local coefficient transfer')
 for r in D['general_master_controls']:
  Q=r['Q'];A=r['distinguished_child'];T,off=blow(Q,[A,(2,0),(0,)]);U,_=blow(op(Q),[A,(2,0),(0,)]);W=flip(T,(1<<len(A))-1);outer=r['outer'];old=[blow(outer,[V]+[(0,)]*(len(outer)-1))[0] for V in(T,U,W)];need(old==list(map(tuple,r['ambient_rows'])),'general master physical triples');ext=sum(1<<(len(T)+j-1) for j in range(1,len(outer)) if outer[j]&1);q=0 if len(A)==1 else 1;X=3<<off[1];Y=1<<off[2];co=[0]*3
  for bits,z in enumerate(r['queries']):
   a=ext|(1 if bits&1 else 0)|(X if bits&2 else 0);c=ext|((1<<q) if bits&4 else 0)|(Y if bits&8 else 0);need((a,c)==(z['a'],z['c']),'generic saturated masks')
   for j,V in enumerate(old):query(V,a,c,z['values'][j]);co[j]+=(-1)**bits.bit_count()*z['values'][j]
  need(co==r['quartics'] and co[1]==co[2] and co[0]-co[1]==r['contrast'],'generic exact quotient/local transfer')
  if len(A)==1:need(co[0]==co[1],'singleton zero boundary')
 totals={'families':0,'strict':0,'equal':0}
 for r in D['ambient_root_families']:
  T,off=blow(r['Q'],r['blocks']);need(T==tuple(r['core_source']),'root core');mods=[M for M in range(1,(1<<len(T))-1) if M.bit_count()>1 and module(T,M)];need(mods==r['proper_modules'] and len(r['families'])==1<<len(mods),'complete original root families');outer=r['outer'];S,_=blow(outer,[T]+[(0,)]*(len(outer)-1));need(S==tuple(r['ambient_source']),'source strong ambient');ext=sum(1<<(len(T)+j-1) for j in range(1,len(outer)) if outer[j]&1);internal=list(map(tuple,r['internal_battery']));lifted=[(a|ext,c|ext) for a,c in internal];need(lifted==list(map(tuple,r['lifted_battery'])) and len(lifted)<=len(r['Q'])+20,'lifted battery and fixed bound');at={q:i for i,q in enumerate(internal)};base=r['baseline_values'];Z=((1<<len(r['blocks'][r['witness'][0]]))-1)<<off[r['witness'][0]];zi=mods.index(Z);hs=ham(T)
  for(a,c),v in zip(lifted,base):query(S,a,c,v)
  for f in r['families']:
   U=op(T)
   for j,M in enumerate(mods):
    if f['selection']>>j&1:U=flip(U,M)
   V,_=blow(outer,[U]+[(0,)]*(len(outer)-1));need(U==tuple(f['core_target']) and V==tuple(f['ambient_target']),'literal root XOR and fixed context');z=f['values'];need(z!=base,'ambient root detection');gap=ham(U)-hs;need(gap==f['core_scalar_increase']>=0,'independent core scalar branch');totals['families']+=1;totals['strict' if gap else 'equal']+=1
   for(a,c),v in zip(lifted,z):query(V,a,c,v)
   if gap:need(z[at[(0,0)]]-base[at[(0,0)]]==f['lifted_zero_gap']>0,'strict branch lifts at actual nonzero exterior mask')
   elif f['selection']>>zi&1:need(sum(s*(base[at[(a,c)]]-z[at[(a,c)]]) for a,c,s in r['quartic'])==f['quartic_gap']>0,'restored triangle grouped quartic')
   else:need([(base[at[tuple(q)]]-z[at[tuple(q)]])%4 for q in r['anchor']]==f['anchor_residues'] and 2 in f['anchor_residues'],'unrestored ambient anchor')
 cube=D['strong_ambient_cube'];need(len(cube['states'])==16 and len({tuple(s['values']) for s in cube['states']})==16,'complete distinct strong ambient cube')
 for r in cube['states']:
  S,_=blow((2,4,1),[r['core_rows'],(0,),(0,)]);need(S==tuple(r['ambient_rows']),'cube ambient physical rows')
  for(a,c),v in zip(cube['queries'],r['values']):query(S,a,c,v)
 with tempfile.TemporaryDirectory(prefix='quartic-cover-audit-') as td:
  td=Path(td)
  for name in('literal','word_profiles','profile_parity'):subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/(name+'.cpp')),'-o',str(td/name)],check=True,capture_output=True)
  wordwire=str(len(P))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T)) for T in P)+'\n';got=[list(map(int,s.split())) for s in subprocess.check_output([str(td/'word_profiles')],input=wordwire.encode()).decode().splitlines()];need(got==list(P.values()),'independent complete word/cut profile bank')
  wire=str(len(H))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in H)+'\n';got=list(map(int,subprocess.check_output([str(td/'literal')],input=wire.encode()).split()));need(got==list(H.values()),'independent predecessor DP all actual ambient counts');parity=json.loads(subprocess.check_output([str(td/'profile_parity')]));need(parity['status']=='PASS','exhaustive extra all-q check')
 need(totals=={'families':D['root_family_count'],'strict':D['strict_families'],'equal':D['equal_families']},'branch totals');need(pins=={n:sha(SRC/n) for n in names},'source freeze')
 print(json.dumps({'status':'PASS','all_q_reversal_exhaustive':parity,'independent_retained_full_profile_queries':len(P),'retained_profile_quartic_controls':len(D['cover_quartic_controls']),'ambient_local_quotient_controls':len(D['ambient_local_quotient_controls']),'generic_master_controls':len(D['general_master_controls']),'root_families':totals,'ambient_cube_states':len(cube['states']),'independent_actual_ambient_counts':len(H),'source_sha256':pins,'audit_sha256':{p.name:sha(p) for p in (Path(__file__),HERE/'literal.cpp',HERE/'word_profiles.cpp',HERE/'profile_parity.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
