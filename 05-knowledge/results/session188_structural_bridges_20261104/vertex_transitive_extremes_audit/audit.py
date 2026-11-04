from pathlib import Path
from itertools import combinations,permutations
import hashlib,json,subprocess,tempfile
HERE=Path(__file__).resolve().parent;SRC=HERE.parent/'vertex_transitive_extremes'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def op(T):return tuple(((1<<len(T))-1)^r^(1<<u) for u,r in enumerate(T))
def module(T,M):return all(all(bool(T[u]>>v&1)==bool(T[u]>>w&1) for v in range(len(T)) if M>>v&1 for w in range(len(T)) if M>>w&1) for u in range(len(T)) if not M>>u&1)
def mods(T):return [M for M in range(1,1<<len(T)) if M.bit_count()>1 and module(T,M)]
def induce(T,V):return tuple(sum(1<<j for j,v in enumerate(V) if T[u]>>v&1) for u in V)
def flip(T,M):return tuple(row^(M^(1<<u) if M>>u&1 else 0) for u,row in enumerate(T))
def isaut(T,p):return sorted(p)==list(range(len(T))) and all(bool(T[u]>>v&1)==bool(T[p[u]]>>p[v]&1) for u in range(len(T)) for v in range(len(T)))
def main():
 names=('PROOF.md','verify.py','count.cpp','normal.json','optimized.json');pins={n:digest(SRC/n) for n in names};D=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer bytes');need(D['status']=='PASS','producer status');need(all(digest(HERE.parent/p)==h for p,h in D['source_sha256'].items()),'producer source pins')
 regular=0
 for n in range(1,6):
  pairs=list(combinations(range(n),2))
  for bits in range(1<<len(pairs)):
   T=[0]*n
   for k,(u,v) in enumerate(pairs):
    if bits>>k&1:T[u]|=1<<v
    else:T[v]|=1<<u
   if len({x.bit_count() for x in T})!=1:continue
   regular+=1;ms=mods(T)
   for A in ms:
    U=induce(T,[u for u in range(n) if A>>u&1]);need(len({r.bit_count() for r in U})==1,'induced regularity')
   need(all(not(A&B) or A&B in(A,B) for A,B in combinations(ms,2)),'regular module laminarity')
 for rec in D['vertex_transitive_controls']+[D['nested_child']]:
  T=tuple(rec['rows']);need(len(rec['automorphisms_from_zero'])==len(T),'one automorphism per target');need(all(isaut(T,p) and p[0]==j for j,p in enumerate(rec['automorphisms_from_zero'])),'literal automorphism witness')
 profiles={};graphs={}
 def prof(T,expected=None):
  T=tuple(T)
  if T not in profiles:profiles[T]=expected
  elif expected is not None:need(profiles[T] in(None,expected),'profile consistency');profiles[T]=expected
 def query(T,a,c,z):
  k=(tuple(T),a,c);need(k not in graphs or graphs[k]==z,'literal query consistency');graphs[k]=z
 baseline=None;child=D['nested_child'];T=tuple(child['rows']);ms=mods(T);need(ms==child['modules'],'complete nested modules');seen=set()
 for f in child['all_family_marked_profiles']:
  U=T
  for M in f['selected_original_modules']:need(M in ms,'selected original nested module');U=flip(U,M)
  need(U==tuple(f['target_rows']),'nested XOR physical rows');expected=[f['q'],f['E'],f['S']];prof(U,expected);seen.add(tuple(f['selected_original_modules']));baseline=expected if baseline is None else baseline;need(expected==baseline,'retained complete marked invariance')
 need(len(seen)==1<<len(ms),'all nested families')
 for f in D['balanced_module_ambient_controls']:prof(f['T']);prof(f['U']);need(flip(tuple(f['T']),7)==tuple(f['U']),'ambient exact module flip')
 totals={'families':0,'scalar_strict':0,'equal_hub_selected':0,'equal_hub_unselected':0};restorations=0
 for rec in D['root_families']:
  T=tuple(rec['T']);n=len(T);m=len(rec['quotient']);children=rec['children'];off=[sum(map(len,children[:i])) for i in range(m)];blocks=[((1<<len(A))-1)<<off[i] for i,A in enumerate(children)];ms=[M for M in mods(T) if M!=(1<<n)-1];need(ms==rec['proper_modules'],'complete original root module list');need(all(any(M&B==M for B in blocks) for M in ms),'proper-module containment');at={tuple(q):k for k,q in enumerate(rec['queries'])};base=rec['baseline_values'];need(len(at)==len(base)<=m+20,'actual battery cardinality');need((0,0) in at,'scalar in quartic battery')
  for (a,c),z in zip(rec['queries'],base):query(T,a,c,z)
  familysets=set()
  for f in rec['families']:
   U=op(T)
   for M in f['selected_proper_modules']:need(M in ms,'original module selection');U=flip(U,M)
   need(U==tuple(f['target_rows']),'literal root XOR');familysets.add(tuple(f['selected_proper_modules']));vals=f['values'];need(vals!=base and vals[at[(0,0)]]>=base[at[(0,0)]],'detect and nondecrease');totals['families']+=1;totals[f['branch']]+=1
   for (a,c),z in zip(rec['queries'],vals):query(U,a,c,z)
   for j in(1,m-1):
    V=induce(U,range(off[j],off[j]+len(children[j])));prof(V);prof(children[j]);restorations+=1
   if f['branch']=='equal_hub_selected':
    need(blocks[0] in f['selected_proper_modules'],'selected branch');delta=sum(s*(base[at[(a,c)]]-vals[at[(a,c)]]) for a,c,s in rec['quartic']);need(delta==f['quartic_gap']>0,'quartic strictly positive')
   elif f['branch']=='equal_hub_unselected':
    need(blocks[0] not in f['selected_proper_modules'],'unselected branch');need([(base[at[tuple(q)]]-vals[at[tuple(q)]])%4 for q in rec['anchor']]==f['anchor_residues'] and 2 in f['anchor_residues'],'literal anchor residues')
  need(len(familysets)==1<<len(ms),'complete root family bank')
 with tempfile.TemporaryDirectory(prefix='vertex-transitive-audit-') as td:
  td=Path(td)
  for name in('profiles','literal'):subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/(name+'.cpp')),'-o',str(td/name)],check=True,capture_output=True)
  wire=str(len(profiles))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T)) for T in profiles)+'\n';lines=subprocess.check_output([str(td/'profiles')],input=wire.encode()).decode().splitlines();need(len(lines)==len(profiles),'profile outputs');actual={}
  for T,line in zip(profiles,lines):
   n=len(T);z=list(map(int,line.split()));need(len(z)==(1+2*n)*(n+1),'complete q/E/S shape');p=[z[:n+1],[z[(u+1)*(n+1):(u+2)*(n+1)] for u in range(n)],[z[(n+1+u)*(n+1):(n+2+u)*(n+1)] for u in range(n)]];actual[T]=p
   if profiles[T] is not None:need(p==profiles[T],'independent all-word/binomial-cut marked profiles')
  for f in D['balanced_module_ambient_controls']:need(actual[tuple(f['T'])]==actual[tuple(f['U'])],'independent ambient complete profile invariance')
  for rec in D['root_families']:
   children=rec['children'];off=[sum(map(len,children[:i])) for i in range(len(children))]
   for f in rec['families']:
    for j in(1,len(children)-1):need(actual[induce(f['target_rows'],range(off[j],off[j]+len(children[j])))]==actual[tuple(children[j])],'independent marked restoration in every root family')
  wire=str(len(graphs))+'\n'+'\n'.join(str(len(T))+' '+' '.join(map(str,T))+f' {a} {c}' for T,a,c in graphs)+'\n';out=list(map(int,subprocess.check_output([str(td/'literal')],input=wire.encode()).split()));need(out==list(graphs.values()),'independent predecessor-DP all actual physical counts')
 need(totals==D['totals'],'branch totals');need(pins=={n:digest(SRC/n) for n in names},'source freeze')
 print(json.dumps({'status':'PASS','regular_controls':regular,'independent_complete_word_cut_profiles':len(profiles),'nested_families':len(seen),'marked_family_restorations':restorations,'root_family_totals':totals,'literal_counts':len(graphs),'source_sha256':pins,'audit_sha256':{p.name:digest(p) for p in (Path(__file__),HERE/'profiles.cpp',HERE/'literal.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
