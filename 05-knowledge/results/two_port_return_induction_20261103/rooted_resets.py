#!/usr/bin/env python3
"""Exact rooted reset interfaces, amplification and one bounded finite grammar.
All proof scopes and the virtual final birth are documented in ROOTED_RESETS.md.
"""
from pathlib import Path
import importlib.util,json,sympy as S,itertools
HERE=Path(__file__).resolve().parent
p=HERE.parent/'adaptive_state_induction_20261103';sp=importlib.util.spec_from_file_location('L',p/'live_port.py');L=importlib.util.module_from_spec(sp);sp.loader.exec_module(L)
states,Ms,seed,f=L.compile_carrier();ix={s:i for i,s in enumerate(states)};d=json.loads((p/'live_port.json').read_text());As=[S.Matrix(m) for m in d['rank']['query_matrices']]
def query(suffix):
 v=f
 for a in reversed(suffix):v=L.advance(v,tuple((j,i,c) for i,j,c in Ms[a]))
 return v
C=S.Matrix.hstack(*[S.Matrix(query(w)) for w in d['rank']['suffixes']]);old=[((-1,-1),),((-1,0),),((0,-1),),((-1,-1),(0,0)),((-1,-1),(-1,0))]
B=S.zeros(5,24)
for i,s in enumerate(old):
 for t,c in L.kernel(s,0,1):B[i,ix[t]]+=c
J=S.Matrix([[1,0,0,0,0],[-1,0,1,0,0],[0,0,0,1,-1],[0,0,0,0,1],[-1,1,0,0,-1]])
F=(J*B*C)[:4,:];P=S.Matrix([[1,0,0,1],[0,0,1,1],[1,1,1,1],[0,0,0,1]]);F=P*F
swap=S.zeros(24)
for i,state in enumerate(states):
 sig=tuple(sorted((1-a if a>=0 else a,1-b if b>=0 else b) for a,b in state));swap[i,ix[sig]]=1
rows=list(C.T.rref()[1]);AS=C[rows,:].inv()*(swap*C)[rows,:];L.need(C*AS==swap*C, 'slot permutation intertwining')
cols=list(F.rref()[1]);inv=F[:,cols].inv()
def macro(w):
 A=S.eye(17)
 for a in w:A=A*As[a]
 if w[-1]==3:A=A*AS
 FA=F*A;M=FA[:,cols]*inv;L.need(M*F==FA, 'macro factorization')
 return M

def ints(m):
 L.need(all(x.q == 1 for x in m), 'integer matrix')
 return tuple(tuple(int(x) for x in row) for row in m.tolist())

def step(v, m):
 return tuple(sum(v[i]*m[i][j] for i in range(len(v))) for j in range(len(m[0])))

def positive(z):
 h,s,q,f=z
 return (h-s,q-s,s,f-h-q+s)

def apply_word(T, ports, word):
 ports=list(ports)
 for a in word:
  r,drop=divmod(a,3);n=len(T);T=L.birth(T,ports[r])
  if drop<2:ports[drop]=n
 return T,tuple(ports)

def compiler(entry, words):
 # Controller slot renaming is compiled into letters, never a graph move.
 raw=[];flip=0
 for word in (entry,)+tuple(words):
  for a in word:
   r,drop=divmod(a,3)
   raw.append(3*(r^flip)+(drop^flip if drop<2 else 2))
  if word[-1]==3:flip^=1
 child,ports,roots=L.decode(raw)
 root=ports[flip]
 L.need(ports[1^flip]==len(child)-1, 'final virtual newborn')
 parent,_,parentroots=L.decode(raw[:-1])
 L.need(parent==tuple(r&((1<<len(parent))-1) for r in child[:-1]), 'omit final birth')
 L.need(root==roots[-1], 'root retained at boundary')
 return parent,root,tuple(raw),parentroots

def amplification():
 f=S.eye(17)[:,0];gs=[(As[a]**2-As[a])*f for a in (2,5)]
 cols=[];words=[]
 def add(v,word):
  if S.Matrix.hstack(*(cols+[v])).rank()>len(cols):cols.append(v);words.append(word)
 for r,v in enumerate(gs):add(v,(r,()))
 j=0
 while j<len(cols):
  for a,A in enumerate(As):add(A*cols[j],(words[j][0],(a,)+words[j][1]))
  j+=1
 G=S.Matrix.hstack(*cols)
 L.need(G.shape==(17,17) and G.det()==128, 'full amplification certificate')
 return {'tests':[{'root_slot':r,'prefix':w} for r,w in words],
         'matrix':[[str(x) for x in row] for row in G.tolist()],
         'determinant':int(G.det()),'maximum_prefix_length':max(len(w) for _,w in words)}

def controls():
 sp=importlib.util.spec_from_file_location('root_profile_source',p.parent/'adaptive_metagraph_reformulation_20261103/adaptive_probe.py')
 old=importlib.util.module_from_spec(sp);sp.loader.exec_module(old)
 sources=[]
 for n in range(1,4):
  for code in range(1<<(n*(n-1)//2)):sources.append(old.decode(n,code))
 sources += [L.birth((2,4,1),0), (6,12,24,17,3)]
 checked=0;actual_profiles=[]
 for T in sources:
  for r,z in enumerate(old.profiles(T)):
   v=positive(z);actual_profiles.append(v);pred=step(v,ints(F));child=L.birth(T,r);ports=(r,len(T))
   for i,w in enumerate(d['rank']['suffixes']):
    graph,_=apply_word(child,ports,w)
    L.need(L.count(graph)==pred[i], 'literal root-child continuation')
    checked+=1
  old.profiles.cache_clear();old.table.cache_clear()
 L.need(S.Matrix(actual_profiles).rank()==4,'actual rooted profiles span four dimensions')
 return {'parent_sources':len(sources),'actual_profile_rank':4,'literal_continuation_counts':checked,
         'includes_singleton_boundary':True,'maximum_parent_order':5}

def grammar():
 from collections import deque
 from hashlib import sha256
 models={}
 for k in range(3):
  for pre in itertools.product((0,2,4,5),repeat=k):
   for a in (1,3):
    w=pre+(a,);M=macro(w);models.setdefault(tuple(M),(w,ints(M)))
 models=list(models.values())
 L.need(len(models)==28 and all(x>=0 for _,M in models for row in M for x in row), '28 positive macro maps')
 seedX=S.Matrix([d['rank']['query_seed']]);results=[]
 sp=importlib.util.spec_from_file_location('compiler_profile_control',p.parent/'adaptive_metagraph_reformulation_20261103/adaptive_probe.py')
 checker=importlib.util.module_from_spec(sp);sp.loader.exec_module(checker)
 for radius in (1,2):
  seeds={}
  for k in range(radius+1):
   for pre in itertools.product((0,2,4,5),repeat=k):
    X=seedX
    for a in pre:X=X*As[a]
    for a in (1,3):
     Y=X*As[a]
     if a==3:Y=Y*AS
     v=Y[:,cols]*inv;L.need(v*F==Y,'entry factor')
     seeds.setdefault(tuple(int(x) for x in v),pre+(a,))
  compiled_profiles=0
  if radius==2:
   for v,entry in seeds.items():
    for words,target in [((),v)]+[((w,),step(v,M)) for w,M in models]:
     graph,root,_,_=compiler(entry,words)
     L.need(positive(checker.profiles(graph)[root])==target,'literal compiled rooted profile')
     compiled_profiles+=1;checker.profiles.cache_clear();checker.table.cache_clear()
  bound=999;seen={s:(entry,()) for s,entry in seeds.items()};queue=deque(seen)
  witnesses={};transitions=0
  while queue:
   v=queue.popleft();entry,words=seen[v];h=v[0]+v[2]
   witnesses.setdefault(h,(v,entry,words))
   for w,M in models:
    transitions+=1;u=step(v,M);hu=u[0]+u[2]
    L.need(hu>h, 'strict parent count growth')
    if hu>bound:continue
    if u not in seen:seen[u]=(entry,words+(w,));queue.append(u)
  missing=[h for h in range(3,bound+1,2) if h not in (7,21) and h not in witnesses]
  L.need(7 not in witnesses and 21 not in witnesses,'hostile exclusions')
  literal=[]
  if radius==2:
   L.need(not missing,'declared bounded coverage')
   for h,(v,entry,words) in sorted(witnesses.items()):
    T,r,raw,roots=compiler(entry,words)
    L.need(L.count(T)==h,'actual scalar witness')
    literal.append({'H':h,'profile':v,'entry':entry,'macros':words,'raw_with_final_virtual_birth':raw,'parent_roots':roots,'root':r,'rows':T})
  results.append({'entry_radius':radius,'entry_states':[{'profile':v,'entry':w} for v,w in seeds.items()],
                  'cutoff_H':bound,'states':len(seen),'tested_transitions':transitions,
                  'state_sha256':sha256(json.dumps(sorted(seen),separators=(',',':')).encode()).hexdigest(),
                  'missing_allowed':missing,'compiled_profile_controls':compiled_profiles,'literal_witnesses':literal})
 return {'alphabet':[{'word':w,'matrix':M,'determinant':int(S.Matrix(M).det()),'rank':S.Matrix(M).rank()} for w,M in models],
         'saturations':results,'scope':'finite source grammar with entry radius1/2 and excursion length at most3; no universal target claim'}

def main():
 from fractions import Fraction
 L.need(all(x==0 for x in (J*B*C)[4,:]),'old end-count cancels')
 L.need(F.rank()==4 and all(x>=0 for x in F),'positive full-rank sealing observer')
 denominators=sorted({Fraction(str(x)).denominator for A in As for x in A})
 L.need(denominators==[1,2],'dyadic observer maps')
 result={'scope':'Exact universal root-child factor and amplification certificates; one bounded constructive reset grammar',
         'coordinate_order':['h-s','q-s','s','F-h-q+s'],'sealing_factor':ints(F),
         'sealing_rank':4,'five_signature_birth_rank':int(B.rank()),
         'query_denominators':denominators,'amplification':amplification(),
         'literal_factor_controls':controls(),'grammar':grammar()}
 print(json.dumps(result,indent=2,sort_keys=True))

if __name__=='__main__':main()
