"""Build exact input, check the representation, and run the exhaustive certificate."""
import json,subprocess,sys,tempfile,itertools,random
from pathlib import Path
HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'ribbon'))
from transfer import automaton,advance,graph,held_karp,count
D=json.loads((HERE.parent/'ribbon-periodic'/'representation.json').read_text())
B=D['basis'];Rs=D['matrices'];f=D['output'];states,mat=automaton(3);N=len(states);r=len(f)
assert N==119 and r==62
for i,st in enumerate(states):
 assert sum(B[i][j]*f[j]for j in range(r))==int(len(st)==1)
 for a in range(8):
  for j in range(r):
   assert sum(c*B[k][j]for k,c in mat[i][a])==sum(B[i][k]*Rs[a][k][j]for k in range(r))
# The exact cap-four automaton recognizes all words with H<=3.
final=[i for i,s in enumerate(states)if len(s)==1];v=(1,)+(0,)*(N-1);vectors=[v];index={v:0};dfa=[]
for v in vectors:
 row=[]
 for a in range(8):
  z=advance(v,mat,a,4)
  if sum(z[i]for i in final)>3:row.append(-1);continue
  if z not in index:index[z]=len(vectors);vectors.append(z)
  row.append(index[z])
 dfa.append(row)
assert len(dfa)==288
# Independent literal four-vertex base fact used by the low-cell lemma.
strong4=0
for bits in range(64):
 a=[[False]*4 for _ in range(4)];pairs=list(itertools.combinations(range(4),2))
 for k,(i,j)in enumerate(pairs):a[i][j]=bool(bits>>k&1);a[j][i]=not a[i][j]
 reach=[{i}|{j for j in range(4)if a[i][j]}for i in range(4)]
 for k in range(4):
  for i in range(4):
   if k in reach[i]:reach[i]|=reach[k]
 if all(len(s)==4 for s in reach):assert held_karp(a)==5;strong4+=1
assert strong4==24
lines=['62 8']+[' '.join(map(str,row))for A in Rs for row in A]
lines +=[str(len(dfa))]+[' '.join(map(str,row))for row in dfa]
results=[]
with tempfile.TemporaryDirectory(prefix='s197-small-affine-')as tmp:
 exe=Path(tmp)/'verify'
 subprocess.run(['g++','-std=c++17','-O3',str(HERE/'verify.cpp'),'-o',str(exe)],check=True)
 for multiplier in [3,4]:
  run=subprocess.run([str(exe),str(multiplier)],input='\n'.join(lines)+'\n',text=True,capture_output=True,check=True)
  result=json.loads(run.stdout)
  assert result['multiplier']==multiplier
  assert [x['classes']for x in result['short']]==[8,36,176,1044,6560,43800]
  assert [x['low_words']for x in result['long']]==[62976,199616,604160,1768448,5147392]
  assert [x['active_words']for x in result['long']]==[29952,39872,29440,24576,9216]
  assert [x['filtered_classes']for x in result['long']]==[417,1091,921,675,329]
  assert result['unresolved']==0
  results.append(result)
payload={'multipliers':results,'strong_four_vertex_graphs_checked':strong4}
(HERE/'VALIDATION.json').write_text(json.dumps(payload,indent=2)+'\n');print(json.dumps(payload,indent=2))
