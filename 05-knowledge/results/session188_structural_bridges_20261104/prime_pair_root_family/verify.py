"""Full literal original-module XOR batteries for the marked-pair root theorem."""
from pathlib import Path
from itertools import product
import sys,hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;BASE=ROOT.parent
sys.path.insert(0,str(BASE/'arbitrary_third_root_family'))
from verify import need,blow,extend,module,flip

def quotient(r,kind):
 n=2*r+1
 if kind=='W':return tuple(sum(1<<j for j in range(n) if (i<2*r and j<2*r and i<j) or(i==2*r and j<2*r and j%2==0) or(i<2*r and j==2*r and i%2==1)) for i in range(n))
 Q=[sum(1<<((i+d)%n) for d in range(1,r+1)) for i in range(n)]
 if kind=='U':
  for i in range(r):
   for j in range(i+1,r):Q[i]^=1<<j;Q[j]^=1<<i
 return tuple(Q)
def makecase(label,Q,children,z,I,J):
 T=blow(Q,children);n=len(T);full=(1<<n)-1;offs=[sum(map(len,children[:i])) for i in range(len(children)+1)];blocks=[((1<<len(C))-1)<<offs[i] for i,C in enumerate(children)];p,q=offs[z],offs[z]+1
 need(children[z]==(2,0),'physical root pair');X=sum(blocks[i] for i in I);Y=sum(blocks[i] for i in J)
 quartet=[]
 for bits in product((0,1),repeat=4):quartet.append(((bits[0]<<q)|(X if bits[1] else 0),(bits[2]<<p)|(Y if bits[3] else 0),(-1)**(4-sum(bits))))
 compressed=[B for i,B in enumerate(blocks) if i!=z]+[1<<q];d=sum(bool(T[p]&B) for B in compressed);anchor=[]
 for B in compressed:anchor.append((B|((1<<p) if T[p]&B else 0),1<<p))
 if d%2:anchor.append((1<<p,1<<p))
 if (len(compressed)+d)%2:anchor.append((0,1<<p))
 anchor.append((1<<p,0));need(len(set(anchor))==len(anchor),'anchor distinct masks')
 queries=list(dict.fromkeys([(a,c) for a,c,s in quartet]+anchor));need(len(queries)<=len(Q)+19,'fixed query bound');mods=[M for M in range(1,full) if M.bit_count()>1 and module(T,M)];need(all(any(M&B==M for B in blocks) for M in mods),'proper modules remain inside children')
 cases=[]
 for bits in range(1<<len(mods)):
  U=flip(T,full);selected=[]
  for i,M in enumerate(mods):
   if bits>>i&1:U=flip(U,M);selected.append(M)
  cases.append({'selected_proper_modules':selected,'pair_selected':blocks[z] in selected,'target_rows':U})
 return {'label':label,'quotient':Q,'children':children,'root_pair_child':z,'marked_exterior_child_sets':[I,J],'source_rows':T,'proper_modules':mods,'queries':queries,'quartic':quartet,'anchor':anchor,'families':cases}
def main():
 specs=[]
 for kind in ('T','U','W'):
  Q=quotient(2,kind);I,J=([2],[3]) if kind in('T','U') else([1],[3]);z=0 if kind in('T','U') else 2
  specs.append(makecase(kind+'5_all_pairs',Q,[(2,0)]*5,z,I,J))
 specs.append(makecase('U5_crossing_TT3',quotient(2,'U'),[(2,0),(6,4,0),(2,0),(2,0),(2,0)],0,[2],[3]))
 specs.append(makecase('W5_cyclic_child',quotient(2,'W'),[(2,4,1),(2,0),(2,0),(2,0),(2,0)],2,[1],[3]))
 specs.append(makecase('T7_local_size_boundary',quotient(3,'T'),[(2,0),(0,),(0,),(2,0),(0,),(0,),(0,)],0,[3],[4]))
 graphs=[];index={}
 def add(G):
  G=tuple(G)
  if G not in index:index[G]=len(graphs);graphs.append(G)
  return index[G]
 for rec in specs:
  queries=rec['queries'];rec['baseline_ids']=[add(extend(rec['source_rows'],a,c)) for a,c in queries]
  for f in rec['families']:f['value_ids']=[add(extend(f['target_rows'],a,c)) for a,c in queries]
 with tempfile.TemporaryDirectory() as td:
  binary=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(binary)],check=True)
  wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(binary)],input=wire.encode()).split()));need(len(values)==len(graphs),'all literal counts returned')
 totals={'families':0,'scalar_strict':0,'equal_pair_selected':0,'equal_pair_not_selected':0}
 for rec in specs:
  base=[values[i] for i in rec.pop('baseline_ids')];rec['baseline_values']=base;lookup={tuple(q):i for i,q in enumerate(rec['queries'])};scalar=lookup[(0,0)]
  for f in rec['families']:
   vals=[values[i] for i in f.pop('value_ids')];f['values']=vals;need(vals[scalar]>=base[scalar],'scalar monotonicity');need(vals!=base,'root-family collision');totals['families']+=1
   if vals[scalar]>base[scalar]:totals['scalar_strict']+=1;f['branch']='scalar_strict'
   elif f['pair_selected']:
    gap=sum(s*(base[lookup[a,c]]-vals[lookup[a,c]]) for a,c,s in rec['quartic']);need(gap>0,'selected pair positive quartic');f['branch']='equal_pair_selected';f['quartic_gap']=gap;totals['equal_pair_selected']+=1
   else:
    residues=[(base[lookup[a,c]]-vals[lookup[a,c]])%4 for a,c in rec['anchor']];need(2 in residues,'unselected pair uniform anchor');f['branch']='equal_pair_not_selected';f['anchor_residues_mod4']=residues;totals['equal_pair_not_selected']+=1
 paths=[ROOT/x for x in ('PROOF.md','verify.py','count.cpp')]+[BASE/d/'PROOF.md' for d in ('module_decomposition','marked_pair_quotient','ordered_module_skew','critical_quotient_synthesis')]+[BASE/'returns_next'/'UNIFORM.md',BASE/'arbitrary_third_root_family'/'verify.py']
 print(json.dumps({'status':'PASS','totals':totals,'distinct_literal_augmented_tournaments':len(graphs),'cases':specs,'source_sha256':{str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths}},indent=2,sort_keys=True))
if __name__=='__main__':main()
