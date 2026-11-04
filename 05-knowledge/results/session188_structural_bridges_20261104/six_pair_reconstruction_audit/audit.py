"""Independent six-pair signatures and exact predecessor-DP response audit."""
from pathlib import Path
from itertools import combinations,permutations
import struct,json,hashlib,subprocess,tempfile,sys
HERE=Path(__file__).resolve().parent;BASE=HERE.parent;SRC=BASE/'six_pair_reconstruction'
def need(x,msg):
 if not x:raise RuntimeError(msg)
def graph(n,code):
 R=[0]*n
 for k,(i,j) in enumerate(combinations(range(n),2)):
  u,v=(i,j) if code>>k&1 else(j,i);R[u]|=1<<v
 return tuple(R)
def op(R):return tuple(((1<<len(R))-1)^row^(1<<i) for i,row in enumerate(R))
def induce(R,S):return tuple(sum(1<<j for j,v in enumerate(S) if R[u]>>v&1) for u in S)
def feature5(R):
 n=len(R);v=[0]*n
 for size in (2,4):
  for S in combinations(range(n),size):
   for p in permutations(S):
    if all(R[i]>>j&1 for i,j in zip(p,p[1:])):v[p[0]]^=1;v[p[-1]]^=1
 return v
def inflate(code,flips):
 Q=graph(6,code);R=[0]*12
 for i in range(12):
  for j in range(i+1,12):
   yes=(not(flips>>(i//2)&1)) if i//2==j//2 else bool(Q[i//2]>>(j//2)&1)
   R[i if yes else j]|=1<<(j if yes else i)
 return tuple(R)
def attach(R,a,c):
 n=len(R);G=list(R)+[1<<(n+1),0]
 for i in range(n):
  for v,m in((n,a),(n+1,c)):
   if m>>i&1:G[i]|=1<<v
   else:G[v]|=1<<i
 return tuple(G)
def expand(m):return sum(3<<(2*i) for i in range(6) if m>>i&1)
def ambient(R):
 n=len(R);G=list(R)+[1<<(n+1),0]
 for i in range(n):G[i]|=1<<n;G[n+1]|=1<<i
 return tuple(G)
def signature(values):
 out=0;k=0
 for i in range(6):
  for j in range(6):
   if i==j:continue
   total=0
   for aa in (0,1<<i,1<<j,(1<<i)|(1<<j)):
    for c in(0,1<<j):total+=(-1)**(3-aa.bit_count()-c.bit_count())*(values[aa,c]-values[63^c,63^aa])
   need(total%4==0,'shared contrast divisibility');out|=((total//4)&1)<<k;k+=1
 return out

def main():
 sources=[SRC/n for n in ('PROOF.md','verify.py','literal.cpp')]+[BASE/d/'PROOF.md' for d in ('even_pair_refinement','mod8_transfer')];pins={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources};data=json.loads((SRC/'literal.json').read_text());raw=(SRC/'h_vectors.bin').read_bytes();vectors=struct.unpack('<32768I',raw);classes={}
 for code in range(1024):
  Q=graph(5,code);need(feature5(Q)==[r.bit_count()%2 for r in Q],'literal permutation endpoint lemma')
 for code in range(32768):
  Q=graph(6,code);want=sum(((Q[i].bit_count()+(Q[i]>>j&1))&1)<<k for k,(i,j) in enumerate(( (i,j) for i in range(6) for j in range(6) if i!=j)))
  need(vectors[code]==want and want!=0,'complete retained signature bank');classes.setdefault(want,[]).append(code)
 need(len(classes)==16384 and all(len(v)==2 and v[0]^v[1]==32767 for v in classes.values()),'exact converse signature classes')
 low={(a,0) for a in range(64) if a.bit_count()<=2}
 low|={(a,1<<j) for j in range(6) for a in range(64) if a.bit_count()<=1 or(a.bit_count()==2 and a>>j&1)}
 battery=sorted(low|{(63^c,63^a) for a,c in low});need(len(low)==94 and len(battery)==188 and battery==list(map(tuple,data['coarse_battery'])),'independent fixed battery')
 graphs=[];idx={}
 def add(R):
  if R not in idx:idx[R]=len(graphs);graphs.append(R)
  return idx[R]
 records=[]
 for rec in data['states']:
  ids=[add(attach(inflate(rec['quotient_code'],rec['pair_flip_mask']),expand(a),expand(c))) for a,c in battery];records.append((rec['quotient_code'],ids,rec['exact_values']))
 pairids=[];Q=graph(6,81);masks=[(1<<(2*i))|sum(3<<(2*j) for j in range(6) if Q[j]>>i&1) for i in range(6)];need(masks==data['pair_decoder_physical_masks'],'independent physical decoder masks')
 for flips in range(64):pairids.append([add(attach(inflate(81,flips),m,m)) for m in masks])
 extra=[]
 for code in (17493,20261):extra.append((code,[add(attach(inflate(code,42),expand(a),expand(c))) for a,c in battery]))
 ambids=[add(ambient(attach(inflate(81,0),expand(a),expand(c)))) for a,c in battery]
 control=data['strong_ambient_control'];ctrlids=[add(ambient(attach(inflate(code,0),expand(control['coarse_a']),expand(control['coarse_c'])))) for code in (81,32767-81)]
 with tempfile.TemporaryDirectory() as td:
  binary=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(HERE/'count.cpp'),'-o',str(binary)],check=True)
  wire=str(len(graphs))+'\n'+'\n'.join(str(len(G))+' '+' '.join(map(str,G)) for G in graphs)+'\n';values=list(map(int,subprocess.check_output([str(binary)],input=wire.encode()).split()));need(len(values)==len(graphs),'all exact values')
 for code,ids,want in records:
  vs=[values[i] for i in ids];need(vs==want,'independent full physical source record');need(signature(dict(zip(battery,vs)))==vectors[code],'independent full signature recovery')
 for flips,ids in enumerate(pairids):need([values[i] for i in ids]==data['pair_decoder_exact_values'][flips],'independent all64 pair decoder states')
 for code,ids in extra:need(signature(dict(zip(battery,[values[i] for i in ids])))==vectors[code],'additional quotient signatures')
 need(signature(dict(zip(battery,[values[i] for i in ambids])))==vectors[81],'complete strong ambient signature and reflection convention')
 need([values[i] for i in ctrlids]==[control['ambient_source'],control['ambient_target']],'retained ambient witness')
 need(pins=={str(p.relative_to(BASE)):hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},'source frozen')
 print(json.dumps({'status':'PASS','independent_literal_five_name_lemmas':1024,'exact_six_name_signatures':32768,'converse_classes':16384,'independent_exact_physical_tournaments':len(graphs),'retained_coarse_states_checked':len(records),'additional_coarse_states':len(extra),'ambient_full_battery_queries':len(ambids),'pair_states_checked':64,'source_sha256':pins,'data_sha256':{n:hashlib.sha256((SRC/n).read_bytes()).hexdigest() for n in ('literal.json','h_vectors.bin')},'audit_sha256':{n:hashlib.sha256((HERE/n).read_bytes()).hexdigest() for n in ('audit.py','count.cpp')}},indent=2,sort_keys=True))
if __name__=='__main__':main()
