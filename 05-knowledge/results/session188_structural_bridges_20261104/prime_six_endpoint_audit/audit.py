"""Independent complete integer word-completion census and literal physical audit."""
from pathlib import Path
from itertools import combinations,permutations,product
from array import array
from collections import Counter
import hashlib,json,subprocess,tempfile
ROOT=Path(__file__).resolve().parent;SRC=ROOT.parent/'prime_six_reconnaissance';PAIRS=list(combinations(range(6),2))
def need(x,m):
 if not x:raise RuntimeError(m)
def rows(code):
 R=[0]*6
 for b,(i,j)in enumerate(PAIRS):
  if code>>b&1:R[i]|=1<<j
  else:R[j]|=1<<i
 return R
def op(R):return [((1<<len(R))-1)^v^(1<<i)for i,v in enumerate(R)]
def census():
 S=array('H',[0])*(32768*6);E=array('H',[0])*(32768*6);ix={e:k for k,e in enumerate(PAIRS)};completions=0
 for w in permutations(range(6)):
  fixed=0;ones=0
  for u,v in zip(w,w[1:]):
   bit=1<<ix[min(u,v),max(u,v)];fixed|=bit
   if u<v:ones|=bit
  free=32767^fixed;M=free
  while True:
   code=M|ones;S[6*code+w[0]]+=1;E[6*code+w[-1]]+=1;completions+=1
   if not M:break
   M=(M-1)&free
 need(completions==720*1024,'every forced-edge completion once per Hamilton word')
 return S,E,completions
def module_witness(R):
 n=len(R);full=(1<<n)-1;incoming=[sum(1<<w for w in range(n)if R[w]>>v&1)for v in range(n)];D=[[((incoming[u]^incoming[v])&~((1<<u)|(1<<v)))for v in range(n)]for u in range(n)]
 for M in range(1,full):
  if M.bit_count()<2:continue
  u=(M&-M).bit_length()-1;outside=full^M
  if all(not(D[u][v]&outside)for v in range(n)if M>>v&1):return M
 return 0
def relabel_code(Q,p):
 R=[0]*6
 for i in range(6):
  for j in range(6):
   if Q[i]>>j&1:R[p[i]]|=1<<p[j]
 return sum((bool(R[i]>>j&1))<<b for b,(i,j)in enumerate(PAIRS))
def tt(n):return [sum(1<<j for j in range(i+1,n))for i in range(n)]
def blow(Q,B):
 off=[];n=0
 for b in B:off.append(n);n+=len(b)
 out=[0]*n
 for i,b in enumerate(B):
  for u in range(len(b)):
   for j,c in enumerate(B):
    for v in range(len(c)):
     if(i==j and b[u]>>v&1)or(i!=j and Q[i]>>j&1):out[off[i]+u]|=1<<(off[j]+v)
 return out,off
def expand(mask,B,off):return sum(1<<(off[i]+v)for i,b in enumerate(B)for v in range(len(b))if mask>>i&1)
def aug(R,a,c):
 n=len(R);out=list(R)+[1<<(n+1),0]
 for i in range(n):
  for v,M in [(n,a),(n+1,c)]:
   if M>>i&1:out[i]|=1<<v
   else:out[v]|=1<<i
 return out
def outer(R,a,c):
 out,off=blow([2,4,1],[R,[0],[0]]);incoming=1<<off[2];return aug(out,a|incoming,c|incoming)
def subs(M):
 S=M
 while True:
  yield S
  if not S:break
  S=(S-1)&M
def main():
 r=json.loads((SRC/'normal.json').read_text());need((SRC/'normal.json').read_bytes()==(SRC/'optimized.json').read_bytes(),'paired producer bytes');need(all(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h for f,h in r['source_pins'].items()),'frozen source pins')
 table=[list(map(int,s.split()))for s in(SRC/'endpoint_table.tsv').read_text().splitlines()];need(len(table)==32768 and all(len(v)==15 and v[0]==i for i,v in enumerate(table)),'exact complete table schema and coverage')
 S,E,words=census();primecodes=set();hist=Counter();phist=Counter();endpointbytes=bytearray()
 for code,row in enumerate(table):
  ss=list(S[6*code:6*code+6]);ee=list(E[6*code:6*code+6]);p=sum(((a+b)&1)<<i for i,(a,b)in enumerate(zip(ss,ee)));need(ss==row[3:9]and ee==row[9:15]and p==row[2]and p,'every exact S/E coordinate and nonzero endpoint mask');w=module_witness(rows(code));need(w==row[1],'independent distinguishing-pair module test');hist[p.bit_count()]+=1
  if w==0:primecodes.add(code);phist[p.bit_count()]+=1
  endpointbytes.extend(ss+ee+[p])
 need(hist==Counter({2:13920,4:12480,6:6368})and phist==Counter({2:720,4:5760,6:3840})and len(primecodes)==10320,'all finite histograms')
 used=set();orbit_sizes=[]
 for o in r['prime_orbits']:
  Q=rows(o['code']);need(Q==o['rows'],'orbit physical code');orb={relabel_code(Q,p)for p in permutations(range(6))};need(orb==set(o['codes'])and not(used&orb),'every full orbit code and disjointness');used|=orb;orbit_sizes.append(len(orb));need(o['S']==table[o['code']][3:9]and o['E']==table[o['code']][9:15]and o['p']==table[o['code']][2],'orbit endpoint labels')
 need(used==primecodes and sorted(orbit_sizes)==[240]+[720]*14,'complete fifteen prime orbits')
 boundary=r['prime8_boundary'];Q8=boundary['rows'];need(Q8==[238,4,240,118,195,146,34,74],'specified prime8 boundary');ss=[0]*8;ee=[0]*8
 for w in permutations(range(8)):
  if all(Q8[u]>>v&1 for u,v in zip(w,w[1:])):ss[w[0]]+=1;ee[w[-1]]+=1
 need(ss==boundary['S']and ee==boundary['E']and all((a+b)%2==0 for a,b in zip(ss,ee))and module_witness(Q8)==0,'prime8 scope boundary')
 graphs=[];expected=[];kinds=[]
 def add(G,values,kind):graphs.append(G);expected.append(values);kinds.append(kind)
 need({(x['code'],x['name'])for x in r['actual_profile_controls']}=={(o['code'],name)for o in r['prime_orbits']for name in ['single','pairs','mixed','mixed_ambient']}and len(r['actual_profile_controls'])==60,'all sixty stated physical comparisons')
 for cr in r['actual_profile_controls']:
  Q=rows(cr['code']);B=cr['child_rows'];i=cr['anchor'];need(table[cr['code']][2]>>i&1,'actual endpoint anchor');need(B==([tt(1)]*6 if cr['name']=='single' else [tt(2)]*6 if cr['name']=='pairs' else[tt(1),tt(2),[2,4,1],tt(1),tt(2),tt(1)]),'specified actual child models');need([(x-y)%4 for x,y in zip(*cr['profiles'])]==[2,0,0],'full first-three physical contrast')
  for side,QQ in enumerate([Q,op(Q)]):
   core,off=blow(QQ,B);a=expand(1<<i,B,off);G=outer(core,a,0)if cr['outer_cyclic']else aug(core,a,0);need(G==cr['graphs'][side]and cr['core_masks'][side]==[a,0],'profile graph labels and root incidences');add(G,cr['profiles'][side],cr['name'])
 Q=r['prime_orbits'][0]['rows'];bank=r['complete4096_named_response_bank'];need([v[:2]for v in bank]==[[a,c]for a,c in product(range(64),repeat=2)],'complete physical4096 bank')
 for a,c,t,u in bank:add(aug(Q,a,c),[t],'full_bank');add(aug(op(Q),a,c),[u],'full_bank')
 B=[tt(2)]*6;T,off=blow(Q,B);p=r['prime_orbits'][0]['p'];anchor=(p&-p).bit_length()-1;a=expand(1<<anchor,B,off);add(aug(T,a,0),[r['root_pair_family_source_count']],'family_source')
 need([v[0]for v in r['all64_original_root_pair_families']]==list(range(64)),'all64 root/pair subsets')
 for bits,z in r['all64_original_root_pair_families']:
  need((r['root_pair_family_source_count']-z)%4==2,'every pair-family separator');U=op(T)
  for i in range(6):
   if bits>>i&1:U[2*i]^=1<<(2*i+1);U[2*i+1]^=1<<(2*i)
  add(aug(U,a,0),[z],'pair_family')
 B=[tt(3)]+[tt(1)]*5;T,off=blow(Q,B);U,_=blow(op(Q),B);V=op(T)
 for M in [3,6]:
  need(all(T[v]&M in(0,M)for v in range(len(T))if not M>>v&1),'original crossing interval module')
  for u,v in combinations(range(len(T)),2):
   if M>>u&1 and M>>v&1:V[u]^=1<<v;V[v]^=1<<u
 aa=expand(1<<anchor,B,off);MM=(1<<len(T))-1;need([v[:2]for v in r['crossing_ambient_queries']]==[[aa,0],[MM,MM^aa]],'both physical dual crossing queries');baseline=[];targetgaps=[]
 for aa,cc,values in r['crossing_ambient_queries']:
  for G,z in zip([T,U,V],values):add(outer(G,aa,cc),z,'crossing_ambient')
  baseline.append(values[0][0]-values[1][0]);targetgaps.append(values[0][0]-values[2][0])
 need(baseline==r['baseline_gaps']and baseline[0]==-baseline[1]and targetgaps[r['sign_selected_query_index']]<0,'exact negative-dual family scope')
 with tempfile.TemporaryDirectory()as td:
  exe=Path(td)/'count';subprocess.run(['g++','-std=c++17','-O2','-Wall','-Wextra','-Werror',str(ROOT/'count.cpp'),'-o',str(exe)],check=True);wire=str(len(graphs))+'\n'+''.join(str(len(G))+' '+' '.join(map(str,G))+'\n'for G in graphs);values=[list(map(int,s.split()))for s in subprocess.check_output([str(exe)],input=wire.encode()).decode().splitlines()]
 need(len(values)==len(graphs)==r['physical_graphs'],'every requested physical graph');need(all(v[:len(e)]==e for v,e in zip(values,expected)),'all exact physical values and first-three cover coordinates')
 d={(a,c):[t,u]for a,c,t,u in bank};half=[M for M in range(64)if M.bit_count()<=2];coeff=[];hist=Counter()
 for a,c in product(half,repeat=2):
  z=[sum((-1)**(a.bit_count()+c.bit_count()-x.bit_count()-y.bit_count())*d[x,y][side]for x in subs(a)for y in subs(c))for side in range(2)];coeff.append([a,c,*z]);hist[a.bit_count()+c.bit_count()]+=z[0]!=z[1]
 need(coeff==r['complete484_named_coefficient_bank'],'every cubic/quartic coefficient and zero')
 for a,c,t,u in bank:need([sum(v[side]for x,y,*v in coeff if x&a==x and y&c==y)for side in range(2)]==[t,u],'all4096 actual queries reconstructed from the full degree basis')
 need(all(hashlib.sha256((SRC/f).read_bytes()).hexdigest()==h for f,h in r['source_pins'].items()),'source unchanged during audit')
 return {'status':'PASS','source_pins':r['source_pins'],'auditor_pins':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest()for f in ['audit.py','count.cpp']},'complete_labelled_codes':32768,'literal_word_completions':words,'full_endpoint_data_sha256':hashlib.sha256(endpointbytes).hexdigest(),'endpoint_support_histogram':{str(k):v for k,v in Counter(row[2].bit_count()for row in table).items()},'prime_count':len(primecodes),'prime_endpoint_support_histogram':dict(phist),'prime_orbit_sizes':orbit_sizes,'physical_graph_occurrences':len(graphs),'distinct_physical_graphs':len({tuple(G)for G in graphs}),'verified_physical_coordinates':sum(len(e)for e in expected),'all_physical_values_sha256':hashlib.sha256(json.dumps(values,separators=(',',':')).encode()).hexdigest(),'physical_control_kinds':dict(Counter(kinds)),'full_named_coefficients':len(coeff),'difference_histogram_by_degree':dict(hist),'prime8_boundary':boundary,'scope':'Complete six-name endpoint certificate and exact fixed-ambient physical controls; no claim for every even quotient order.'}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
