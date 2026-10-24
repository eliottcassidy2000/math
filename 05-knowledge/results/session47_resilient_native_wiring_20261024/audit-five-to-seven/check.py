#!/usr/bin/env python3
"""Independent native5-to7 obstruction audit, with exhaustive 7! canonization."""
from pathlib import Path
from itertools import product,permutations
from collections import Counter
from random import Random
import csv,hashlib,json,subprocess,argparse

HERE=Path(__file__).resolve().parent
PACKET=None

def resolve_packet(producer_dir=None):
 candidates=([Path(producer_dir).expanduser().resolve()] if producer_dir else
             [HERE.parent/'carrier-transfer'/'five-to-seven', HERE.parent/'session47-carrier-transfer'/'five-to-seven'])
 for path in candidates:
  if (path/'normal.json').is_file():return path
 raise FileNotFoundError('Cannot find five-to-seven/normal.json; pass --producer-dir with the five-to-seven packet directory.')
COUNT=0

def need(x,msg):
 global COUNT
 COUNT+=1
 if not x:raise RuntimeError(msg)

def matrix(n,w):
 A=[[0]*n for _ in range(n)];k=0
 for i in range(n):
  for j in range(i+1,n):
   a=1 if j==i+1 else (w>>k)&1
   k+=int(j>i+1);A[i][j]=a;A[j][i]=1-a
 return A

def encode(A):
 c=0;k=0
 for i in range(len(A)):
  for j in range(i+1,len(A)):c|=A[i][j]<<k;k+=1
 return c

def python_canonical(A):
 n=len(A);lo=1<<(n*(n-1)//2);hi=0
 for p in permutations(range(n)):
  c=0
  for i in range(n):
   for j in range(i+1,n):c=(c<<1)|A[p[i]][p[j]]
  lo=min(lo,c);hi=max(hi,c)
 conv=((1<<(n*(n-1)//2))-1)^hi
 return lo,conv,min(lo,conv)

def code_canon(inputs):
 p=subprocess.run([str(HERE/'full_canonical'),'codes'],input=''.join(f'{n} {code}\n' for n,code in inputs),text=True,capture_output=True,check=True)
 return [tuple(map(int,row.split())) for row in p.stdout.splitlines()]

def graph(n,stored):
 path=HERE/f'canonical{n}.tsv'
 subprocess.run([str(HERE/'full_canonical'),str(n),str(path)],check=True,stdout=subprocess.DEVNULL,stderr=subprocess.DEVNULL)
 canon=[tuple(int(r[k]) for k in ('ordinary','converse','merged')) for r in csv.DictReader(path.open(),delimiter='\t')]
 words=1<<((n-1)*(n-2)//2);full=words-1
 need(len(canon)==len(stored['owner'])==words,('all literal words',n))
 ordinary={r[0] for r in canon};classes={r[2] for r in canon};selfconv={r[0] for r in canon if r[0]==r[1]}
 mapping={};reverse={}
 for w,(o,c,merged) in enumerate(canon):
  label=stored['owner'][w]
  if merged in mapping:need(mapping[merged]==label,('canonical class split',n,w))
  if label in reverse:need(reverse[label]==merged,('canonical classes collapsed',n,w))
  mapping[merged]=label;reverse[label]=merged
 need(set(reverse)==set(range(len(classes))) and len(stored['keys'])==len(classes),('complete node labels',n))
 inputs=[];q=n*(n-1)//2
 for key in stored['keys']:
  inputs.append((n,sum(((key>>(q-1-k))&1)<<k for k in range(q))))
 for i,c in enumerate(code_canon(inputs)):need(c[2]==reverse[i],('stored key authenticates node',n,i))
 # Separate Python all-permutation check of the accelerated Gray-code canonizer.
 probes=list(range(words)) if n==5 else sorted(Random(5747).sample(range(words),12))
 for w in probes:need(python_canonical(matrix(n,w))==canon[w],('Python/C++ canonical cross-check',n,w))
 chords=[(i,j) for i in range(n) for j in range(i+2,n)];ix={p:k for k,p in enumerate(chords)}
 colors=Counter();L=[[0]*len(classes) for _ in classes]
 for w in range(words):
  reflection=sum(((w>>k)&1)<<ix[n-1-j,n-1-i] for k,(i,j) in enumerate(chords))
  need(canon[reflection][0]==canon[w][1],('actual ordinary converse under reflection',n,w))
  if (w>>ix[0,n-1])&1:continue
  u,v=sorted((stored['owner'][w],stored['owner'][w^full]));color='blue' if reflection==w else 'black';colors[u,v,color]+=1
  if u!=v:L[u][u]+=1;L[v][v]+=1;L[u][v]-=1;L[v][u]-=1
 pairs=sorted({(i,j) for i,j,c in colors})
 need([[i,j,colors[i,j,'blue']+colors[i,j,'black']] for i,j in pairs]==stored['edges'],('every literal edge weight',n))
 need(L==stored['L'],('complete literal Laplacian',n))
 need(all(colors[i,j,'black']%2==0 for i,j in pairs),('reflected black pair multiplicity',n))
 need(sum(colors.values())==words//2,('complete complement-line mass',n))
 seen={0}
 while True:
  new=seen|{j for i in seen for j,v in enumerate(L[i]) if i!=j and v<0}
  if new==seen:break
  seen=new
 need(len(seen)==len(classes),('connected whole carrier',n))
 summary={'order':n,'words':words,'permutations_per_word':120 if n==5 else 5040,'ordinary_classes':len(ordinary),'self_converse_classes':len(selfconv),'merged_nodes':len(classes),'literal_lines':sum(colors.values()),'blue_lines':sum(v for (i,j,c),v in colors.items() if c=='blue'),'black_lines':sum(v for (i,j,c),v in colors.items() if c=='black'),'blue_loops':sum(v for (i,j,c),v in colors.items() if c=='blue' and i==j),'black_loops':sum(v for (i,j,c),v in colors.items() if c=='black' and i==j),'edge_support':len(pairs),'connected':True,'colored_edges':[[i,j,colors[i,j,'blue'],colors[i,j,'black']] for i,j in pairs]}
 return summary,mapping,canon

def mv(A,x,p):return tuple(sum(a*b for a,b in zip(row,x))%p for row in A)

def binary_rank(rows):
 basis={}
 for r in rows:
  while r:
   j=r.bit_length()-1
   if j in basis:r^=basis[j]
   else:basis[j]=r;break
 return len(basis)

def bitrows(A):return [sum((a&1)<<j for j,a in enumerate(row)) for row in A]

def determinant_mod(A,p):
 # Column pivoting, distinct from the producer's row-pivot echelon routine.
 A=[[a%p for a in row] for row in A];n=len(A);det=1
 for k in range(n):
  j=next((j for j in range(k,n) if A[k][j]),None)
  if j is None:return 0
  if j!=k:
   for row in A:row[k],row[j]=row[j],row[k]
   det=-det
  pivot=A[k][k];det=det*pivot%p;inverse=pow(pivot,-1,p)
  for i in range(k+1,n):
   factor=A[i][k]*inverse%p
   if factor:
    for j in range(k+1,n):A[i][j]=(A[i][j]-factor*A[k][j])%p
   A[i][k]=0
 return det%p

def det_bareiss(A):
 A=[row[:] for row in A];n=len(A);last=1;sign=1
 for k in range(n-1):
  if not A[k][k]:
   j=next(i for i in range(k+1,n) if A[i][k]);A[k],A[j]=A[j],A[k];sign=-sign
  pivot=A[k][k]
  for i in range(k+1,n):
   for j in range(k+1,n):
    v=A[i][j]*pivot-A[i][k]*A[k][j];need(v%last==0,'exact Bareiss division');A[i][j]=v//last
  for i in range(k+1,n):A[i][k]=0
  last=pivot
 return sign*A[-1][-1]

def field_rank(A,p):
 # Column elimination on the transpose supplies an independent rank path.
 rows=[[a%p for a in col] for col in zip(*A)];rank=0
 for j in range(len(rows[0])):
  i=next((i for i in range(rank,len(rows)) if rows[i][j]),None)
  if i is None:continue
  rows[i],rows[rank]=rows[rank],rows[i];inv=pow(rows[rank][j],-1,p);rows[rank]=[a*inv%p for a in rows[rank]]
  for i in range(rank+1,len(rows)):
   q=rows[i][j]
   if q:rows[i]=[(a-q*b)%p for a,b in zip(rows[i],rows[rank])]
  rank+=1
  if rank==len(rows):break
 return rank

def main():
 raw=(PACKET/'normal.json').read_bytes();data=json.loads(raw)
 subprocess.run(['g++','-std=c++17','-O3','-Wall','-Wextra','-pedantic',str(HERE/'full_canonical.cpp'),'-o',str(HERE/'full_canonical')],check=True,capture_output=True)
 old,oldmap,oldcanon=graph(5,data['old']);new,newmap,newcanon=graph(7,data['new'])
 need((old['ordinary_classes'],old['merged_nodes'])==(12,10),'all old tournament classes')
 need((new['ordinary_classes'],new['merged_nodes'])==(456,272),'all new tournament classes')
 # Rebuild every source and target readout from actual adjacency matrices.
 inputs=[];records=[]
 for w in range(64):
  B=matrix(5,w);S=[[0]*7 for _ in range(7)]
  for i in range(5):
   for j in range(5):S[i][j]=B[i][j]
   S[5][i]=1;S[i][6]=1
  S[5][6]=1;path=(5,0,1,2,3,4,6);arcs={frozenset((a,b)) for a,b in zip(path,path[1:])}
  need(all(S[a][b] for a,b in zip(path,path[1:])),'actual suspended Hamiltonian path')
  T=[[0 if i==j else S[i][j] if frozenset((i,j)) in arcs else 1-S[i][j] for j in range(7)] for i in range(7)]
  need([row[:5] for row in T[:5]]==matrix(5,w^63),'literal inner complement')
  inputs.extend(((7,encode(S)),(7,encode(T))));records.append((w,B,S,T))
 answers=code_canon(inputs);readouts=set();source_readouts=set();record_summary=[]
 for (w,B,S,T),source,target in zip(records,answers[::2],answers[1::2]):
  ns=newmap[source[2]];nt=newmap[target[2]];os=data['old']['owner'][w];ot=data['old']['owner'][w^63]
  source_readouts.add((ns,os));readouts.update(((ns,os),(nt,ot)));record_summary.append({'word':w,'old_source':os,'old_target':ot,'new_source':ns,'new_target':nt})
 need(sorted(readouts)==[tuple(r) for r in data['readouts']],'all source and target readouts')
 need(len(readouts)==48,'48 distinct readouts')
 need(all(len({o for n,o in readouts if n==node})==1 for node,_ in readouts),'no endpoint collision')
 need(data['new_root']==0 and data['new_reduced_coordinates']==list(range(1,272)),'grounding includes all new rows')
 L5=data['old']['L'];L7=data['new']['L'];A=[row[1:] for row in L5[1:]];N=[row[1:] for row in L7[1:]]
 coords=list(range(1,272));R=[[int(j==node) for j in coords] for node,owner in sorted(readouts)];S=[[int(j==owner) for j in range(1,10)] for node,owner in sorted(readouts)]
 M=N+R;J=[[0]*9 for _ in N]+S
 old_bank=[x for x in product((0,1),repeat=9) if not any(mv(A,x,2))]
 need(len(old_bank)==8,'exact old mod2 bank has dimension3')
 rows=bitrows(M);rank=binary_rank(rows);success=[];failure=[]
 for x in old_bank:
  rhs=mv(J,x,2);aug=[r|(v<<271) for r,v in zip(rows,rhs)]
  (success if binary_rank(aug)==rank else failure).append(x)
 need(len(success)==len(failure)==4,'exact mod2 extendible subspace dimension2')
 need(sorted(success)==sorted(tuple(x) for x in data['tiny_obstruction']['extendible_old_vectors']),'producer identifies complete successful subspace')
 for test in data['tests']:
  p=test['prime']
  for witness in test['successful_generators']:
   x=witness['old_vector'];y=witness['new_vector'];need(not any(mv(A,x,p)),'successful old vector harmonic');need(mv(M,y,p)==mv(J,x,p),'full producer extension literal check')
  for witness in test['failures']:
   x=witness['old_vector'];z=witness['dual'];need(not any(mv(A,x,p)),'failed old vector harmonic')
   need(all(sum(z[i]*M[i][j] for i in range(len(M)))%p==0 for j in range(271)),'full left dual annihilates constraint columns')
   pair=sum(a*b for a,b in zip(z,mv(J,x,p)))%p;need(pair==witness['pairing'] and pair!=0,'dual pairing detects old direction')
 rowset=(10,41,196,244);fullsum=[sum(L7[i][j] for i in rowset)%2 for j in range(272)]
 need(fullsum==[int(j in (10,41)) for j in range(272)],'ungrounded four-row native relation')
 h=tuple(int(i in (3,4,6,7)) for i in range(10));need(not any(mv(L5,h,2)) and h[0]==0,'actual old full harmonic')
 need(h[3]!=h[9],'old character distinguishes required source labels')
 for word,oldnode,newnode in ((55,3,10),(0,9,41)):
  rec=record_summary[word];need(rec['old_source']==oldnode and rec['new_source']==newnode,'literal source-only obstruction')
  need((newnode,oldnode) in source_readouts,'source-only required readout present')
 need(data['old']['keys'][3]==16 and data['old']['keys'][9]==258 and data['new']['keys'][10]==256 and data['new']['keys'][41]==8208,'verified canonical source keys')
 olddet=det_bareiss(A);need(olddet==41096==8*11*467,'complete old determinant')
 odd=[]
 for p,expected in ((11,8),(467,397)):
  rankold=field_rank(A,p);detnew=determinant_mod(N,p)
  need(rankold==8,'one old odd-prime harmonic direction');need(detnew==expected,'new grounded determinant nonzero at odd prime')
  odd.append({'prime':p,'old_rank':rankold,'old_nullity':9-rankold,'new_rank':271,'new_determinant_mod_prime':detnew})
 source_constraints=N+[[int(j==node) for j in coords] for node,owner in sorted(source_readouts)]
 source_rhs=tuple([0]*271+[h[owner] for node,owner in sorted(source_readouts)])
 sr=bitrows(source_constraints);need(binary_rank([r|(v<<271) for r,v in zip(sr,source_rhs)])>binary_rank(sr),'source-only full carrier extension fails')
 return {'status':'PASS: genuine complete native5-to7 source-suspension failure; not8-to10','source_receipt_sha256':hashlib.sha256(raw).hexdigest(),'old_graph':old,'new_graph':new,'source_records':record_summary,'readouts':sorted(readouts),'source_only_readouts':sorted(source_readouts),'old_determinant':olddet,'old_mod2_bank':old_bank,'extendible_old_mod2_bank':success,'nonextendible_old_mod2_bank':failure,'new_full_laplacian_row_relation':{'rows':rowset,'sum_support_mod2':[10,41]},'old_harmonic_support':[3,4,6,7],'old_incompatible_source_nodes':[3,9],'source_readout_witnesses':[record_summary[55],record_summary[0]],'odd_prime_controls':odd,'checks':COUNT}

if __name__=='__main__':
 parser=argparse.ArgumentParser()
 parser.add_argument('--producer-dir',type=Path,help='Directory containing the five-to-seven normal.json receipt; defaults to either supported sibling publication layout.')
 args=parser.parse_args();PACKET=resolve_packet(args.producer_dir)
 print(json.dumps(main(),sort_keys=True,indent=2))
