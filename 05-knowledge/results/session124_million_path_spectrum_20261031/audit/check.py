"""Independent seed witnesses, product closure and adjacency preparation; no source imports."""
import argparse,itertools,json,math,pathlib
base=pathlib.Path(__file__).resolve().parent.parent
parser=argparse.ArgumentParser()
parser.add_argument('--work-dir',type=pathlib.Path,required=True,help='independent.cpp output directory')
args=parser.parse_args();dst=args.work_dir;src=base/'constructors'
def require(b,msg):
 if not b: raise RuntimeError(msg)
def seed(n,code,cut):
 A=[[False]*n for _ in range(n)]
 for k,(i,j) in enumerate(itertools.combinations(range(n),2)):
  A[i][j]=bool(code&(1<<k));A[j][i]=not A[i][j]
 return A,[bool(cut&(1<<i)) for i in range(n)]
def literal(A,bits):
 n=len(A); h=f=q=t=0
 for p in itertools.permutations(range(n)):
  if all(A[u][v] for u,v in zip(p,p[1:])):
   h+=1;t+=not bits[p[-1]]
  for s in range(n+1):
   P,Q=p[:s],p[s:]
   if not all(A[u][v] for part in (P,Q) for u,v in zip(part,part[1:])):continue
   right=not Q or bits[Q[0]]
   q+=bool(Q) and right
   f+=(not P or not bits[P[-1]]) and right
 return h,f,q,t
roots=json.loads((src/'seed_roots.json').read_text())
independent={tuple(map(int,l.split())) for l in (dst/'seeds.tsv').read_text().splitlines()}
require({tuple(r['state']) for r in roots}==independent,'seed universe mismatch')
for r in roots:require(literal(*seed(*r['witness']))==tuple(r['state']),('bad seed witness',r))
print('all 273 stored actual seed witnesses and independently enumerated universe PASS')
support=(dst/'support.bin').read_bytes()
require(len(support)==1000001 and set(support)<={0,1},'invalid support bitmap')
require(sum(support)==499980 and not support[7] and not support[21],'wrong direct support')
factors={34869:[3,59,197],85885:[5,89,193],116587:[23,37,137],153451:[173,887],245923:[31,7933],535405:[5,13,8237],595543:[13,61,751],609873:[3,11,18481],618723:[49,9,23,61],767135:[5,153427],782103:[3,260701],824119:[149,5531],826373:[353,2341],868367:[157,5531]}
for m,fs in factors.items():require(math.prod(fs)==m and all(support[f] for f in fs),('bad factor repair',m))
repairs={103061,594793,774073,949777}
miss={i for i in range(1,1000000,2) if not support[i] and i not in (7,21)}
require(miss==set(factors)|repairs,'miss partition')
print('14 SCC product repairs checked against independent support PASS')
for p in (14723,594793,774073,949777):require(all(p%d for d in range(2,math.isqrt(p)+1)),('not prime',p))
rows=[json.loads(l) for l in (src/'block_target_witnesses.jsonl').read_text().splitlines()]
require({r['target'] for r in rows}==repairs,'repair set')
for r in rows:
 A,bits=seed(*r['seed'])
 for block in r['blocks']:
  size=3 if block[0]=='C' else 1;old=len(A)
  A=[row+[True]*size for row in A]+[[False]*(old+size) for _ in range(size)]
  if size==3:
   for i in range(3):A[old+i][old+(i+1)%3]=True
   bits.extend(i<int(block[1]) for i in range(3))
  else:bits.append(block=='1')
 n=len(A);A=[row+[not b] for row,b in zip(A,bits)]+[bits+[False]]
 require(all(not A[i][i] for i in range(n+1)) and all(A[i][j]!=A[j][i] for i in range(n+1) for j in range(i+1,n+1)),'not tournament')
 (dst/f"adjacency_{r['target']}.txt").write_text(str(n+1)+' '+str(r['target'])+'\n'+'\n'.join(''.join('1' if x else '0' for x in row) for row in A)+'\n')
 require((dst/f"adjacency_{r['target']}.txt").read_bytes()==(base/'audit'/f"adjacency_{r['target']}.txt").read_bytes(),'stored adjacency differs')
print('four literal full adjacency matrices prepared; independent DP pending')
