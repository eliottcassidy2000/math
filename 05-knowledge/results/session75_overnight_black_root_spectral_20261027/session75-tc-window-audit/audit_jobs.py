"""Independent Cartesian-product reconstruction of each group and all orbitals.
Uses explicit fail checks, not assert; identical under Python -O.
"""
import itertools as it, pathlib, re, json, sys, hashlib
ROOT=pathlib.Path(__file__).resolve().parent

def need(ok, message):
 if not ok: raise RuntimeError(message)

def f(n):
 k,r=divmod(n,3);return (3*5**(k-1),5**k,9*5**(k-1))[r]

def expected(ns):
 ans={}
 for n in ns:
  for r,p in ((0,0),(1,0),(2,0),(0,1)):
   for c in range(n//3+1):
    k=n-5*r-7*p-3*c
    if k<0:continue
    size=3**c*5**r*21**p
    if 21*size<f(n):continue
    prefix=f'n{n}_K{k}_C{c}_R{r}_P{p}_'
    ans[prefix+'P.job']=(n,k,c,r,p,None)
    if 21*size<3*f(n):continue
    for coeff in it.product(range(3),repeat=p+c):
     if not any(coeff) or next(x for x in coeff if x)!=1:continue
     # The only nontransitive kernel projection is singleton support at a C3 coordinate.
     if sum(x!=0 for x in coeff)==1 and next(i for i,x in enumerate(coeff) if x)>=p:continue
     ans[prefix+'ker'+''.join(map(str,coeff))+'.job']=(n,k,c,r,p,coeff)
 return ans

def group(n,k,c,r,p,coeff):
 blocks=[];choices=[];offs=0
 for size,num,kind in ((7,p,'P'),(5,r,'R'),(3,c,'C'),(1,k,'K')):
  for _ in range(num):
   block=list(range(offs,offs+size));blocks.append(block);offs+=size
   if kind=='P':choices.append([(tuple(block[(a*j+b)%7] for j in range(7)),z) for z,a in enumerate((1,2,4)) for b in range(7)])
   elif kind in ('R','C'):choices.append([(tuple(block[(j+b)%size] for j in range(size)),b if kind=='C' else None) for b in range(size)])
   else:choices.append([((block[0],),None)])
 G=[]
 for terms in it.product(*choices):
  labels=[z for _,z in terms if z is not None]
  if coeff is not None and sum(a*b for a,b in zip(coeff,labels))%3:continue
  G.append(tuple(x for perm,z in terms for x in perm))
 return blocks,G

def parse_arcs(line):return {tuple(map(int,t.split(','))) for t in line.split()}

def audit(dirname,ns):
 folder=ROOT/dirname; exp=expected(ns); actual={x.name for x in folder.glob('*.job')}
 need(actual==set(exp),f'{dirname}: job-name mismatch')
 configs=set(); totals={'jobs':0,'group_elements':0,'pair_orbitals':0,'modules':0,'nonmodules':0}
 for name,args in sorted(exp.items()):
  lines=(folder/name).read_text().splitlines(); n,K,gs=map(int,lines[0].split())
  blocks,G=group(*args);configs.add(args[:5]);need(n==args[0] and gs==len(G),name+': group size')
  need([list(map(int,t.split('.'))) for t in lines[1].split()]==blocks,name+': block list')
  need(len(lines)==K+3,name+': lines')
  fixed=parse_arcs(lines[2]); expectedfixed=set()
  for b in blocks:
   steps={1:( ),3:(1,),5:(1,2),7:(1,2,4)}[len(b)]
   expectedfixed|={(b[i],b[(i+d)%len(b)]) for i in range(len(b)) for d in steps}
  need(fixed==expectedfixed,name+': fixed arcs')
  for b in blocks:need({g[b[0]] for g in G}==set(b),name+': nontransitive projection')
  allpairs={frozenset((u,v)) for u in range(n) for v in range(u+1,n)}
  used={frozenset(pr) for pr in fixed};need(len(used)==len(fixed),name+': duplicated fixed inverse')
  rect=True
  blockof={v:b for b in blocks for v in b}
  for line in lines[3:]:
   orb=parse_arcs(line);need(bool(orb),name+': empty orbital');u,v=min(orb)
   trueorb={(g[u],g[v]) for g in G};need(orb==trueorb,name+': invalid orbital')
   unor={frozenset(pr) for pr in orb};need(len(unor)==len(orb) and used.isdisjoint(unor),name+': inverse/repeated orbital')
   need(blockof[u]!=blockof[v],name+': noncross orbital');used|=unor
   rect &=len(orb)==len(blockof[u])*len(blockof[v])
  need(used==allpairs,name+': missing pair')
  totals['jobs']+=1;totals['group_elements']+=len(G);totals['pair_orbitals']+=K;totals['modules' if rect else 'nonmodules']+=1
 totals['configs']=len(configs)
 totals['job_tree_sha256']=hashlib.sha256(b''.join(name.encode()+b'\0'+(folder/name).read_bytes() for name in sorted(exp))).hexdigest()
 return totals

out={'window':audit('jobs-normal',list(range(13,21))+[22]),'n12_alternate':audit('jobs-n12',[12])}
for a,b in [('jobs-normal','jobs-optimized'),('jobs-n12','jobs-n12-opt')]:
 for p in (ROOT/a).glob('*.job'):need(p.read_bytes()==(ROOT/b/p.name).read_bytes(),p.name+': -O mismatch')
print(json.dumps(out,sort_keys=True,indent=2))
