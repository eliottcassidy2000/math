"""Exact fixed-prefix/monochrome-tail templates, with algebraic coefficient filter."""
import subprocess,json
from pathlib import Path
from head_lift_probe import ends,converse
ns={'__name__':'inherited'};exec(subprocess.check_output(['git','show','f92e070445:05-knowledge/results/session251_creative_completion_20261111/ternary/compiler.py']).decode(),ns)
g=ns['g'];count=ns['count']
records=json.loads(subprocess.check_output(['git','show','f92e070445:05-knowledge/results/actual_hamilton_multiplicities_20261021/session24-multiplicities/extended_census.json']))['records']
bank=[(0,),(2,0)]
for r in records:
 C=tuple(r['adjacency_outmasks']);bank.append(C)
 if r['ordinary_multiplicity']==2:bank.append(tuple(converse(C)))

def signatures(C):
 n=len(C);M=(1<<n)-1;h=[sum(v) for v in ends(C)];h[0]=1;P=[[0]*n for _ in C];Q=[[0]*n for _ in C]
 def walk(mask,start,end):
  w=h[M^mask];P[start][end]+=w;Q[start][end]+=w*(-1 if mask.bit_count()%2 else 1)
  avail=C[end]&~mask&M
  while avail:
   b=avail&-avail;avail-=b;walk(mask|b,start,b.bit_length()-1)
 for i in range(n):walk(1<<i,i,i)
 kp=sum(map(sum,P));km=sum(map(sum,Q))
 for O in range(M+1):
  I=[i for i in range(n) if not O>>i&1];J=[i for i in range(n) if O>>i&1]
  ip=sum(P[i][j] for i in range(n) for j in I);im=sum(Q[i][j] for i in range(n) for j in I)
  op=sum(P[i][j] for i in J for j in range(n));om=sum(Q[i][j] for i in J for j in range(n))
  jp=sum(P[i][j] for i in J for j in I);jm=sum(Q[i][j] for i in J for j in I)
  yield O,(h[M],kp,km,ip,im,op,om,jp,jm)

hits=[];seen=set();total=0
for C in bank:
 for O,sig in signatures(C):
  if sig in seen:continue
  seen.add(sig);d,kp,km,ip,im,op,om,jp,jm=sig;Q2=d+kp
  for side in ('right','left'):
   linear=op if side=='right' else ip
   for odd in (1,5,7,35):
    u=35//odd
    for t in range(13):
     L=odd*2**t
     if L<linear or (L-linear)%Q2:continue
     fixed=(L-linear)//Q2
     if side=='right':
      const=d*(g(fixed)+g(u-1)-1)-L+ip*fixed-im*(fixed%2)+(d+km)*(fixed%2)-om+jp+jm
     else:
      const=d*(g(u-1)+g(fixed)-1)-L+op*fixed-om*(fixed%2)+(d+km)*(fixed%2)-im+jp+jm
     if const not in (9,17,33):continue
     s={9:3,17:4,33:5}[const];start=max(1,t-s+1)
     hit={'core':C,'cut':O,'signature':sig,'side':side,'fixed':fixed,'u':u,'t':t,'s':s,'start_r':start}
     hits.append(hit);print(json.dumps(hit),flush=True)
     for r in range(start,start+3):
      moving=u*2**(r+s-t)-1;A,B=(fixed,moving) if side=='right' else(moving,fixed)
      assert count(sig,A,B)==1+2**s+35*2**(r+s)
Path(__file__).with_name('one_tail_receipt.json').write_text(json.dumps({'unique_signatures':len(seen),'hits':hits},indent=2)+'\n')
print('DONE',len(seen),'signatures',len(hits),'hits',flush=True)
