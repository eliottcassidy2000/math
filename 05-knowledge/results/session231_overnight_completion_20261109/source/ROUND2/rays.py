import json,math
from pathlib import Path
Y=[tuple(map(int,s.split())) for s in Path('../bank.txt').read_text().splitlines()]
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def scale(x,n):return tuple(a*n for a in x)
def shift(x):
 if x[-1]:raise RuntimeError('degree')
 return (0,)+x[:-1]
def fixed(state,w):
 v,g,p=state
 for b in w:
  if b=='0':g=add(g,add(v,(p,0,0)))
  v=add(scale(v,2),(int(b),0,0));p=int(b)
 return v,g,p
def run(st,b):
 v,g,p=st
 if b=='0':return shift(v),add(g,add(add(shift(v),scale(v,-1)),(p,0,0))),0
 return add(shift(v),(-1,1,0)),g,1
rat={}
for y in Y:
 a,b,c=y;d=math.gcd(a,c);rat.setdefault((c//d,a//d),[]).append(y)
words=['']+[f'{r:0{k}b}' for k in range(1,5) for r in range(2**k)]
found={};templates=0
for pref in ['']+[bin(x)[2:] for x in range(1,16)]:
 for r1 in '01':
  st1=run(fixed(((0,0,0),(1,0,0),0),pref),r1)
  for mid in words:
   for r2 in '01':
    st2=run(fixed(st1,mid),r2)
    for suff in words:
     v,g,p=fixed(st2,suff);templates+=1
     if not v[2]:continue
     if v[1]==0:
      if g[1]:continue
      ys=Y
     else:
      c0,a0=-g[1],v[1]
      if a0<0:a0,c0=-a0,-c0
      if c0<=0:continue
      d=math.gcd(a0,c0);ys=rat.get((c0//d,a0//d),[])
     for a,b,c in ys:
      A=a*g[2]+c*v[2];B=a*g[1]+c*v[1];D=a*g[0]+c*(v[0]+1)-b*(1-p)
      if B or A<=0 or A&(A-1) or abs(D)>255 or D%2==0:continue
      j=A.bit_length()-1;key=(j%2,D)
      row={'bank':[a,b,c],'template':[pref,r1,mid,r2,suff],'coefficient_power':j,'constant':D,'v_poly':v,'g_poly':g,'parity':p,'minimum_exponent':j+2}
      score=(row['minimum_exponent'],a+c,sum(map(len,[pref,mid,suff])))
      if key not in found or score<found[key][0]:found[key]=(score,row)
rows=[found[k][1] for k in sorted(found)]
Path('RAYS.json').write_text(json.dumps(rows,indent=2)+'\n')
print('templates',templates,'rays',len(rows))
for parity in [0,1]:
 missing=[d for d in range(-255,256,2) if (parity,d) not in found]
 print('parity',parity,'covered',256-len(missing),'missing',missing)
print('paired offsets',[d for d in range(-255,256,2) if all((p,d) in found for p in [0,1])])
