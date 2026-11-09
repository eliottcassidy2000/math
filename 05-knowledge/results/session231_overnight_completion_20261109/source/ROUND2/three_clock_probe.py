import json,math
from pathlib import Path
Y=[tuple(map(int,s.split())) for s in Path('../bank.txt').read_text().splitlines() if int(s.split()[0])>1]
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def mul(x,c):return tuple(c*a for a in x)
def z(x):
 if x[-1]:raise RuntimeError('degree')
 return (0,)+x[:-1]
def fixed(st,w):
 v,g,p=st
 for bit in w:
  if bit=='0':g=add(g,add(v,(p,0,0,0)))
  v=add(mul(v,2),(int(bit),0,0,0));p=int(bit)
 return v,g,p
def run(st,b):
 v,g,p=st
 if b=='1':return add(z(v),(-1,1,0,0)),g,1
 return z(v),add(g,add(add(z(v),mul(v,-1)),(p,0,0,0))),0
rat={}
for a,b,c in Y:
 h=math.gcd(a,c);rat.setdefault((c//h,a//h),[]).append((a,b,c))
words=['','0','1','00','01','10','11'];found=[];count=0
for pref in ['','1','10','11']:
 for b1 in '01':
  st=run(fixed(((0,0,0,0),(1,0,0,0),0),pref),b1)
  for mid1 in words:
   for b2 in '01':
    st2=run(fixed(st,mid1),b2)
    for mid2 in words:
     for b3 in '01':
      st3=run(fixed(st2,mid2),b3)
      for suff in words:
       v,g,p=fixed(st3,suff);count+=1
       ys=Y
       for i in [1,2]:
        if v[i]==0:
         if g[i]:ys=[];break
        else:
         cc,aa=-g[i],v[i]
         if aa<0:aa,cc=-aa,-cc
         if cc<=0:ys=[];break
         h=math.gcd(aa,cc);ys=rat.get((cc//h,aa//h),[]);break
       for a,b,c in ys:
        if any(a*g[i]+c*v[i] for i in [1,2]):continue
        A=a*g[3]+c*v[3]
        if A<=0 or A&(A-1):continue
        # Require three distinct variable blocks not trivially merged.
        genuine=(b1!=b2 or any(t!=b1 for t in mid1)) and (b2!=b3 or any(t!=b2 for t in mid2))
        if genuine:found.append({'row':[a,b,c],'template':[pref,b1,mid1,b2,mid2,b3,suff],'A':A,'D':a*g[0]+c*(v[0]+1)-b*(1-p),'v_poly':v,'g_poly':g})
print(json.dumps({'templates':count,'genuinely_separated_actual_bank_hits':found},indent=2))
