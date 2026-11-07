from inherited import *
import json

def odd(n):
 k=0
 while n%2==0:n//=2;k+=1
 return n,k
banks=[]
for pair,js,root,a,b in PAIRS:
 for side,rows,j in zip('AB',(a,b),js):
  init=profile(remark(rows,root));seen={init:''};level={init:''};bank={}
  for dep in range(11):
   nxt={}
   for p,w in level.items():
    x=base(p,j);h,s,e,L,R,d,C,b,T=p;coef=R+C;offset=x-2*L
    o,k=odd(coef)
    key=(offset,o)
    if key not in bank or ('J' in w and 'J' not in bank[key][2]):bank[key]=(pair,side,w,p,j,k)
    if dep<10:
     for c in 'OZJ':
      q=move(p,c)
      if q not in seen: seen[q]=w+c;nxt[q]=w+c
   level=nxt
  banks.append(bank)
  print(pair,side,len(bank),flush=True)
results=[]
for i in range(0,8,2):
 for k in range(1,8,2):
  for key in banks[i].keys()&banks[k].keys():
   a,b=banks[i][key],banks[k][key]
   if 'J' in a[2]+b[2]:results.append((key,a,b))
results.sort(key=lambda z:(max(len(z[1][2]),len(z[2][2])),z[0]))
print('FAMILIES',len(results));print(results[:10]);json.dump(results,open('jrays.json','w'),indent=2)
