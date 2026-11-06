from itertools import product
from probe import T,state,step

def G(s):return s[0]+2*s[2]
def F(s):return s[0]-s[1]+s[2]
def source_words(n,readout=G):
 result=[];stack=[(T,'')]
 while stack:
  s,w=stack.pop();v=readout(s)
  if v==n:result.append(w)
  if v<n:
   for c in 'OZJ':
    q=step(s,c)
    if readout(q)<=n:stack.append((q,w+c))
 return result

def runs(w):
 out=[]
 for c in w:
  if out and out[-1][0]==c:out[-1]=(c,out[-1][1]+1)
  else:out.append((c,1))
 return out

if __name__=='__main__':
 W=source_words(227)
 print('words',len(W),W[:30],flush=True)
 matches=[];seen=set();tests=0
 for w in W:
  rs=runs(w);vary=[i for i,r in enumerate(rs)if r[0]!='J']
  for increments in product(range(7),repeat=len(vary)):
   if not 4<=sum(increments)<=8:continue
   off=dict(zip(vary,increments))
   make=lambda r:''.join(c*(k+off.get(i,0)*(r-1))for i,(c,k)in enumerate(rs))
   tests+=1
   if G(state(make(2)))!=3971:continue
   if G(state(make(3)))!=65027:continue
   if G(state(make(4)))!=1046531:continue
   matches.append((w,increments))
 print('tests',tests,'matches',matches)
