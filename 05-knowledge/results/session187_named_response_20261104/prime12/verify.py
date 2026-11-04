"""Independent integer subset DP, exact ordered-cover convolution and literal queries."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def need(x,msg):
 if not x:raise RuntimeError(msg)
def op(a):return [((1<<len(a))-1)^x^(1<<u) for u,x in enumerate(a)]
def endpoints(a):
 n=len(a);inc=op(a);d=[[0]*n for _ in range(1<<n)]
 for m in range(1,1<<n):
  z=m
  while z:
   bit=z&-z;u=bit.bit_length()-1;r=m^bit
   if not r:d[m][u]=1
   else:
    x=r&inc[u]
    while x:
     b=x&-x;v=b.bit_length()-1;d[m][u]+=d[r][v];x^=b
   z^=bit
 return d
def hp(a):return sum(endpoints(a)[-1])
def ex(a,x,c):
 n=len(a);b=list(a)+[1<<(n+1),0]
 for u in range(n):
  if x>>u&1:b[u]|=1<<n
  else:b[n]|=1<<u
  if c>>u&1:b[u]|=1<<(n+1)
  else:b[n+1]|=1<<u
 return b
def mods(a):
 f=(1<<len(a))-1;out=[]
 for m in range(1,f):
  if m.bit_count()<2:continue
  if all(a[u]&m in (0,m) for u in range(len(a)) if not m>>u&1):out.append(m)
 return out
def main():
 d=json.loads((ROOT/'cubic.json').read_text());a=d['rows'];n=len(a);f=(1<<n)-1
 E=endpoints(a);S=endpoints(op(a));need(sum(E[f])==84599==d['H'],'Hamilton count')
 need(all((x-y)%2==0 for x,y in zip(E[f],S[f])),'endpoint parity');need(not mods(a),'prime');need(sorted(x.bit_count() for x in a)!=sorted(x.bit_count() for x in op(a)),'nonselfconverse score certificate')
 C=[[0]*n for _ in range(n)]
 for m in range(1,f):
  for u in range(n):
   if not m>>u&1:continue
   for v in range(n):
    if m>>v&1:continue
    C[u][v]+=E[m][u]*S[f^m][v]
 bad=[[u,v,C[u][v]-C[v][u]] for u in range(n) for v in range(u+1,n) if (C[u][v]-C[v][u])%4]
 need(bool(bad),'cover-skew witness')
 coeff={(x,y):v for x,y,v in d['coefficients']};need(len(coeff)==1885,'complete cubic coordinate count')
 need(all(v==0 for (x,y),v in coeff.items() if (x==0 and y.bit_count()==2) or (y==0 and x.bit_count()==2)),'pure quadratics')
 cache={}
 def delta(x,y):
  if (x,y) not in cache:cache[x,y]=hp(ex(a,x,y))-hp(ex(op(a),x,y))
  return cache[x,y]
 checked=[]
 for x,y in [(1,3),(3,4),(3,3)]:
  c=0;xx=x
  while True:
   yy=y
   while True:
    c+=(-1)**((x^xx).bit_count()+(y^yy).bit_count())*delta(xx,yy)
    if yy==0:break
    yy=(yy-1)&y
   if xx==0:break
   xx=(xx-1)&x
  need(c==coeff.get((x,y),0),'independent literal finite difference');checked.append([x,y,c])
 need(coeff[1,3]==-56686,'named modular separator')
 for x,y in [(341,2730),(f,f)]:
  z=sum(v for (i,j),v in coeff.items() if i&x==i and j&y==j)
  need(z==delta(x,y),'independent full assignment')
 return {'status':'pass','rows':a,'H':sum(E[f]),'E':E[f],'S':S[f],'prime':True,'full_C':C,'cover_skew_mod4_failures':bad,'literal_coefficients':checked,'full_named_cubic_coordinate_count':len(coeff),'nonzero_exact_coefficients':sum(v!=0 for v in coeff.values()),'nonzero_modular_cubics':sum(v%4!=0 for (x,y),v in coeff.items() if x.bit_count()+y.bit_count()==3),'literal_delta_queries':[[x,y,v] for (x,y),v in sorted(cache.items())]}
if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
