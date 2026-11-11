"""Discovery probe: exact affine-in-2^k intersections of two physical twin languages."""
from functools import cache
from collections import defaultdict
from pathlib import Path
import json
HERE=Path(__file__).resolve().parent
@cache
def g(n):
 if n<2:return 1
 q=n//2;return g(q)+(q+(q&1) if not(n&1) else 0)
def rev(n,l):return int(format(n,'0%db'%l)[::-1],2) if l else 0
def append_coeff(r,l):
 c=g(r)-g(0)
 a=(g((2<<l)+r)-g(2)-c)//2
 b=g((1<<l)+r)-g(1)-a-c
 return a,b,c

def families(maxprefix=512,maxsuffix=5,kind='old'):
 d,o,e,h=(13,32,6,18) if kind=='old' else (25,122,16,48)
 ubase=d+o
 for p in range(maxprefix+1):
  lp=p.bit_length();rp=rev(p,lp)
  for l in range(maxsuffix+1):
   for r in range(1<<l):
    a,b,c=append_coeff(r,l)
    for digit in (0,1):
     if p==0 and digit==0:continue
     parity=(r&1) if l else digit
     if kind=='new' and parity!=1:continue
     if digit==0:
      addrA=p*(1<<l);addrB=r
      gA=(1+a)*p;gB=g(p)-p+(p&1)+c
      UA=(ubase*(1<<l)-d*rev(r,l))*(1<<lp);UB=-d*rp
     else:
      addrA=(p+1)*(1<<l);addrB=r-(1<<l)
      gA=a*(p+1);gB=g(p)-a+b+c
      UA=(ubase*(1<<l)-d*(rev(r,l)+1))*(1<<lp);UB=d*((1<<lp)-rp)
     F_A=d*gA+o*addrA;F_B=d*gB+o*addrB+e*parity+h
     rec={'p':p,'loop':digit,'suffix':r,'suffix_length':l,'address_coeff':[addrA,addrB],'g_coeff':[gA,gB],'H_coeff':[F_A,F_B],'U_coeff':[UA,UB]}
     # exact consequence controls; loop family starts k>=1
     for k in (1,2,3):
      B=addrA*(1<<k)+addrB
      if g(B)!=gA*(1<<k)+gB:raise RuntimeError(('g',rec,k))
      if ubase*(1<<B.bit_length())-d*rev(B,B.bit_length())!=UA*(1<<k)+UB:raise RuntimeError(('U',rec,k))
     yield rec

def main():
 old=defaultdict(list)
 for f in families():old[tuple(f['H_coeff'])].append(f)
 hits=[]
 for new in families(kind='new'):
  for first in old.get(tuple(new['H_coeff']),[]):
   da=first['U_coeff'][0]-new['U_coeff'][0];db=first['U_coeff'][1]-new['U_coeff'][1]
   if da>=0 and 2*da+db>0:
    hits.append({'old':first,'new':new,'defect_drop_coeff':[da,db]})
 # deduplicate address parametrizations
 uniq={tuple(h['old']['address_coeff']+h['new']['address_coeff']):h for h in hits}
 ans=sorted(uniq.values(),key=lambda h:(h['old']['H_coeff'][0]*2+h['old']['H_coeff'][1],sum(h['old']['address_coeff']+h['new']['address_coeff'])))
 report={'prefix_max':512,'suffix_length_max':5,'universe':'p digit^k suffix, 0<=p<=512 (empty prefix only for one loops), fixed suffix length<=5, k>=1; both loop digits; odd new addresses; exact affine coefficient equality and U drop for everyk>=1','positive_intersections':len(ans),'first':ans[:40]}
 (HERE/'intersections.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report))
if __name__=='__main__':main()
