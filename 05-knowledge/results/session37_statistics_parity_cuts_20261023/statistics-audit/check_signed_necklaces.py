#!/usr/bin/env python3
"""Independent weighted-word rotation and PBW Hilbert-product audit."""
import json
from fractions import Fraction
from math import gcd
from pathlib import Path
P=Path(__file__).resolve().parent

def check(ok,msg):
 if not ok:raise AssertionError(msg)
def conv(a,b):
 c=[0]*9
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   if i+j<=8:c[i+j]+=x*y
 return c

def main():
 a=[r['types']for r in json.loads((P/'characters.json').read_text())['rows'][:9]]
 s=[0]*9;e=[0]*9;logder=[0]*9
 for n in range(1,9):
  s[n]=a[n]-sum(s[k]*a[n-k]for k in range(1,n))
  logder[n]=n*a[n]-sum(logder[k]*a[n-k]for k in range(1,n))
  num=logder[n]-sum(k*e[k]for k in range(1,n)if n%k==0)
  check(num%n==0,'Euler integrality');e[n]=num//n
 words=[[()]]+[[]for _ in range(8)];rows=[];prod=[1]+[0]*8
 for n in range(1,9):
  for k in range(1,n+1):
   for j in range(s[k]):
    for w in words[n-k]:words[n].append(((k,j),)+w)
  check(len(words[n])==a[n],('word count',n))
  necklaces={min(w[i:]+w[:i]for i in range(len(w)))for w in words[n]}
  d=primitive=0;burn=Fraction(0)
  for w in words[n]:
   wt=0;val=0
   for i in range(len(w)):
    if w[i:]+w[:i]==w:val+=(-1)**(wt*(n-wt))
    wt+=w[i][0]
   burn+=Fraction(val,len(w))
  for w in necklaces:
   length=len(w);period=next(k for k in range(1,length+1)if length%k==0 and w==w[:k]*(length//k))
   weight=sum(k for k,j in w[:period]);r=length//period
   if r==1:primitive+=1
   good=not(weight%2 and r%2==0)
   wt=0;bad=False
   for i in range(length):
    if w[i:]+w[:i]==w and wt*(n-wt)%2:bad=True
    wt+=w[i][0]
   check(good==(not bad),('root sign',n,w));d+=good
  c=sum(e[k]for k in range(1,n+1)if n%k==0)
  expected=c-(sum(e[k]for k in range(1,n+1)if n%k==0 and k%2)if n%2==0 else 0)
  formula=sum(Fraction(sum(gcd(i,k)==1 for i in range(1,k+1)),k)*Fraction(logder[n//k],n//k)*(-1)**((k-1)*(n//k))for k in range(1,n+1)if n%k==0)
  check(primitive==e[n] and len(necklaces)==c and d==expected==burn==formula,('necklace totals',n))
  beta=e[n]+(e[n//2]if n%4==2 else 0)if n%2==0 else 0
  phi=e[n]if n%2 else 0
  factor=[0]*9;factor[0]=1
  if phi:
   factor[n]=1
   for _ in range(phi):prod=conv(prod,factor)
  if beta:
   for j in range(n,9,n):factor[j]=1
   for _ in range(beta):prod=conv(prod,factor)
  rows.append({'n':n,'strong_types':s[n],'ordinary_primitive_necklaces':e[n],'ordinary_necklaces':c,'signed_necklaces':d,'even_Lie_dimension':beta,'odd_Lie_dimension':phi})
 check(prod==a,'PBW Hilbert series')
 receipt={'status':'PASS','rows':rows,'signed_rotation_Burnside':True,'signed_log_formula':True,'root_parity_predicate':True,'PBW_Hilbert_product':True,'word_count':sum(a[1:])}
 (P/'signed_necklaces.json').write_text(json.dumps(receipt,indent=2)+'\n');print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
