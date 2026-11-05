#!/usr/bin/env python3
from itertools import product
from fractions import Fraction
import json

def need(p,s):
 if not p:raise RuntimeError(s)
def trans(n):return tuple(sum(1<<j for j in range(i+1,n))for i in range(n))
def birth(r,p,op):
 side,mask,drop=op;n=len(r);a=list(r)+[0]
 for v in range(n):
  wins=(v==p and mask==1)
  if side:wins=not wins
  if wins:a[n]|=1<<v
  else:a[v]|=1<<n
 return tuple(a),n if drop==0 else p
def order(r):
 n=len(r);w=sorted(range(n),key=lambda v:-r[v].bit_count())
 if all(r[w[i]]>>w[j]&1 for i in range(n)for j in range(i+1,n)):return w
 return None
def category(r,p):
 w=order(r)
 if w is None:return 'X'
 if len(r)==1:return 'U'
 return 'B'if p in(w[0],w[-1])else'I'
def hp(r):
 n=len(r);d=[[0]*n for _ in range(1<<n)]
 for v in range(n):d[1<<v][v]=1
 for S in range(1,1<<n):
  for v,c in enumerate(d[S]):
   if not c:continue
   for w in range(n):
    if not S>>w&1 and r[v]>>w&1:d[S|1<<w][w]+=c
 return sum(d[-1])
A=list(product(range(2),range(2),range(2)))
def main():
 transitions=exits=0
 for n in range(1,10):
  for i in range(n):
   r=trans(n);counts={s:0 for s in'BIX'}
   for op in A:
    out,p=birth(r,i,op);c=category(out,p);counts[c]+=1;transitions+=1
    if c=='X':
     need(op[1]==1,'exit selected port');d=i if op[0]else n-1-i
     need(d>=1 and hp(out)==2**d+1,'exact first-exit physical H');exits+=1
   expected={'B':8,'I':0,'X':0}if n==1 else({'B':4,'I':2,'X':2}if i in(0,n-1)else{'B':2,'I':2,'X':4})
   need(counts==expected,'physical transition multiplicities')
 fib=[0,1]
 for j in range(405):fib.append(fib[-1]+fib[-2])
 b=i=1
 for m in range(201):
  need(b==2**m*fib[2*m+2]and i==2**m*fib[2*m+1],'Fibonacci closed forms')
  b,i=4*b+2*i,2*b+2*i
 prefix=[]
 for name,r,p in [('U',trans(1),0),('B',trans(2),0),('I',trans(3),1)]:
  frontier=[(r,p)];counts=[1]
  for depth in range(1,6):
   nxt=[]
   for r,p in frontier:
    for op in A:
     out,q=birth(r,p,op)
     if order(out)is not None:nxt.append((out,q))
   frontier=nxt;counts.append(len(frontier))
   want=2**depth*fib[2*depth+(2 if name=='B'else 1)]if name!='U'else 2**(depth+2)*fib[2*depth]
   need(len(frontier)==want,'complete word-prefix count')
  prefix.append({'initial_category':name,'counts_through_depth5':counts})
 z=Fraction(1,8);den=1-6*z+4*z*z
 need(1/den==Fraction(16,5)and(1-2*z)/den==Fraction(12,5)and 1+8*z/den==Fraction(21,5),'exact expectations')
 witnesses=0
 for n in range(1,5):
  for root in range(n):
   for d in range(1,6):
    r,p=trans(n),root;r,p=birth(r,p,(0,0,0));r,p=birth(r,p,(0,1,0))
    for _ in range(d-1):r,p=birth(r,p,(0,0,1))
    need(order(r)is not None,'pre-exit still transitive');r,p=birth(r,p,(0,1,1))
    need(hp(r)==2**d+1,'all first-exit values constructed');witnesses+=1
 print(json.dumps({'status':'PASS','physical_transitions':transitions,'literal_first_exit_counts':exits,'complete_word_prefixes':prefix,'Fibonacci_indices_checked':201,'first_exit_converse_witnesses':witnesses,'exact_mean_first_exit_times':{'B':'16/5','I':'12/5','U':'21/5'}},sort_keys=True,indent=2))
if __name__=='__main__':main()
