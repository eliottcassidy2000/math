"""Independent subset-DP/permutation controls for corridor/PROOF.md."""
import itertools,json
from collections import Counter
from pathlib import Path
from charge_probe import charge,signature,count
from c5probe import form,run

def H(t):
 n=len(t);dp=[[0]*n for _ in range(1<<n)]
 for i in range(n):dp[1<<i][i]=1
 for mask in range(1,1<<n):
  for i in range(n):
   if not(mask>>i&1):continue
   z=dp[mask][i];rem=t[i]&~mask
   while rem:b=rem&-rem;rem^=b;dp[mask|b][b.bit_length()-1]+=z
 return sum(dp[-1]) if n else 1

def hp(t):return [p for p in itertools.permutations(range(len(t))) if all(t[u]>>v&1 for u,v in zip(p,p[1:]))]
def rows(l,r,a):
 n=l+5+r;t=[0]*(n+1);cy=[6,12,24,17,3]
 for i in range(n):
  for j in range(i+1,n):
   f=((cy[i-l]>>(j-l))&1) if l<=i<j<l+5 else 1
   if f:t[i]|=1<<j
   else:t[j]|=1<<i
  if a>>i&1:t[n]|=1<<i
  else:t[i]|=1<<n
 return t

def ttapex(n,a):
 t=[((1<<n)-1)^((1<<(i+1))-1) for i in range(n)]+[a]
 for i in range(n):
  if not(a>>i&1):t[i]|=1<<n
 return t

def fan(t,m,n):
 t=t[:]
 for _ in range(m):
  k=len(t);t.append(((1<<k)-1)^1);t[0]|=1<<k
 for _ in range(n):
  k=len(t);t.append(1)
  for i in range(1,k):t[i]|=1<<k
 return t

def main():
 checks=Counter();tab={};mins=[]
 for k in range(6):
  vals=[]
  for a in range(32):
   if a.bit_count()!=k:continue
   h=len(hp(rows(0,0,a)));assert h==H(rows(0,0,a));vals.append(h);checks['core_permutation_inputs']+=1
  tab[k]=sorted(set(vals))
 assert tab=={0:[15],1:[33],2:[41,43],3:[41,43],4:[33],5:[15]}
 for l in range(3):
  for r in range(3-l):
   hs=form(l,r)
   for a,h in enumerate(hs):assert H(rows(l,r,a))==h;checks['quadratic_form_direct_DP']+=1
   got=run(l,r);assert got['count_range']==[41,43] and got['component_states']==20
   checks['island_components']+=1
 # Everycore cut, both tail polarities, all distances through3, directly.
 for a in range(32):
  k=a.bit_count();base=H(rows(0,0,a))
  for d in range(4):
   left=(a<<(d+1))|1
   right=a|(((1<<d)-1)<<5)
   assert H(rows(d+1,0,left))==base+(120-24*k if d==0 else (1<<d)*(120-21*k))
   assert H(rows(0,d+1,right))==base+(24*k if d==0 else (1<<d)*(21*k+15))
   checks['anomaly_formula_direct_DP']+=2
 # Explicit boundary minima, recording actual cuts.
 for l,r in [(0,0),(1,0),(0,1)]:
  n=l+5+r;hs=form(l,r);island={(a<<l)|(((1<<r)-1)<<(l+5)) for a in range(32) if a.bit_count() in (2,3)}
  all_exit=[];tail_exit=[]
  moves=[1<<i for i in range(n)]+[(1<<i)|(1<<j) for i in range(n) for j in range(i+1,n)]
  for a in sorted(island):
   for d in moves:
    b=a^d
    if b in island:continue
    item=(abs(hs[b]-hs[a]),a,b,hs[a],hs[b])
    all_exit.append(item)
    if d&~(31<<l):tail_exit.append(item)
  assert min(all_exit)[0]==8
  if tail_exit:assert min(tail_exit)[0]==14
  mins.append({'left':l,'right':r,'all_exit_min':min(all_exit),'tail_exit_min':min(tail_exit) if tail_exit else None})
 t41=[56,13,25,48,34,6];t43=[52,13,24,49,34,6];t45=[52,13,24,33,42,6]
 assert list(map(lambda t:len(hp(t)),[t41,t43,t45]))==[41,43,45]
 frames=[]
 for a,b in [(t41,t43),(t43,t45)]:
  common=next(p for p in hp(a) if all(b[u]>>v&1 for u,v in zip(p,p[1:])))
  frames.append(common)
 for m in range(1,4):
  assert H(fan(t45,m,0))-H(fan(t43,m,0))==2**(m+1)-2
  checks['source_escape_difference']+=1
 for m,n in [(1,1),(1,2),(2,1),(2,2)]:
  assert H(fan(t45,m,n))-H(fan(t43,m,n))==-16*2**(m+n)+18*2**m+12*2**n-16
  checks['bifan_escape_difference']+=1
 for n in range(1,8):
  groups={}
  for a in range(1<<n):
   assert H(ttapex(n,a))==count(n,a);checks['one_apex_formula_direct_DP']+=1
   c=charge(n,a);sig=signature(n,a)
   assert c not in groups or groups[c]==sig
   groups[c]=sig
  assert len(groups)==3**(n//2)
 receipt={'checks':dict(checks),'core_cut_table':tab,'exact_exit_minima':mins,'common_frames':frames,'status':'PASS; finite controls separate from all-order proofs'}
 print(json.dumps(receipt,indent=2))
if __name__=='__main__':main()
