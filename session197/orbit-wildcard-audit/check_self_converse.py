#!/usr/bin/env python3
import hashlib
import itertools
import json
import math
from pathlib import Path
from check_graph import check,graph,h_walk_inclusion_exclusion,automorphisms,verify_certificate,trial_prime

SOURCE=Path('/workspace/scratch/07d5532f2fce/session197-orbit-wildcard')
HASHES={'SELF_CONVERSE.md':'79c2fa74478e485df184c7b4b9831b08f14c67008d42bcea5c31afec0ff553a9',
        'verify_self_converse.py':'6af2f980b72d638e5cb167035fdf9842e5d1643f8db93d88a53508da1f099698',
        'self-converse-normal.json':'c35d8ba1b4da4157d5521ca545be0010cd68b1d73e85803762c5440347beedd7',
        'self-converse-witness.json':'6c52d4b98277eb0e83927ba33ac0bab2b5e508f17e38894d07704c080ab3df17'}

def quotient_substitution(Q,U):
 # Root0 is replaced by U, then quotient vertices1,2,... are appended.
 m=len(U);n=m+len(Q)-1;rows=U[:]+[0]*(len(Q)-1)
 for i in range(len(Q)):
  for j in range(len(Q)):
   if not Q[i]>>j&1:continue
   left=range(m) if i==0 else [m+i-1]
   right=range(m) if j==0 else [m+j-1]
   for x in left:
    for y in right:rows[x]|=1<<y
 return rows

def run_counts(Q):
 n=len(Q);a=[0]*(n+1)
 for perm in itertools.permutations(range(1,n)):
  for gaps in range(1,1<<n):
   word=[]
   for k in range(n):
    if gaps>>k&1:word.append(0)
    if k<n-1:word.append(perm[k])
   if all(Q[v]>>w&1 for v,w in zip(word,word[1:])):a[gaps.bit_count()]+=1
 return a

def ordered_covers(U):
 n=len(U);answer=[0]*4
 for perm in itertools.permutations(range(n)):
  failures=sum(not(U[v]>>w&1) for v,w in zip(perm,perm[1:]))
  for k in range(1,4):
   if failures<=k-1:answer[k]+=math.comb(n-1-failures,k-1-failures)
 return answer

def fixed_words(rows,tau):
 n=len(rows);check(sorted(tau)==list(range(n)),'anti permutation')
 check(all(tau[tau[v]]==v for v in range(n)),'anti involution')
 check(all(((rows[u]>>v)&1)==((rows[tau[v]]>>tau[u])&1)
           for u in range(n) for v in range(n) if u!=v),'anti arcs')
 fixed=[v for v in range(n) if tau[v]==v]
 check(len(fixed)==n%2,'correct central fixed point count')
 pairs=[];done=set(fixed)
 for v in range(n):
  if v not in done:pairs.append([v,tau[v]]);done.update([v,tau[v]])
 k=len(pairs);pair_of={v:i for i,P in enumerate(pairs) for v in P}
 central=(lambda v: bool(rows[v]>>fixed[0]&1)) if fixed else (lambda v: bool(rows[v]>>tau[v]&1))
 dp=[[0]*n for _ in range(1<<k)]
 for v in pair_of:dp[1<<pair_of[v]][v]=1
 for M in range(1,1<<k):
  for v in pair_of:
   if not dp[M][v]:continue
   for w in pair_of:
    b=1<<pair_of[w]
    if not M&b and rows[v]>>w&1:dp[M|b][w]+=dp[M][v]
 F=sum(dp[-1][v] for v in pair_of if central(v))
 IE=0
 for M in range(1,1<<k):
  vertices=[v for v in pair_of if M>>pair_of[v]&1]
  vector=[int(v in vertices) for v in range(n)]
  edges=[(v,w) for v in vertices for w in vertices if rows[v]>>w&1]
  for _ in range(k-1):
   nextv=[0]*n
   for v,w in edges:nextv[w]+=vector[v]
   vector=nextv
  IE+=(-1 if (k-M.bit_count())&1 else 1)*sum(vector[v] for v in vertices if central(v))
 check(F==IE,'pair-subset DP and pair-color walk count')
 return F,fixed

def check_word(rows,tau,word):
 check(sorted(word)==list(range(len(rows))),'fixed word labels')
 check(all(rows[v]>>w&1 for v,w in zip(word,word[1:])),'fixed word arcs')
 check(word==[tau[v] for v in word[::-1]],'fixed word action')

def main():
 for name,h in HASHES.items():check(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==h,('frozen extension hash',name))
 data=json.loads((SOURCE/'self-converse-witness.json').read_text())
 Q=[6,24,26,1,9];U=[sum(1<<((i+d)%7) for d in [1,2,3]) for i in range(7)]
 rows=quotient_substitution(Q,U)
 check(rows==data['eleven']['rows'],'independent11 construction')
 H,_=h_walk_inclusion_exclusion(rows);check(H==56483==7*8069,'independent11 H')
 check(trial_prime(8069),'11 quotient primality')
 certs={};verify_certificate(data['eleven']['prime_certificate'],certs)
 run=run_counts(Q);covers=ordered_covers(U)
 check(run==[0,11,17,4,0,0],'complete quotient run counts')
 check(covers[1:]==[175,1638,6678],'literal ordered carousel covers')
 check(sum(run[k]*covers[k] for k in range(1,4))==H,'physical substitution profile formula')
 auto,nodes=automorphisms(rows);check(len(auto)==7,'full11 automorphism order')
 rotations=[[((i+k)%7 if i<7 else i) for i in range(11)] for k in range(7)]
 check(sorted(auto)==sorted(rotations),'full11 rotation group')
 reach=[rows[v]|1<<v for v in range(11)]
 for k in range(11):
  for v in range(11):
   if reach[v]>>k&1:reach[v]|=reach[k]
 check(all(x==(1<<11)-1 for x in reach),'strong11 graph')
 tq=[0,4,3,2,1]
 check(all(((Q[i]>>j)&1)==((Q[tq[j]]>>tq[i])&1) for i in range(5) for j in range(5) if i!=j),'quotient anti-involution')
 F11,fix11=fixed_words(rows,data['eleven']['anti_involution'])
 check(F11==111==data['eleven']['fixed_count'] and fix11==[0],'eleven fixed words')
 check_word(rows,data['eleven']['anti_involution'],data['eleven']['fixed_word'])
 R=graph();F14,fix14=fixed_words(R,data['fourteen']['anti_involution'])
 check(F14==6111==data['fourteen']['fixed_count'] and fix14==[],'fourteen fixed words')
 check_word(R,data['fourteen']['anti_involution'],data['fourteen']['fixed_word'])
 tau=data['fourteen']['anti_involution'];rho=[7*(i//7)+(i%7+1)%7 for i in range(14)]
 rinv=[7*(i//7)+(i%7-1)%7 for i in range(14)]
 check([tau[rho[tau[v]]] for v in range(14)]==rinv,'dihedral relation')
 check((15144269+7*F14)//14==1084789==data['fourteen']['dihedral_word_orbits'],'Burnside count')
 print(json.dumps({'status':'PASS','frozen_extension_hashes':HASHES,'eleven_H':H,
                   'prime_quotient':8069,'prime_trial_bound':math.isqrt(8069),'Lucas_nodes':certs,
                   'quotient_run_counts':run[1:],'carousel_ordered_covers':covers[1:],
                   'eleven_automorphism_order':len(auto),'eleven_automorphism_search_nodes':nodes,
                   'eleven_fixed_Hamilton_words':F11,'fourteen_fixed_Hamilton_words':F14,
                   'dihedral_word_orbits_fourteen':1084789,
                   'scope':'Explicit self-converse witnesses and fixed-word counts independently checked; minimum-order proof reviewed separately.'},
                  indent=2,sort_keys=True))

if __name__=='__main__':main()
