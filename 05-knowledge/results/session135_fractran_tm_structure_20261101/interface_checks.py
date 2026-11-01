"""Small exact cross-interface controls, not new all-order existence proofs."""
from pathlib import Path
from fractions import Fraction
from itertools import product,permutations
from math import isqrt,gcd
import json
ROOT=Path(__file__).resolve().parent

def need(x,m):
 if not x:raise RuntimeError(m)
def literal(rows,N):
 steps=0
 while True:
  for f in rows:
   if N%f.denominator==0:N=N//f.denominator*f.numerator;steps+=1;break
  else:return N,steps
  need(steps<10000,'literal guard')
def good(t):return all(min((j*t)%1,1-(j*t)%1)>=Fraction(1,8) for j in range(1,7))
def main():
 fractions=[Fraction(s.strip()) for s in (ROOT/'square_machines/square.fractions.txt').read_text().strip().split(',')]
 data=json.loads((ROOT/'square_machines/witnesses.json').read_text())
 cyc=next(x['cycle'] for x in data['cycles'] if x['n']==32);n=len(cyc);total=0;R=0
 need(sorted(cyc)==list(range(1,n+1)),'not a Hamiltonian witness')
 for u,v in zip(cyc,cyc[1:]+cyc[:1]):
  s=u+v;k=isqrt(s);need(k*k==s,'bad square edge')
  final,steps=literal(fractions,7*2**s)
  need(final==31*5**(2*k-1) and steps==5*s+2,'literal square cost')
  total+=steps;R+=k
 need(total==5*n*(n+1)+2*n and R==180,'cycle cost identities')
 tm=5*n*n+14*n+11+3*cyc[0]+6*R;need(tm==6662,'inherited literal TM cost')
 eighths=[]
 for b in (1,3,5,7):
  t=Fraction(b,8);owners=[]
  for j in range(1,7):
   residue=(j*t)%1
   if residue==Fraction(1,8):owners.append([j,'lower'])
   if residue==Fraction(7,8):owners.append([j,'upper'])
  left=good(t-Fraction(1,10000));right=good(t+Fraction(1,10000))
  need(good(t) and (left,right)==({1:(False,True),3:(False,False),5:(False,False),7:(True,False)}[b]),'eighth stratum')
  eighths.append({'t':str(t),'tight_owners':owners,'left_good':left,'right_good':right})
 circles=[];word_checks=0;path_checks=0
 for q in (3,5,7):
  for bits in product((0,1),repeat=(q-1)//2):
   S={i if bit==0 else q-i for i,bit in enumerate(bits,1)}
   accepted=0
   for w in product(sorted(S),repeat=q-1):
    path=[0]
    for d in w:path.append((path[-1]+d)%q)
    distinct=len(set(path))==q;local=all((v-u)%q in S for u,v in zip(path,path[1:]))
    need(local,'word local validity');accepted+=distinct;word_checks+=1
   h=0
   for P in permutations(range(q)):
    h+=all((v-u)%q in S for u,v in zip(P,P[1:]));path_checks+=1
   need(h==q*accepted,'circulant translation bijection')
   units=sum(gcd(d,q)==1 for d in S)
   need(h>=q*units,'constant unit paths')
   circles.append({'q':q,'S':sorted(S),'H':h,'words_starting_0':accepted,'constant_unit_words':units})
 q=9;S={1,2,3,4};w=[3]*8;path=[0]
 for d in w:path.append((path[-1]+d)%q)
 need(all(x in S for x in w) and len(set(path))==3,'composite collision hostile')
 return {'status':'PASS','square_cycle32':{'root_weight':R,'independent_edge_FRACTRAN_steps':total,'full_named_mark_TM_steps':tm,'scope':'edge checks alone versus full witness verifier; not an equal-work speed comparison'},
  'odd_eighths_in_G6':eighths,'circulant_word_checks':word_checks,'circulant_literal_path_checks':path_checks,
  'circulants':circles,'composite_constant_step_hostile':{'q':q,'S':sorted(S),'word':w,'partial_sums':path},
  'scope':'Elementary named-vertex/word correspondences; no new universal Hamiltonicity, LRC or conjecture-search proof.'}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
