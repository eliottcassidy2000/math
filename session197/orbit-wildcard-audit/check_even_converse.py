#!/usr/bin/env python3
from pathlib import Path
from itertools import combinations,permutations
import hashlib,json
from check_graph import check,h_walk_inclusion_exclusion,literal_h
from check_converse_mass import normalize_universe,degree_canonical
from check_converse_spectrum import valid,converse
SOURCE=Path(__file__).resolve().parents[1]/'session197-orbit-wildcard'
HASHES={'EVEN_CONVERSE.md':'0fc0e7ca35a1d48626ffd1b49251b64dcc0564f66594bfd0e70307287c4c4cbb','verify_even_converse.py':'ee8c15a3a27a2b5f1540d17eb04f9e71084fc3fec612659903c53bed2a413e77','even-converse-normal.json':'dd1e6716e9fa29d7cf2764c013c4ac75b2531cfac4098f36da6e479d50edb4df','even-converse-witness.json':'08c74a5182999514685de72edd9c776970aa49324fa1e20816c93f6ceb5cd2cc'}
EXPECTED={1,3,5,9,11,13,15,17,25,27,29,31,33,37,41,45,49,51,57,61,65,67,69,71,73}

def oddcycles(rows):
 n=len(rows)
 for k in range(3,n+1,2):
  for vertices in combinations(range(n),k):
   start=vertices[0]
   for rest in permutations(vertices[1:]):
    word=(start,)+rest
    if all(rows[a]>>b&1 for a,b in zip(word,word[1:]+word[:1])):yield word

def cycle_image(word,tau):
 image=tuple(tau[v] for v in reversed(word));i=image.index(min(image))
 return image[i:]+image[:i]

def main():
 for name,x in HASHES.items():check(hashlib.sha256((SOURCE/name).read_bytes()).hexdigest()==x,('frozen',name))
 support=set();orders=[];cycle_checks=even_checks=0
 for n in range(1,9):
  vals=set();count=0;tau=tuple(reversed(range(n)))
  for rows in normalize_universe(n):
   valid(rows);H=h_walk_inclusion_exclusion(rows)[0];count+=1
   if H<=73:vals.add(H);support.add(H)
   if not(n%2):check(H%4==1,'even all-order congruence control');even_checks+=1
   if n<=7:
    cs=set(oddcycles(rows));fixed=[w for w in cs if cycle_image(w,tau)==w]
    check({cycle_image(w,tau) for w in cs}==cs,'cycle involution closure')
    check(H%4==(1+2*len(cs))%4==(1+2*len(fixed))%4,'independent invariant-cycle congruence')
    if n%2:check(all(n//2 in w for w in fixed),'invariant odd cycles contain unique fixed vertex')
    else:check(not fixed and not(len(cs)%2),'fixed-point-free odd-cycle pairing')
    cycle_checks+=1
  check(count==1<<((n-1)**2//4),'complete normalized count')
  orders.append({'n':n,'presentations':count,'support_through73':sorted(vals)})
 check(support==EXPECTED,'full finite support')
 check(5**8>73**3,'strong central cutoff')
 check({d*d*h for d in [1,3,5,7] for h in support if d*d*h<=73}==support,'all-order square closure')
 witness=json.loads((SOURCE/'even-converse-witness.json').read_text());escape={}
 for value,item in witness['self_converse'].items():
  rows=item['rows'];tau=item['anti_involution'];valid(rows)
  check(h_walk_inclusion_exclusion(rows)[0]==int(value),'physical positive witness')
  check(all(tau[tau[i]]==i for i in range(len(rows))),'witness involution')
  check(all(((rows[i]>>j)&1)==((rows[tau[j]]>>tau[i])&1) for i in range(len(rows)) for j in range(len(rows))),'witness anti arcs')
 for value,item in witness['new_unrestricted_escapes'].items():
  rows=item['rows'];word=item['word'];valid(rows);n=len(rows)
  build=[0]*n;build[0]=2
  for k,c in enumerate(word):
   v=k+2
   for a in [0,1]:
    if c>>a&1:build[a]|=1<<v
    else:build[v]|=1<<a
   for j in range(k):build[j+2]|=1<<v
  check(rows==build,'independent physical two-apex construction')
  H=h_walk_inclusion_exclusion(rows)[0];check(H==int(value)==literal_h(rows),'independent escape counts')
  check(H not in support and degree_canonical(rows)!=degree_canonical(converse(rows)),'actual asymmetric escape')
  escape[value]={'H':H,'order':n}
 check(h_walk_inclusion_exclusion([2,4,1])[0]%4==3,'odd-order hostile')
 check(h_walk_inclusion_exclusion([2,8,3,5])[0]==5,'mod8 hostile')
 print(json.dumps({'status':'PASS','source_hashes':HASHES,'orders':orders,'invariant_cycle_controls':cycle_checks,'even_mod4_controls':even_checks,'all_order_support_through73':sorted(support),'escape_checks':escape},indent=2,sort_keys=True))
if __name__=='__main__':main()
