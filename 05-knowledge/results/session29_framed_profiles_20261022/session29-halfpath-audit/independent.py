from itertools import permutations,product,combinations
from functools import lru_cache
from fractions import Fraction
from math import factorial
import json

def require(x,m):
 if not x:raise RuntimeError(m)
def matrix(masks):return tuple(tuple(bool(row>>j&1) for j in range(len(masks))) for row in masks)
def autom(T,p):return all(T[i][j]==T[p[i]][p[j]] for i in range(len(T)) for j in range(len(T)))
def anti(T,p):return all(T[i][j]==T[p[j]][p[i]] for i in range(len(T)) for j in range(len(T)))
def allpaths(T):
 n=len(T)
 def walk(path,used):
  if len(path)==n:yield path;return
  for j in range(n):
   if not used>>j&1 and (not path or T[path[-1]][j]):yield from walk(path+(j,),used|1<<j)
 yield from walk((),0)
def half(T,s):
 n=len(T);pair=[i for i in range(n) if s[i]!=i];fixed=[i for i in range(n) if s[i]==i];m=len(pair)//2
 if not m:return 1,0
 @lru_cache(None)
 def count(mask,last,rooted):
  if mask==1<<last:return int(not rooted or T[fixed[0]][last])
  prev=mask^(1<<last)
  return sum(count(prev,j,rooted) for j in pair if prev>>j&1 and T[j][last])
 hp=cy=0
 for mask in range(1,1<<n):
  if any(mask>>v&1 for v in fixed) or any(mask>>i&1 and mask>>s[i]&1 for i in pair):continue
  size=mask.bit_count()
  for x in pair:
   if not mask>>x&1:continue
   if size==m and T[x][fixed[0] if fixed else s[x]]:hp+=count(mask,x,False)
   if fixed and T[x][s[x]]:cy+=count(mask,x,True)
 return hp,cy
def cycles_fixed(T,s):
 n=len(T);answer=0
 def normalized(p):
  j=p.index(min(p));return p[j:]+p[:j]
 for start in range(n):
  def walk(path,used):
   nonlocal answer
   if len(path)>=3 and len(path)%2 and T[path[-1]][start]:
    answer+=normalized(tuple(s[v] for v in path[::-1]))==path
   for j in range(start+1,n):
    if not used>>j&1 and T[path[-1]][j]:walk(path+(j,),used|1<<j)
  walk((start,),1<<start)
 return answer
def inspect(T,s=None):
 n=len(T);autos=[p for p in permutations(range(n)) if autom(T,p)]
 if s is None:s=next(p for p in permutations(range(n)) if all(p[p[i]]==i for i in range(n)) and anti(T,p))
 paths=list(allpaths(T));f=sum(p==tuple(s[v] for v in p[::-1]) for p in paths)
 tiles={tuple(tuple(T[p[i]][p[j]] for j in range(n)) for i in range(n)) for p in paths}
 b=sum(anti(U,tuple(range(n-1,-1,-1))) for U in tiles)
 c=sum(all(a[s[i]]==s[a[i]] for i in range(n)) for a in autos)
 hf,cf=half(T,s)
 require(hf==f and cf==cycles_fixed(T,s),'independent fixed counts')
 require(f==c*b and len(paths)==len(autos)*len(tiles),'quotients')
 require(f%2==1,'parity')
 if n>=2:require(f<=2**(n//2-1)*factorial(n//2),'bound')
 return {'H':len(paths),'aut':len(autos),'centralizer':c,'fixed_HP':f,'blue':b,'tiles':len(tiles),'fixed_cycles':cf},min(tiles)
def invariant(n):
 s=tuple(range(n-1,-1,-1));edges=list(combinations(range(n),2));orbits=[];done=set()
 for e in edges:
  if e in done:continue
  orb={e,(s[e[1]],s[e[0]])};orbits.append(orb);done|=orb
 for bits in product((False,True),repeat=len(orbits)):
  T=[[False]*n for _ in range(n)]
  for bit,orb in zip(bits,orbits):
   for a,b in orb:T[a][b]=bit;T[b][a]=not bit
  yield tuple(map(tuple,T)),s

def main():
 named={};data={'C3':matrix([2,4,1]),'C5':tuple(tuple((j-i)%5 in (1,2) for j in range(5)) for i in range(5)),'P7':tuple(tuple((j-i)%7 in (1,2,4) for j in range(7)) for i in range(7)),'A':matrix([16,9,3,5,14]),'B':matrix([18,8,3,5,14])}
 for name,T in data.items():named[name]=inspect(T)[0]
 require(named['P7']=={'H':189,'aut':21,'centralizer':3,'fixed_HP':9,'blue':3,'tiles':9,'fixed_cycles':12},'P7')
 census=[];controls=0
 for n in range(1,7):
  classes={};total=0;size=0;signed_to_spine={};orbit_sizes={}
  for T,s in invariant(n):
   row,key=inspect(T,s);classes[key]=row;total+=row['fixed_HP'];size+=1;controls+=1
   commuting=[p for p in permutations(range(n)) if all(p[s[i]]==s[p[i]] for i in range(n))]
   pairkey=min(tuple(tuple(T[p[i]][p[j]] for j in range(n)) for i in range(n)) for p in commuting)
   require(pairkey not in signed_to_spine or signed_to_spine[pairkey]==key,'signed pair isomorphism')
   signed_to_spine[pairkey]=key;orbit_sizes[key]=orbit_sizes.get(key,0)+1
  m=n//2;o=m*m+(m if n%2 else 0)
  require(total==size*factorial(m),'mean')
  require(len(signed_to_spine)==len(classes),'signed pair class count')
  require(all(orbit_sizes[k]*r['centralizer']==2**m*factorial(m) for k,r in classes.items()),'class orbit sizes')
  mass=sum(Fraction(1,r['centralizer']) for r in classes.values());blue=sum(r['blue'] for r in classes.values())
  require(mass==Fraction(2**(o-m),factorial(m)) and blue==2**(o-m),'census')
  census.append({'n':n,'labelled':size,'classes':len(classes),'mean':total//size,'centralizer_mass':str(mass),'blue':blue})
 print(json.dumps({'status':'PASS','independent_full_controls':controls,'named':named,'complete_census_n1_to_6':census,'methods':'Bitmask recursion with partner exclusion; directed path/cycle DFS; full vertex-permutation automorphisms; directly normalized spine tiles.'},indent=2))
if __name__=='__main__':main()
