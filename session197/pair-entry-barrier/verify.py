from pathlib import Path
from itertools import combinations
import json
from core import *
P=Path(__file__).resolve().parent

def main():
 bank=json.loads((P/'control_graphs.json').read_text());checked=0;excesses=set()
 for n in range(3,7):
  for z in bank[str(n)]:
   T=z['rows']
   if cycles(T)!=1:continue
   p=inventories(T);h=p['H']
   for a in range(n):
    if T[a].bit_count()!=1:continue
    b=T[a].bit_length()-1;ea,eb,Ra,Rb=p['e'][a],p['e'][b],p['R'][a],p['R'][b]
    require(Rb==2*(h-ea)-eb,'adjacent response identity')
    U,ports=birth(T,[a,b],3,1);f=n;q=inventories(U)
    excess=q['R'][f]-2*(q['H']-q['e'][f]);require(excess==Ra-h+eb>=1,'first clone outside cone')
    require(q['R'][a]<=2*(q['H']-q['e'][a]),'last clone inside cone');excesses.add(excess);checked+=1
 strong=[]
 for n in (3,4):
  pairs=list(combinations(range(n),2));values=set();count=0
  for bits in range(1<<len(pairs)):
   T=[0]*n
   for j,(u,v) in enumerate(pairs):
    if bits>>j&1:u,v=v,u
    T[u]|=1<<v
   if cycles(T):values.add(literal_h(T));count+=1
  require(values==({3} if n==3 else {5}),'small strong H classification');strong.append({'order':n,'labelled_strong':count,'H_values':sorted(values)})
 T=[18,20,25,3,8];require(root_profile(T,0)==(13,2,18) and inventories(T)['s'][0]==2 and cycles(T)==2,'inside-cone witness')
 require(5**4>7**3,'H7 strong cutoff')
 return {'status':'PASS','physical_first_last_clone_pairs':checked,'positive_excesses':sorted(excesses),'small_strong_checks':strong,'inside_cone_profile':[13,2,18],'impossible_singleton_inverse':[7,5,11]}
if __name__=='__main__':print(json.dumps(main(),indent=2,sort_keys=True))
