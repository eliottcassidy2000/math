#!/usr/bin/env python3
from verify import *
from fractions import Fraction

A=[16,9,3,5,14]
B=[24,17,3,6,12]
C=[36,1,2,7,15,30]
D=[44,1,2,6,15,30]
TRI=[4,1,2]
FIVE=[8,1,3,6]


def aut_order(rows):
 return sum(induced(rows,p)==rows for p in itertools.permutations(range(len(rows))))


def ordered_join(blocks):
 n=sum(map(len,blocks));rows=[];offset=0
 for block in blocks:
  later=((1<<n)-1)^((1<<(offset+len(block)))-1)
  rows.extend((x<<offset)|later for x in block);offset+=len(block)
 return rows


def literal_runs(rows,r):
 exterior=[x for x in range(len(rows)) if x!=r];out=[0]*(len(rows)+1)
 for perm in itertools.permutations(exterior):
  for selected in range(1,1<<(len(exterior)+1)):
   word=[]
   for gap in range(len(exterior)+1):
    if selected>>gap&1:word.append(r)
    if gap<len(exterior):word.append(perm[gap])
   if all(rows[x]>>y&1 for x,y in zip(word,word[1:])):
    out[selected.bit_count()]+=1
 return out


def matching(rows,r):
 O=[x for x in range(len(rows)) if rows[r]>>x&1]
 I=[x for x in range(len(rows)) if rows[x]>>r&1]
 states={0:0}
 for x in O:
  new=states.copy()
  for used,card in states.items():
   for j,y in enumerate(I):
    if not used>>j&1 and rows[x]>>y&1:
     z=used|(1<<j);new[z]=max(new.get(z,0),card+1)
  states=new
 return max(states.values())


def main():
 labelled_counts={n:0 for n in range(1,7)}
 strong_counts={n:0 for n in range(1,7)}
 seen=set();found=[];strong6_min=10**100
 small={3:set(),5:set()}
 for n in range(1,7):
  for rows in labelled(n):
   h=profiles(rows)[0][0][0]
   strong=len(components(rows))==1
   if n==6 and strong:strong6_min=min(strong6_min,h)
   if n<=5 and strong and h in small:small[h].add(canonical(rows))
   if h!=15:continue
   labelled_counts[n]+=1
   if not strong:continue
   strong_counts[n]+=1
   if tuple(rows) not in seen:
    found.append(rows)
    for perm in itertools.permutations(range(n)):seen.add(tuple(induced(rows,perm)))
 expected={canonical(rows) for rows in [A,B,C,D]}
 require({canonical(rows) for rows in found}==expected,'all strong H15 types')
 require(strong6_min==15,'order6 strong floor')
 require(labelled_counts=={1:0,2:0,3:0,4:0,5:64,6:2208},'all labelled H15 small counts')
 require(strong_counts=={1:0,2:0,3:0,4:0,5:64,6:1440},'strong labelled H15 counts')
 require(small[3]=={canonical(TRI)} and small[5]=={canonical(FIVE)},'unique H3/H5 types')
 autos=[aut_order(rows) for rows in [A,B,C,D]]
 require(autos==[3,5,1,1],'strong automorphism orders')
 require(aut_order(TRI)==3 and aut_order(FIVE)==1,'product automorphism orders')

 formulas=[];constructed=0
 for n in range(5,10):
  items=[]
  for core,aut in zip([A,B,C,D],autos):
   if len(core)>n:continue
   for before in range(n-len(core)+1):
    after=n-len(core)-before
    items.append(([[0]]*before+[core]+[[0]]*after,aut))
  if n>=7:
   for left,right in [(TRI,FIVE),(FIVE,TRI)]:
    for a in range(n-7+1):
     for b in range(n-7-a+1):
      c=n-7-a-b
      items.append(([[0]]*a+[left]+[[0]]*b+[right]+[[0]]*c,3))
  keys=set();mass=Fraction(0)
  for blocks,aut in items:
   key=tuple(canonical(x) for x in blocks)
   require(key not in keys,'duplicate SCC type');keys.add(key)
   T=ordered_join(blocks)
   require(profiles(T)[0][0][0]==15,'constructed SCC count')
   require(len(components(T))==len(blocks),'constructed SCC list')
   if n<=6:require(aut_order(T)==aut,'full automorphism padding check')
   mass+=Fraction(1,aut);constructed+=1
  unlabelled=(n-3)*(n-4)
  number=math.factorial(n)*(5*n*n-17*n-32)//15
  require(len(items)==unlabelled,'unlabelled formula')
  require(math.factorial(n)*mass==number,'labelled formula')
  if n<=6:require(number==labelled_counts[n],'independent labelled formula control')
  formulas.append({'order':n,'unlabelled':unlabelled,'labelled':number})

 functions={(0,0,15)};details=[];clone_checks=0
 for rows in [A,B,C,D]:
  roots=[]
  for r in range(len(rows)):
   a=literal_runs(rows,r);nu=matching(rows,r)
   require(nu<=2 and max(k for k,x in enumerate(a) if x)==1+nu,'matching/run rank')
   require(a[1]==15 and not any(a[4:]),'literal coefficient support')
   coeff=(a[2]-3*a[3],a[3],15-2*a[2]+3*a[3]);functions.add(coeff)
   values=[]
   for m in range(1,6):
    h=profiles(clone(rows,r,m))[0][0][0]
    require(h==coeff[0]*(1<<m)+coeff[1]*3**m+coeff[2],'full clone expansion')
    values.append(h);clone_checks+=1
   roots.append({'root':r,'run_coefficients':a[1:4],'matching_number':nu,
                 'exponential_coefficients':coeff,'actual_counts_m1_to5':values})
  details.append({'rows':rows,'automorphism_order':aut_order(rows),'roots':roots})
 expected_functions={(0,0,15),(9,0,-3),(5,0,5),(6,0,3),(11,0,-7),
                     (5,2,-1),(7,2,-5),(4,6,-11)}
 require(functions==expected_functions,'complete clone function list')
 for small,expected_small in [(TRI,{(1,1)}),(FIVE,{(3,-1),(2,1)})]:
  got=set()
  for r in range(len(small)):
   run=literal_runs(small,r)
   require(matching(small,r)==1 and not any(run[3:]),'small component binary rank')
   got.add((run[2],run[1]-2*run[2]))
  require(got==expected_small,'small component exact clone functions')
 for blocks in [[TRI,FIVE],[FIVE,TRI]]:
  T=ordered_join(blocks)
  for r in range(len(T)):
   vals=[profiles(clone(T,r,m))[0][0][0] for m in range(1,4)]
   require(any(vals==[a*(1<<m)+b*3**m+c for m in range(1,4)]
               for a,b,c in expected_functions),'weak product clone family')
   clone_checks+=3
 T=ordered_join([A,[0]])
 require(profiles(clone(T,5,3))[0][0][0]==15,'constant clone family')
 possible_two=sorted({a*4+b*9+c for a,b,c in functions})
 require(possible_two==[15,25,27,33,37,41,59],'two-clone scalar spectrum')
 require(all(a*8+b*27+c>=45 for a,b,c in functions if (a,b)!=(0,0)),
         'later clone lower bound')
 print(json.dumps({'status':'PASS','strong_order6_minimum':strong6_min,
                   'labelled_H15_counts':labelled_counts,'strong_labelled_H15_counts':strong_counts,
                   'strong_types_and_clone_data':details,'SCC_count_formula_checks':formulas,
                   'constructed_padding_classes':constructed,'actual_clone_graph_checks':clone_checks+1,
                   'complete_clone_exponential_coefficients':sorted(functions),
                   'complete_two_clone_count_set':possible_two,
                   'scope':'Complete H15 stratum and transitive cloning responses; full spectrum remains open.'},
                  indent=2,sort_keys=True))


if __name__=='__main__':main()
