#!/usr/bin/env python3
"""Complete normalized fixed-root self-converse ternary search through order 9."""
from verify_minimum_ternary import hp,rank,clone,require
import itertools,json

def normalized(n):
 pairs=list(itertools.combinations(range(n),2));orbits=[];index={}
 for e in pairs:
  if e[1]==e[0]+1 or e in index:continue
  group={e,(n-1-e[1],n-1-e[0])}
  for x in group:index[x]=len(orbits)
  orbits.append(group)
 require(len(orbits)==(n-1)**2//4,'complete independent-bit universe')
 for mask in range(1<<len(orbits)):
  rows=[0]*n
  for i,j in pairs:
   if j==i+1 or mask>>index[i,j]&1:rows[i]|=1<<j
   else:rows[j]|=1<<i
  yield rows

def main():
 expected={1:{0:1},3:{0:1,1:1},5:{0:1,1:2,2:13},7:{0:2,1:6,2:154,3:350},9:{0:8,1:32,2:2544,3:10016,4:52936}}
 records=[];rank2=total=first=second=0
 for n in [1,3,5,7,9]:
  histogram={};one=two=0;root=n//2;tau=list(reversed(range(n)))
  for rows in normalized(n):
   total+=1
   require(all(((rows[i]>>j)&1)==((rows[tau[j]]>>tau[i])&1) for i in range(n) for j in range(i+1,n)),'physical anti-involution')
   k=rank(rows,root);histogram[k]=histogram.get(k,0)+1
   if k!=2:continue
   rank2+=1;H1=hp(rows);H2=hp(clone(rows,root,2))
   require(H2%4==1,'even clone congruence')
   if H2!=3*H1+2:continue
   one+=1;first+=1
   if (H1+1)%3:continue
   A=(H1+1)//3;require(A%4==2,'symmetric pure-clock candidate residue')
   H3=hp(clone(rows,root,3));a2=(H2-H1)//2;a3=(H3-H1-6*a2)//6
   require((H3-H1-6*a2)%6==0,'integral third run coefficient')
   if H3==9*H1+8:two+=1;second+=1
  require(histogram==expected[n],'complete matching-rank histogram')
  records.append({'order':n,'normalized_graphs':sum(histogram.values()),'rank_histogram':histogram,'first_equation_survivors':one,'pure_clock_survivors':two})
 require(second==0,'no fixed-root symmetric ternary clock through9')
 require(9*5**2>221 and 5**11>221**3,'odd strong source cutoff for A<=74')
 print(json.dumps({'status':'PASS','orders':records,'total_normalized_sources':total,'rank2_sources':rank2,'first_equation_survivors':first,'pure_clock_survivors':second,'all_order_amplitude_residue_mod4':2,'all_order_amplitude_lower_bound':78,'scope':'Self-converse source equipped with an anti-involution fixing the cloned root; no attainment at78 claimed.'},indent=2,sort_keys=True))
if __name__=='__main__':main()
