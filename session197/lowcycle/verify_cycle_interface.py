#!/usr/bin/env python3
"""Literal and graph-DP certificate for the five-coordinate cycle interface."""
from itertools import permutations
from pathlib import Path
import json
from verify import require
from verify_rooted import profile,glue,induced

def literal_profile(rows,r):
 n=len(rows);h=e=s=c=R=0
 for p in permutations(range(n)):
  valid=all(rows[a]>>b&1 for a,b in zip(p,p[1:]))
  if valid:
   h+=1;e+=p[-1]==r;s+=p[0]==r
   if p[0]==r and n>=3 and rows[p[-1]]>>r&1:c+=1
  for cut in range(n):
   if p[cut]!=r:continue
   if all(rows[p[i]]>>p[i+1]&1 for i in range(n-1) if i!=cut-1):R+=1
 vs=[u for u in range(n) if u!=r];small=induced(rows,vs)
 d=sum(all(small[a]>>b&1 for a,b in zip(p,p[1:])) for p in permutations(range(n-1)))
 return h,e,R,s,c,d

def determinant(a):
 n=len(a);total=0
 for p in permutations(range(n)):
  inversions=sum(p[i]>p[j] for i in range(n) for j in range(i+1,n));v=(-1)**inversions
  for i in range(n):v*=a[i][p[i]]
  total+=v
 return total

def main():
 collision=[([32,17,11,19,37,14],3),([32,17,11,35,13,22],2)]
 reports=[]
 for A,a in collision:
  p=profile(A,a);lp=literal_profile(A,a);require(lp==p[:6],'literal input profile')
  require((p[0],p[1],p[2],p[5],p[4])==(37,5,71,13,4),'collision signature')
  T,r=glue(A,a,[0],0);out=profile(T,r);require(literal_profile(T,r)==out[:6],'literal output profile')
  require(out[:3]==(103,32,108) and out[4]==p[3],'collision continuation')
  reports.append({'input_rows':A,'input_root':a,'input_profile_H_e_R_s_C_d':list(p[:6]),'output_rows':T,'output_root':r,'output_H':out[0],'output_C':out[4]})
 require([x['output_C'] for x in reports]==[9,8],'distinct cycles')
 seed=[2,4,1];inputs=[(seed,0)]
 for G,b in [([0],0),([0,1],0),([0,1,3],0),([8,1,3,6],1)]:inputs.append(glue(seed,0,G,b))
 expected=[(3,1,3,1,1),(5,2,6,1,1),(9,3,9,1,1),(15,3,15,1,1),(41,6,54,3,2)]
 table=[];input_reports=[]
 for (A,a),want in zip(inputs,expected):
  p=profile(A,a);require(p[:5]==want,'five-state input')
  require(literal_profile(A,a)==p[:6],'literal five-state input')
  S,sr=glue(A,a,[0],0);E,er=glue(A,a,[0,1],0);ps=profile(S,sr);pe=profile(E,er)
  table.append([p[0],ps[0],pe[0],p[4],pe[4]])
  input_reports.append({'rows':A,'root':a,'state_H_e_R_s_C':list(want)})
 require(table==[[3,5,9,1,1],[5,9,17,1,1],[9,15,27,1,1],[15,27,45,1,1],[41,89,149,2,3]],'continuation table')
 require(determinant(table)==-24,'rank-five minor')
 out={'collision':reports,'five_inputs':input_reports,'continuation_table':table,'determinant':-24,'joint_linear_dimension':5}
 text=json.dumps(out,sort_keys=True,indent=2)+'\n';(Path(__file__).resolve().parent/'cycle_interface_certificate.json').write_text(text);print(text,end='')
if __name__=='__main__':main()
