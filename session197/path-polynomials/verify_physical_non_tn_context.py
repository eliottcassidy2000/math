#!/usr/bin/env python3
from verify_extension import *

def run():
 root=cyclic((0,),transitive(3),transitive(2));padded=join(root,(0,));validate(root);validate(padded)
 require(root==(14,60,56,48,33,1) and padded==(78,124,120,112,97,65,0),'physical row lists')
 op=skeleton_operator(root,0);pop=skeleton_operator(padded,0)
 expected={1:[0,0,10,0,12,0,1],3:[0,0,232,0,240,0,36,0,1],5:[0,0,1580,0,2862,0,830,0,60,0,1]};profiles={}
 for m in range(7):
  graph=cyclic(transitive(m),transitive(3),transitive(2));p=bin_polynomial(graph);profiles[m]=p
  require(p==apply_skeleton(op,mon(m)),'actual clone profile versus literal skeleton')
  require(mul([0,1],p)==apply_skeleton(pop,mon(m)),'padded root skeleton')
  if m in expected:require(p==expected[m],'exact displayed physical monomial column')
 actualq=[profiles[1],add(profiles[3],scale(profiles[1],-1)),add(profiles[5],scale(profiles[3],-5),scale(profiles[1],4))]
 want=[[23,17,1],[486,424,49,1],[2880,6220,1571,85,1]];cols=[];pcols=[]
 for k in range(3):
  require(actualq[k]==apply_skeleton(op,basis(k,1)),'virtual column from physical monomials')
  c=coordinates(actualq[k],0);require(c==want[k],'exact positive column');cols.append(c)
  pcols.append(coordinates(mul([0,1],actualq[k]),1))
 M=[[cols[j][i] for j in range(3)] for i in range(3)];PM=[[pcols[j][i] for j in range(3)] for i in range(3)]
 require(M==[[23,486,2880],[17,424,6220],[1,49,1571]] and det(M)==-468310,'negative odd-to-even minor')
 require(PM==[[23,486,2880],[91,2182,27760],[26,865,20359]] and det(PM)==-16859160,'negative parity-preserving minor')
 require(det(PM)==36*det(M),'exact bidiagonal determinant factor')
 print(json.dumps({'status':'PASS','root_rows':root,'padded_rows':padded,'actual_clone_profiles':profiles,'odd_to_even_minor':M,'determinant':det(M),'parity_preserving_minor':PM,'padded_determinant':det(PM)},sort_keys=True,indent=2))
if __name__=='__main__':run()
