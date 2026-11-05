#!/usr/bin/env python3
from verify_retained_root_sewing import *

def sink(A):return tuple(x|(1<<len(A))for x in A)+(0,)
def K(A,a):return sew(A,a,(0,),0)
def s195_state(A,a):
 h,e,_,R,_,d,_,_=state(A,a)[0]
 return(h,e,R,d)
def newA(z):
 h,e,R,d=z;return(h+R-e,h-e,h+R,h)
def oldK(z):
 h,e,R,d=z;return(h+R-e,d,2*R-d,d)
def oldD(z):
 h,e,R,d=z;return(h,0,2*R,d)
def run():
 checked=0
 for n in range(1,5):
  for code in range(1<<(n*(n-1)//2)):
   A=tournament(n,code)
   for a in range(n):
    for m in range(1,6):
     for j in range(m):
      B=A
      for _ in range(j):B=sink(B)
      B=K(B,a)
      for _ in range(m-1-j):B=sink(B)
      require(B==sew(A,a,transitive(m),j),'literal inherited KD decomposition');checked+=1
 zero=(1,1,1,1);older=[zero,newA(zero),newA(newA(zero)),oldD(newA(newA(zero)))]
 newer=[zero,oldD(oldK(zero)),oldK(oldD(zero)),oldD(oldK(oldD(oldK(zero))))]
 require(older==newer,'same source input states')
 oldrows=[[z[0]for z in older],[newA(z)[0]for z in older],[newA(newA(z))[0]for z in older],[newA(oldK(z))[0]for z in older]]
 newrows=[[z[0]for z in newer],[oldK(z)[0]for z in newer],[oldK(oldK(z))[0]for z in newer],[oldK(oldD(z))[0]for z in newer]]
 require(newrows==[oldrows[i]for i in(0,1,3,2)]and det(newrows)==16,'source rank-certificate correspondence')
 print(json.dumps({'status':'PASS','literal_transitive_block_identifications':checked,'common_input_states':older,'older_matrix':oldrows,'retained_matrix':newrows,'retained_determinant':16},sort_keys=True,indent=2))
if __name__=='__main__':run()
