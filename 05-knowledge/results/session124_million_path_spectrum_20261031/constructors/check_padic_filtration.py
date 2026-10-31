from itertools import product
import random
rng=random.Random(124124)
def req(x,m):
 if not x:raise RuntimeError(m)
for p in (2,3,5,7):
 for trial in range(1000):
  F,q,g,h=[rng.randrange(1,50)for _ in range(4)];r=0
  for j in range(30):
   H,E,a,b,c,d,u,v=[rng.randrange(10)for _ in range(8)]
   if j%3!=0:H=p*H;E=p*E
   singular=H%p==0 and E%p==0
   # General integral form of the four-state transfer. The actual form
   # specializes a=p_B,b=H-s,c=Fb-H,d=qB,u=H-tB;v is unused.
   F,q,g,h=H*F+a*q+b*g+c*h,E*q+d*h,u*h,H*h
   r+=singular
   req(h%p**r==0,('h',p,trial,j,r))
   req(q%p**max(0,r-1)==0 and g%p**max(0,r-1)==0,('qg',p,trial,j,r))
   req(F%p**max(0,r-2)==0,('F',p,trial,j,r))
print('p_adic_filtration_120000_general_integral_matrix_steps_PASS')
print('NOTE_finite_hostile_audit_not_substitute_for_inductive_proof')
