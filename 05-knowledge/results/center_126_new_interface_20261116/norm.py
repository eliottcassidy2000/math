from itertools import product,combinations
from verify import dp,prof,op,spacer
def main():
 E=list(combinations(range(5),2))
 P={r:set() for r in range(8,12)}
 for b in product((0,1),repeat=len(E)):
  t=[0]*5
  for (i,j),s in zip(E,b):
   if not s:i,j=j,i
   t[i]|=1<<j
  for u in range(5):
   h,g,r,s=prof(t,u);e=dp(op(t),u)[-1]
   assert s-r==e and e<=g
   if r+s in P:P[r+s].add((h,e))
 expect={8:{65},9:{105,123,125,153,159,161},10:{125,145,165},11:{193,199,201}}
 for a,T in P.items():
  U={a*a+h*j-e*f for h,e in T for j,f in T}
  assert U==expect[a];print(a,sorted(T),sorted(U))
 assert dp(spacer('11','11010',()))[-1]==127
 # Pell norm recipe identities.
 x,y=1,1
 for _ in range(16):
  assert x*x-3*y*y==-2
  x,y=2*x+3*y,x+2*y
if __name__=='__main__':main()
