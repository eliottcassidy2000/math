from atlas import canon, subset, op
from program import prof,dp,diff
from itertools import permutations
def main():
 a=(6,12,56,97,11,83,23)
 b=(10,20,33,70,77,27,39)
 assert canon(a)==canon(op(a)) and canon(b)==canon(op(b))
 assert subset(b,(6,4,5,3,1,2,0))==op(b)
 l=(1<<7)-1
 H,B=dp(a),dp(b)
 assert sorted(H[l^(1<<i)] for i in range(7))==[19,19,37,37,41,43,43]
 assert sorted(B[l^(1<<i)] for i in range(7))==[29,29,29,31,31,37,37]
 A37={canon(subset(a,tuple(i for i in range(7) if i!=j))) for j in range(7) if H[l^(1<<j)]==37}
 B37={canon(subset(b,tuple(i for i in range(7) if i!=j))) for j in range(7) if B[l^(1<<j)]==37}
 assert A37==B37 and len(A37)==2
 raw=[];locked=[]
 for p in permutations(range(7)):
  R=subset(b,p);E=diff(a,R)
  raw.append(len(E))
  if E and len(set(E[0]).intersection(*(set(v) for v in E))):
   for i in set(E[0]).intersection(*(set(v) for v in E)):
    if H[l^(1<<i)]==37:locked.append(len(E))
 assert min(raw)==3 and min(locked)==4
 print('Both H37 parents verified; edit minima 3 unrestricted, 4 preserving a shared H37 parent.')
if __name__=='__main__':main()
