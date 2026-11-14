from itertools import product

def transitions(d):
 for r,p,q,s in product(range(2),range(3),range(3),range(2)):
  a,b=(d>>(1-r))&1,(d>>r)&1
  i=18*r+6*p+2*q+s
  yield i,18*r+6*(a+1)+2*q+s
  if q or a:yield i,18*r+6*p+2*(b+1)+s
  if s or b:yield i,18*r+6*p+2*q+1
E=list(map(list,map(transitions,range(4))))
def push(v,d):
 u=[0]*36
 for i,j in E[int(d)]:u[j]+=v[i]
 return u
def state(w):
 v=[0]*36
 v[0]=v[18]=1
 for d in w:v=push(v,d)
 return v
def H(w,e=1):
 v=state(w)
 return sum(v[18*r+6*a+2+c]+(v[18*r+6*a+c] if e!=r else 0) for r,a,c in product(range(2),range(2),range(2)))
if __name__=='__main__':
 for n in range(1,6):
  print(n,[( ''.join(w),H('11'+''.join(w)+'311')) for w in product('02',repeat=n)])
