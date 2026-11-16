from math import isqrt
def rc(w):return w[::-1].translate(str.maketrans('01','10'))
def h(w):
 a,x,p=1,0,0
 for b in map(int,w):
  if b==0:a+=x+p
  p=b;x=2*x+b
 return a
def wbits(n):return format(n,'b')
def tags(P):
 for m in range(2,max(2,P.bit_length()-4)):
  yield '0'+'1'*(m-2)+'0'
  if m==3:
   yield '001'
   yield '100'
def inverse(P):
 ans=[]
 for W in tags(P):
  L=1<<(4+len(W))
  x=isqrt(P//L)
  if not x:continue
  u=wbits(x);v=rc(u)
  l,h1=u+'01'+W+'10'+v,u+'10'+W+'01'+v
  assert h(h1)==h(l)+2
  if P in (h(l),h(h1)):ans.append((x,W,h(l)))
 x=isqrt(P//16)
 if x>0 and x%2==0:x-=1
 if x>0:
  u=wbits(x);v=rc(wbits(x+1))
  a,b=u+'0110'+v,u+'1001'+v
  assert h(b)==h(a)+2
  if P in (h(a),h(b)):ans.append((x,'',h(a)))
 return ans
def main():
 for p in (10007,10009,12539,12541):
  a=inverse(p);assert a;print(p,a)
 for p in (30011,30013):assert not inverse(p)
 for m in range(2,15):
  W='0'+'1'*(m-2)+'0'
  for x in (1,2,5,17,123456789):
   u=wbits(x);v=rc(u)
   for a in (u+'01'+W+'10'+v,u+'10'+W+'01'+v):
    assert (x,W,h(u+'01'+W+'10'+v)) in inverse(h(a))
if __name__=='__main__':main()
