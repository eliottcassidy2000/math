from itertools import product

def need(x):
 if not x:raise AssertionError
def params(w):
 A,B=1,0
 for j,c in enumerate(w):
  if c:A,B=3*A,3*B+2**j
 return A,B

def main():
 need(pow(5,32,256)==129 and pow(5,64,256)==1)
 need(pow(5,32,641)==pow(2,32,641)==640 and pow(519,4,641)==10 and pow(519,64,641)==640)
 vals=[(s+v,v,v) for s in range(1,8) for v in range(1,4) if 2**s==5**v+3]
 need(vals==[(4,1,1),(10,3,3)])
 for k in range(1,8):
  for w in product((0,1),repeat=k):
   A,B=params(w)
   for N in range(1,11):
    x,z=divmod((N<<k)-B,A)
    if z:continue
    need(x>0)
    for c in w:
     need(x%2==c);x=(3*x+1)//2 if c else x//2
    need(x==N)
 print('Fixed-offset modular certificates and small exact endpoint decoders verified.')
if __name__=='__main__':main()
