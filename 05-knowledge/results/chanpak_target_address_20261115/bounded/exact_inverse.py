from placement import h_word,reverse_complement

def exact_inverse(N):
 if N==1:return (1,0)
 if N==3:return (2,2)
 if N<=0 or not N&1:return None
 k=N.bit_length()-1
 for m in range(max(3,k),k+3):
  r,e=divmod(m,2);L=1<<r;J=r+1
  for a in range(L//2,L):
   lo=max(L//2,(N-J)//((a+1)<<e))
   hi=min(L-1,(N+J)//(a<<e))
   for b in range(lo,hi+1):
    suffix=reverse_complement(b,r)
    for c in range(1<<e):
     w=(a<<(r+e))|(c<<r)|suffix
     if h_word(w,m)==N:return m,w
 return None

if __name__=='__main__':
 import sys
 if len(sys.argv)>1:
  for N in map(int,sys.argv[1:]):
   result=exact_inverse(N)
   print(N,'no one-apex witness' if result is None else f'{result[1]:0{result[0]}b}')
 else:
  brute={1,3}
  for m in range(3,12):
   brute.update(h_word((1<<(m-1))|(a<<1),m) for a in range(1<<(m-2)))
  for N in range(1,1024,2):
   result=exact_inverse(N)
   assert (result is not None)==(N in brute)
   if result:
    m,w=result
    assert h_word(w,m)==N
    if N>1:assert w>>(m-1)==1 and not w&1
  from random import Random
  rng=Random(57742)
  for m in range(12,24):
   for i in range(10):
    w=(1<<(m-1))|(rng.randrange(1<<(m-2))<<1)
    N=h_word(w,m)
    length,wit=exact_inverse(N)
    assert h_word(wit,length)==N
  print('bounded target decisions:',512,'additional larger reconstructions:',120)
