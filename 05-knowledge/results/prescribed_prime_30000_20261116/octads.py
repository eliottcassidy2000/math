from verify import count,flip,h_word
from math import isqrt

def dot(x,y):return sum(a*b for a,b in zip(x,y))
def scale(a,x):return tuple(a*b for b in x)
def add(x,y):return tuple(a+b for a,b in zip(x,y))
def indicator(s):return tuple(int(i in s) for i in range(1,25))
G={0}
for v in [7973<<j for j in range(11)]+[(1<<24)-1]:
 G|={a^v for a in G}
def membership(x):
 p=x[0]&1
 return (all(k%2==p for k in x) and sum(x)%8==4*p and
         sum(1<<i for i,k in enumerate(x) if k%4==p+2) in G)
def signed(pos,plus):
 return tuple((2 if i in plus else -2) if i in pos else 0 for i in range(1,25))
def main():
 A={1,3,6,9,10,11,12,13}
 B={5,7,8,14,15,18,19,20}
 D={2,4,16,17,21,22,23,24}
 R={2,3,4,5,7,8,10,15}
 S={1,2,3,4,6,7,9,14}
 assert len(G)==4096 and all(sum(1<<(i-1) for i in t) in G for t in (A,B,D,R,S))
 assert A|B|D==set(range(1,25)) and not A&B and not A&D and not B&D
 x,y,z=(scale(100,indicator(t)) for t in (A,B,D))
 s=add(add(x,y),z)
 a,b,r,q=(signed(t,p) for t,p in ((A,{1,6,10,11}),(B,{5,7,8,14}),(R,{2,4,5,10}),(S,{1,2,4,6})))
 us=[add(a,scale(-1,r)),add(a,b),add(a,r),add(a,q)]
 assert all(membership(t) for t in (x,y,z,s,a,b,r,q,*us))
 assert [dot(t,t)//8 for t in (x,y,z,s)]==[10000,10000,10000,30000]
 assert all(dot(t,u)==0 for t,u in zip((x,x,s,s),us))
 for t,u,n in zip((x,x,s,s),us,(10006,10008,30010,30012)):
  v=add(t,u);assert membership(v) and dot(v,v)==8*n
 w=(1232,389,1977,979,770,1819,1462,560,1665,67,666)
 e,f=(0,6),(3,4)
 assert [count(flip(w,*es)) for es in ((),(e,),(f,),(e,f))]==[30089,23331,38397,30091]
 for t in (30089,30091):assert all(t%d for d in range(2,isqrt(t)+1))
 hits={}
 for n in range(2,17):
  for i in range(1<<(n-2)):
   v='1'+format(i,f'0{n-2}b')+'0' if n>2 else '10'
   h=h_word(v)
   if h in (30089,30091):hits.setdefault(h,set()).add(v)
 assert hits=={30089:{'1001010101101100','1100100101010110'}}
 print('Four explicit Leech norms: 10006, 10008, 30010, 30012.')
 print('Second physical window: 30089/30091; intermediates 23331 and 38397.')
if __name__=='__main__':main()
