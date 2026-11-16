from itertools import combinations,product
def dp(r,root=None):
 n=len(r);N=1<<n
 D=[[0]*n for _ in range(N)]
 for i in range(n) if root is None else (root,):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,h in enumerate(row):
   if not h:continue
   e=r[i]&~s
   while e:
    b=e&-e;e^=b
    D[s|b][b.bit_length()-1]+=h
 return [1]+[sum(row) for row in D[1:]]
def prof(r,x):
 a=dp(r);b=dp(r,x);m=len(a)-1
 z=[sum(b[s]*a[m^s] for s in range(1,m+1) if s.bit_count()%2==i) for i in (0,1)]
 return a[m],a[m^(1<<x)],*z
def adjoin(r,mask):
 n=len(r)
 return tuple(row|((not(mask>>i&1))<<n) for i,row in enumerate(r))+(mask,)
def join(a,b):
 n,m=len(a),len(b)
 return tuple(i|(((1<<m)-1)<<n) for i in a)+tuple(i<<n for i in b)
def tt(n):return tuple((1<<n)-(1<<(i+1)) for i in range(n))
def flip(r,es):
 t=list(r)
 for i,j in es:t[i]^=1<<j;t[j]^=1<<i
 return tuple(t)
def diff(a,b):
 return tuple((i,j) for i,j in combinations(range(len(a)),2) if (a[i]^b[i])>>j&1)
def one(w):
 r=tt(len(w)+1)
 return flip(r,[(0,i) for i,b in enumerate(w,1) if b=='0'])
def op(r):
 n=len(r)
 return tuple(((1<<n)-1)^a^(1<<i) for i,a in enumerate(r))
def spacer(a,b,m):
 # 'b' is converse module's original word. Both marked roots are index 0.
 x,y=one(a),op(one(b))
 r=join(join(x,m),y)
 return flip(r,[(0,len(x)+len(m))])
def restrict(t,vs):
 return tuple(sum(((t[u]>>v)&1)<<i for i,v in enumerate(vs)) for u in vs)
def use(r,i,j):
 A=dp_end(r);B=dp_end(op(r));N=(1<<len(r))-1
 return sum(a[i]*B[N^s][j] for s,a in enumerate(A) if s>>i&1 and not(s>>j&1))
def dp_end(r):
 n=len(r);N=1<<n;D=[[0]*n for _ in range(N)]
 for i in range(n):D[1<<i][i]=1
 for s,row in enumerate(D):
  for i,h in enumerate(row):
   if not h:continue
   e=r[i]&~s
   while e:
    b=e&-e;e^=b;D[s|b][b.bit_length()-1]+=h
 return D
def main():
 a=(4,5,72,35,15,87,155,63);b=(4,5,72,99,47,71,147,63)
 assert diff(a,b)==((3,6),(4,5))
 assert [prof(t,2) for t in (a,flip(a,((3,6),)),flip(a,((4,5),)),b)]==[
  (125,9,60,69),(141,13,66,79),(105,9,50,59),(127,11,60,71)]
 for k in range(5):
  assert prof(a,2)==(5+120*2**k,9,60*2**k,9+60*2**k)
  assert prof(b,2)==(7+120*2**k,11,60*2**k,11+60*2**k)
  assert diff(a,b)==((3,6),(4,5))
  a,b=adjoin(a,4),adjoin(b,4)
 a,b=(1018,1020,1017,272,288,840,920,824,512,24),(1018,1020,1017,784,288,328,920,824,512,48)
 ea=((0,2),(2,5),(2,9),(3,4),(3,9),(6,7))
 eb=((0,2),(2,5),(3,4),(4,7))
 c,d=flip(a,ea),flip(b,eb)
 for k in range(4):
  q=60*2**k
  assert [prof(t,9) for t in (a,b,c,d)]==[(3+2*q,27,q,27+q),(9+2*q,27,q,27+q),
                                          (5+2*q,9,q,9+q),(7+2*q,11,q,11+q)]
  assert diff(a,c)==ea and diff(b,d)==eb
  a,b,c,d=[adjoin(t,1<<9) for t in (a,b,c,d)]
 lo=(1018,1020,1017,272,288,840,920,824,512,24)
 hi=(1018,1020,1017,784,288,328,920,824,512,48)
 c=(1022,1020,952,768,296,968,412,792,512,80)
 d=(1022,1020,952,768,360,840,396,824,512,80)
 assert diff(lo,c)==((0,2),(2,6),(3,4),(3,9),(5,7),(6,9))
 assert diff(hi,d)==((0,2),(2,6),(3,4),(4,6),(5,9),(6,9))
 assert len(diff(lo,hi))==len(diff(c,d))==2
 assert [prof(z,9) for z in (lo,hi,c,d)]==[(123,27,60,87),(129,27,60,87),
                                         (125,9,60,69),(127,11,60,71)]
 h5=(2,4,9,3)
 t=join(join(h5,h5),h5)
 assert prof(t,9)==(125,25,50,75)
 assert prof(adjoin(t,1<<9),9)==(225,25,100,125)
 assert dp(spacer('1001','1001',()))[-1]==125
 for k in range(1,6):
  assert dp(spacer('100','110',tt(k)))[-1]==25+35*2**k
  assert dp(spacer('100010','',tt(k)))[-1]==27+35*2**k
 for t,m in (((4,9,66,21,39,15,187,63),(4,1)),((16,9,11,49,70,23,175,63),(2,0))):
  assert dp(t)[-1]==127 and use(t,*m)==27
 assert prof((256,1,267,51,359,199,399,287,42),8)==(377,33,235,268)
 assert prof((2,4,25,307,99,199,143,31,247),2)==(379,33,235,268)
 aa=(4,5,72,51,39,455,155,31,223);bb=(4,5,72,179,103,455,139,23,223)
 assert diff(aa,bb)==((3,7),(4,6))
 assert prof(aa,2)==(377,23,186,209) and prof(bb,2)==(379,29,186,215)
 T=(2,64,11,19,39,79,157,63)
 assert dp(T)[-1]==127 and dp(restrict(T,(0,1,2,3,5,6,7)))[-1]==27 and use(T,5,3)==27
 for T,rows in (((4,5,72,35,15,87,155,63),[(31,9,17,22),(65,9,30,39)]),
                ((4,5,72,99,47,71,147,63),[(33,11,17,24),(67,11,30,41)]),
                (aa,[(89,23,51,65),(191,23,93,116)]),
                (bb,[(91,29,51,71),(193,29,93,122)])):
  for j,profile in ((2,rows[0]),(1,rows[1])):
   assert T[j-1]==(4 if j==1 else 5)
   u=restrict(T,range(j,len(T)))
   assert prof(u,2-j)==profile
   assert adjoin(u,1<<(2-j))==restrict(T,tuple(range(j,len(T)))+(j-1,))
 Q=(62,124,120,112,32,64,17)
 assert dp(Q)[-1]==27
 u,v=adjoin(Q,54),adjoin(Q,100)
 assert (dp(u)[-1],dp(v)[-1])==(125,127)
 assert diff(u,v)==((1,7),(4,7),(6,7))
 for u,v,w in product((0,1),repeat=3):
  cut=46+46*u-18*v+4*w-25*u*v-6*u*w-5*v*w
  m=(1<<2)|(1<<5)|(u<<1)|(v<<4)|(w<<6)
  assert dp(adjoin(Q,m))[-1]==27+2*cut
 print('Physical 125/127 +2 windows, both centered trims, slope-100 lower host, spacer family, and literal 100+27 path pools verified.')
if __name__=='__main__':main()
