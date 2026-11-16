from itertools import combinations
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
 print('Physical 125/127 +2 windows, both centered trims, slope-100 lower host, spacer family, and literal 100+27 path pools verified.')
if __name__=='__main__':main()
