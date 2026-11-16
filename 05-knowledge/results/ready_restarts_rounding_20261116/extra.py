from verify import dp,prof,join,adjoin,reorder
def rr(T,y):
 n=len(T); m=(1<<n)-1
 rev=tuple(m^v^(1<<i) for i,v in enumerate(T))
 e=dp(rev,y); h=dp(T)
 return sum(e[s]*h[m^s] for s in range(1,m+1) if s>>y&1)
def main():
 for T,y,h in (((18,4,9,19,6),3,11),((58,4,9,18,38,14),0,29)):
  assert dp(T)[-1]==h and rr(T,y)==h+1
  p=join((0,),T)
  q=adjoin(p,((1<<len(p))-1)^(1<<(y+1)))
  assert dp(q)[-1]==3*h+2
 C=(256,5,9,323,271,287,439,319,6)
 assert prof(C,8)==(125,27,62,89)
 assert prof(reorder(C,(0,1,2,3,4,6,8)),6)==(41,9,20,29)
 for v,H,P in ((6246,10003,30011),(6558,10003,30011),(5272,10029,30089)):
  w=format(v,'b')[::-1].translate(str.maketrans('01','10'))
  T=[]
  for a in range(len(w)+1):
   T.append(sum(1<<b for b in range(a+1,len(w)+1) if a or w[b-1]=='1') | (a>0 and w[a-1]=='0'))
  assert dp(T,0)[-1]==1 and rr(T,0)==v+1
  B=adjoin(T,((1<<len(T))-1)^1)
  assert dp(B)[-1]==H
  assert H==2*(P//6)+1
  p=join((0,),B);q=adjoin(p,((1<<len(p))-1)^(1<<(len(B))))
  assert dp(q)[-1]==P
 print('Verified 11→35 and 29→89 via a different port, and a rooted 41→125 two-vertex lift; also verified rounded ancestors for 30011 and 30089.')
if __name__=='__main__':main()
