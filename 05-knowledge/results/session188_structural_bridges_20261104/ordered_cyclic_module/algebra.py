"""Integer polynomial identity and squarefree-exterior positive functional."""
from math import comb

def symbolic_identity():
 k=10;zero=(0,)*k
 def var(j):a=[0]*k;a[j]=1;return {tuple(a):1}
 def add(*ps):
  out={}
  for p in ps:
   for e,c in p.items():out[e]=out.get(e,0)+c
  return {e:c for e,c in out.items() if c}
 def mul(*ps):
  out={zero:1}
  for p in ps:
   z={}
   for e,c in out.items():
    for f,d in p.items():g=tuple(a+b for a,b in zip(e,f));z[g]=z.get(g,0)+c*d
   out={e:c for e,c in z.items() if c}
  return out
 def sc(c,p):return {e:c*d for e,d in p.items() if c*d}
 one={zero:1};x,y,sg,a,b,al,be,mu,ka,ta=map(var,range(k));h=add(x,y,mul(x,y));D=add(one,sc(-1,mul(sg,h)));U=add(mul(mu,h),x,mul(x,y));V=add(mul(ka,h),y,mul(x,y))
 Nxx=mul(x,add(one,sc(-1,mul(sg,y))));Nyy=mul(y,add(one,sc(-1,mul(sg,x))));Nxy=mul(x,y,add(one,sg));Nyx=mul(x,y,sg);Cx=x;Cy=mul(y,add(one,x));Rx=mul(x,add(one,y));Ry=y
 A=add(mul(al,D),mul(a,add(Cy,mul(ka,h))));B=add(mul(be,D),mul(b,add(Rx,mul(mu,h))));T=add(mul(ta,D),Nxy,mul(mu,Cy),mul(ka,Rx),mul(mu,ka,h));Ex=mul(a,Cx);Ey=mul(a,Cy);Sx=mul(b,Rx);Sy=mul(b,Ry);rx=add(Nxx,mul(mu,Cx));ry=add(Nxy,mul(mu,Cy));lx=add(Nxy,mul(ka,Rx));ly=add(Nyy,mul(ka,Ry))
 num=add(sc(2,mul(A,B,add(Nxy,sc(-1,Nyx)))),sc(-1,mul(add(D,sc(2,T)),add(mul(Ex,Sy),sc(-1,mul(Ey,Sx))))),mul(A,add(mul(rx,Sy),sc(-1,mul(ry,Sx)))),mul(B,add(mul(Ex,ly),sc(-1,mul(Ey,lx)))))
 L=add(sc(2,mul(al,be)),mul(al,b),mul(a,be));J=add(mul(al,b,U),mul(a,be,V),mul(a,b,add(mul(add(ka,mu,sc(2,ta)),h),sc(2,x),sc(2,y),mul(x,y))))
 pos=mul(x,y,add(mul(L,D,D),mul(J,D),sc(2,mul(a,b,U,V))))
 if num!=pos:raise RuntimeError('exact positive rational identity')
 return {'variables':['x','y','sigma','a','b','alpha','beta','mu','kappa','tau'],'numerator_monomials':len(num),'complete_integer_identity':True}

def positive(E,Aout,Mout,bx,by):
 from kernel import matrix
 n=len(E);N=1<<n;nx=len(bx)-1;ny=len(by)-1;zero=(0,0,0);one={zero:1}
 def varx():return {(0,1,0):1}
 def vary():return {(0,0,1):1}
 def add(*ps):
  out={}
  for p in ps:
   for e,c in p.items():out[e]=out.get(e,0)+c
  return {e:c for e,c in out.items() if c}
 def sc(c,p):return {e:c*d for e,d in p.items() if c*d}
 def mul(*ps):
  out=one
  for p in ps:
   z={}
   for (M,i,j),a in out.items():
    for (S,u,v),b in p.items():
     if M&S or i+u>nx or j+v>ny:continue
     key=(M|S,i+u,j+v);z[key]=z.get(key,0)+a*b
   out={e:c for e,c in z.items() if c}
  return out
 G=matrix(E)
 def paths(starts,ends):return {(mask,0,0):sum(G[u][v][mask] for u in starts for v in ends) for mask in range(1,N) if sum(G[u][v][mask] for u in starts for v in ends)}
 ia=[v for v in range(n) if not Aout>>v&1];oa=[v for v in range(n) if Aout>>v&1];im=[v for v in range(n) if not Mout>>v&1];om=[v for v in range(n) if Mout>>v&1]
 al=add(one,paths(range(n),ia));be=add(one,paths(oa,range(n)));a=add(one,paths(range(n),im));b=add(one,paths(om,range(n)));sg=paths(om,im);mu=paths(oa,im);ka=paths(om,ia);ta=paths(oa,ia);x=varx();y=vary();h=add(x,y,mul(x,y));U=add(mul(mu,h),x,mul(x,y));V=add(mul(ka,h),y,mul(x,y));sh=mul(sg,h);powers=[one]
 for r in range(1,n+1):powers.append(mul(powers[-1],sh))
 ds=[None]+[add(*(sc(comb(r+k-1,k-1),p) for r,p in enumerate(powers))) for k in range(1,4)]
 p1=mul(x,y,add(sc(2,mul(al,be)),mul(al,b),mul(a,be)),ds[1]);p2=mul(x,y,add(mul(al,b,U),mul(a,be,V),mul(a,b,add(mul(add(ka,mu,sc(2,ta)),h),sc(2,x),sc(2,y),mul(x,y)))),ds[2]);p3=mul(x,y,sc(2,mul(a,b,U,V)),ds[3]);ps=[p1,p2,p3]
 if any(c<0 for p in ps for c in p.values()):raise RuntimeError('nonnegative functional terms')
 value=lambda p:sum(c*bx[i]*by[j] for (M,i,j),c in p.items() if M==N-1)
 lower=2*bx[1]*by[1]*mul(al,be).get((N-1,0,0),0)
 vals=list(map(value,ps))
 if sum(vals)<lower or lower<=0:raise RuntimeError('strict lower bound')
 return {'nonnegative_terms':vals,'coefficient':sum(vals),'lower_bound':lower,'full_support_monomials':[[i,j,c] for (M,i,j),c in sorted(add(*ps).items()) if M==N-1]}
