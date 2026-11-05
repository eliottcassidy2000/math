import json,math,hashlib

def need(b,s):
 if not b:raise RuntimeError(s)
def g(L):
 need(L>=1,'positive L');v=1
 while L>1:
  if L%2==0:L//=2
  else:L=(L+1)//2;v+=L-L%2
 return v
def f(n):return g(n+1)
def fromword(w,n=0,v=1):
 for b in w:
  if b=='0':v+=n+n%2;n*=2
  else:n=2*n+1
 return n,v
def phi(n):return sum(math.gcd(j,n)==1 for j in range(1,n+1))
def general(M,u,z,extra=1):
 need(M>=2 and M%2==0 and z%2==1,'valid joint residues')
 q=M;k=0
 while q%2==0:q//=2;k+=1
 a=max(2,2*phi(q));D=((1<<a)-1)*((1<<k)-1)+1
 need(math.gcd(D,M)==2,'loop gcd')
 r=(u-1)%M;mod=M//2
 j=0 if mod==1 else (((z-f(r))//2)*pow(D//2,-1,mod))%mod
 j+=extra*mod
 w=('1'*a+'0'*k)*j+format(r,f'0{(M-1).bit_length()}b')
 return w,{'a':a,'k':k,'j':j,'Delta':D}
def dyadic(k,u,z,extra=1):
 M=1<<k;r=(u-1)%M;j=((f(r)-z)//2)%(M//2)
 w=('1'+'0'*k)*extra
 for i in range(k-1):
  if j>>i&1:
   if i: w+='1'*(i+1)+'0'*k
   w+='11'+'0'*k
 w+=format(r,f'0{k}b')
 return w
def graph(L):
 n=L-1;s=n.bit_length();rows=[0]*(s+1)
 for i in range(s):
  for j in range(i+1,s):rows[i]|=1<<j
  if n>>i&1:rows[i]|=1<<s
  else:rows[s]|=1<<i
 return rows
def seed_graph(L):
 G=graph(L);s=len(G);rows=[(((1<<s)-1)^(1<<(s-1)))<<1]
 rows += [r<<1 for r in G];rows[-1]|=1
 return rows
def hp(rows):
 n=len(rows);D=[{} for _ in range(1<<n)]
 for v in range(n):D[1<<v][v]=1
 for S in range(1,1<<n):
  for v,c in D[S].items():
   out=rows[v]&~S
   while out:
    b=out&-out;out-=b;w=b.bit_length()-1
    T=D[S|b];T[w]=T.get(w,0)+c
 return sum(D[-1].values())
def tournament(rows):
 return all(not(rows[i]>>i&1) for i in range(len(rows))) and all(((rows[i]>>j)&1)+((rows[j]>>i)&1)==1 for i in range(len(rows)) for j in range(i))
def main():
 general_checks=dyadic_checks=loop_checks=prefix_checks=physical=0;digest=hashlib.sha256();examples=[]
 for M in range(2,42,2):
  for u in range(M):
   for z in range(1,M,2):
    w,data=general(M,u,z);n,v=fromword(w)
    need(n==int(w,2) and v==g(n+1),'independent left recurrence')
    need((n+1)%M==u and v%M==z,'general joint residues')
    need(len(w)<=(M-1)*(M+data['k'])+(M-1).bit_length(),'general length bound')
    digest.update(f'{M},{u},{z},{n+1},{v}\n'.encode());general_checks+=1
  _,d=general(M,1,1);loop='1'*d['a']+'0'*d['k']
  for t in range(4):
   for value in (1,3,17):
    nn,vv=fromword(loop,M*t,value)
    need(nn%M==0 and (vv-value)%M==d['Delta']%M,'loop from arbitrary state');loop_checks+=1
 for k in range(1,11):
  M=1<<k
  pairs=[(u,z) for u in range(M) for z in range(1,M,2)] if k<=5 else [((73*i+11)%M,(82*i+1)%M) for i in range(128)]
  for u,z in pairs:
   w=dyadic(k,u,z);n,v=fromword(w)
   need(v==g(n+1) and (n+1)%M==u and v%M==z,'short joint residues')
   need(len(w)<=4*k*k+2*k+1,'quadratic word bound')
   w2=dyadic(k,u,z,2);n2,v2=fromword(w2)
   need(n2>n and (n2+1)%M==u and v2%M==z,'arbitrarily large neutral padding')
   dyadic_checks+=1
   if u==1 and z in (1,3,7) and k<=4:examples.append({'k':k,'target_residue':z,'L':n+1,'g':v,'H':v+n,'word_bits':len(w)})
 for k in (16,32,64):
  M=1<<k
  for u,z in ((1,7),(M-3,M-1),(M//2+1,12345)):
   w=dyadic(k,u,z);n,v=fromword(w)
   need(v==g(n+1) and (n+1)%M==u%M and v%M==z%M,'large exact joint residues')
   need(len(w)<=4*k*k+2*k+1,'large quadratic bound');dyadic_checks+=1
 for k in range(1,7):
  M=1<<k
  for h,e,R in ((1,1,1),(3,1,3),(5,1,7),(13,2,18)):
   for N in range(1,M,2):
    u=(3*N)%M;z=(pow(h,-1,M)*(N-R*u+e*(u%2)))%M
    w=dyadic(k,u,z);n,v=fromword(w);L=n+1
    H=h*v+R*L-e*(L%2)
    need(H%M==N,'fixed-prefix target residue');prefix_checks+=1
 # Literal small graphs; these checks do not use the scalar seam recurrence.
 literal_L=set(range(2,65))
 for k in (1,2):
  for u in range(1<<k):
   for z in range(1,1<<k,2):literal_L.add(int(dyadic(k,u,z),2)+1)
 for L in sorted(literal_L):
  G=graph(L);T=seed_graph(L)
  need(tournament(G) and tournament(T),'literal adjacency tournament')
  need(hp(G)==g(L),'literal gadget H')
  need(hp(T)==g(L)+L-L%2,'literal singleton-prefix H')
  cycle=list(range(len(T)-1))+[len(T)-1]
  need(all(T[cycle[i]]>>cycle[(i+1)%len(T)]&1 for i in range(len(T))),'explicit Hamilton cycle')
  need(all(T[i]>>j&1 for i in range(len(T)-1) for j in range(i+1,len(T)-1)),'root deletion transitive');physical+=1
 print(json.dumps({'status':'PASS','general_joint_residue_checks':general_checks,'dyadic_joint_residue_checks':dyadic_checks,'arbitrary_loop_state_checks':loop_checks,'fixed_prefix_residue_checks':prefix_checks,'literal_graph_pairs':physical,'general_constructed_pair_sha256':digest.hexdigest(),'examples':examples,'scope':'Residue universality only; no exact-target, cofinite or infinite-primality assertion.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
