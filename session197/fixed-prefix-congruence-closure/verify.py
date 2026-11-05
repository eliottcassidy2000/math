import math,json,collections

def need(b,s):
 if not b:raise RuntimeError(s)
def g(L):
 v=1
 while L>1:
  if L%2==0:L//=2
  else:L=(L+1)//2;v+=L-L%2
 return v
def constants(L):
 if L==1:return 1,1,1
 h,Q,B=constants((L+1)//2)
 return (h+((L+1)//2-((L+1)//2)%2) if L%2 else h,Q+B if L%2 else Q,2*B)
def bank(B):
 S={(1,1,1)};q=collections.deque(S)
 while q:
  h,e,R=q.popleft()
  for L in range(1,B//R+1):
   a,Q,D=constants(L);x=(a*h-(L%2)*e+L*R,h-(L==1)*e,Q*h+D*R)
   if x[2]<=B and x not in S:S.add(x);q.append(x)
 return sorted(S)
def egcd(a,b):
 if not b:return a,1,0
 d,x,y=egcd(b,a%b);return d,y,x-(a//b)*y
def joint(M,u,z):
 q=M;k=0
 while q%2==0:q//=2;k+=1
 a=2
 while pow(2,a,q)!=1%q:a+=1
 D=((1<<a)-1)*((1<<k)-1)+1;r=(u-1)%M;mod=M//2
 j=0 if mod==1 else ((z-g(r+1))//2*pow(D//2,-1,mod))%mod
 j+=mod;word=('1'*a+'0'*k)*j+format(r,f'0{(M-1).bit_length()}b')
 return int(word,2)+1
def construct(P,m,r):
 h,e,R=P;d,x,y=egcd(h,R);D=math.gcd(d,m)
 N=r if r%2 else r+m
 need(N%2==1,'odd representative')
 t=next(t for t in (0,1) if (N+t*e)%D==0);K=h+t*(R-e);b=(N-K)//2;mod=m//D
 s=0 if mod==1 else (b//D*pow(d//D,-1,mod))%mod
 w,v=x*s,y*s;u=(2*v+t)%(2*m);z=(2*w+1)%(2*m)
 L=joint(2*m,u,z);H=h*g(L)+R*L-e*(L%2)
 need(H%m==r and L%(2*m)==u and g(L)%(2*m)==z,'actual local witness')
 return L,H

def main():
 P=bank(100);checks=cardinal=actual=0;largest=0
 for h,e,R in P:
  need(h%2==1 and (R-e)%2==0,'physical parity');d=math.gcd(h,R)
  for m in range(1,33):
   D=math.gcd(d,m);image={(h+2*h*w+2*R*v+t*(R-e))%m for w in range(m) for v in range(m) for t in (0,1)}
   expected={r for r in range(m) if (m%2 or r%2) and (r%D==0 or (r+e)%D==0)}
   need(image==expected,'complete local cosets');checks+=1
   size=(m//D if m%2 else m//(2*D))*(1 if e%D==0 else 2)
   need(len(image)==size,'exact residue cardinality');cardinal+=1
 for p in ((1,1,1),(3,1,3),(9,3,9),(33,1,33),(13,2,18)):
  h,e,R=p;d=math.gcd(h,R)
  for m in range(1,25):
   D=math.gcd(d,m)
   for r in range(m):
    if (m%2 or r%2) and (r%D==0 or (r+e)%D==0):
     L,H=construct(p,m,r);actual+=1;largest=max(largest,L.bit_length())
 need((33,1,33) in P and (3,1,3) in P,'named prefixes physically generated')
 N=143;need(N%2==1 and all(N%p in (0,p-1) for p in (3,11)) and N%33 not in (0,32),'primewise branch hostile')
 print(json.dumps({'status':'PASS','physical_prefix_states':len(P),'complete_modular_images':checks,'cardinality_checks':cardinal,'actual_constructed_local_witnesses':actual,'maximum_constructed_index_bit_length':largest,'composite_branch_hostile':{'prefix':[33,1,33],'N':143,'mod3':143%3,'mod11':143%11,'mod33':143%33},'scope':'Complete congruence criterion, not exact-target sufficiency.'},sort_keys=True,indent=2))
if __name__=='__main__':main()
