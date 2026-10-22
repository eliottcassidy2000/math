from math import comb,gcd
import json
N=60
def mul(a,b):return [sum(a[j]*b[i-j] for j in range(i+1)) for i in range(N+1)]
def power(a,k):
 r=[1]+[0]*N
 for _ in range(k):r=mul(r,a)
 return r
def es(k,c):return [1]+[c*sum(d**(k-1) for d in range(1,n+1) if n%d==0) for n in range(1,N+1)]
e4,e6=es(4,240),es(6,-504)
e43,e62=power(e4,3),power(e6,2)
delta=[(a-b)//1728 for a,b in zip(e43,e62)]
assert all((a-b)%1728==0 for a,b in zip(e43,e62))
p=[1]+[0]*N
for n in range(1,N+1):
 f=[0]*(N+1)
 for j in range(min(24,N//n)+1):f[n*j]=(-1)**j*comb(24,j)
 p=mul(p,f)
assert delta==[0]+p[:-1]
assert power(e4,2)==es(8,480)
assert mul(e4,e6)==es(10,-264)
assert mul(power(e4,2),e6)==es(14,-24)
inv=[1]+[0]*N
for n in range(1,N+1):inv[n]=-sum(p[k]*inv[n-k] for k in range(1,n+1))
jshift=mul(e43,inv) # q*j(q)
leech=[a-720*b for a,b in zip(e43,delta)]
golayA=[a-672*b for a,b in zip(e43,delta)]
d24=[a+384*b for a,b in zip(e43,delta)]
assert leech[:4]==[1,0,196560,16773120]
assert mul(leech,inv)==[a-(720 if i==1 else 0) for i,a in enumerate(jshift)]
# Ordered positive two-cube representations, search through 1729 inclusive.
cubes={}
for a in range(1,13):
 for b in range(a,13):
  n=a**3+b**3
  if n<=1729:cubes.setdefault(n,[]).append((a,b))
assert min(n for n,v in cubes.items() if len(v)>1)==1729
assert 7*13*19==1729 and all(1728%(p-1)==0 for p in [7,13,19])
assert all(pow(a,1728,1729)==1 for a in range(1,1729) if gcd(a,1729)==1)
# Root enumeration and half-coset score floor count.
roots8=2*8*7+2**7
assert roots8==240 and comb(8,4)==70
# Exact E8 Lie graded character q^(1/3)*Theta_E8/eta^8 first terms.
p8=[1]+[0]*N
for n in range(1,N+1):
 f=[0]*(N+1)
 for j in range(min(8,N//n)+1):f[n*j]=(-1)**j*comb(8,j)
 p8=mul(p8,f)
inv8=[1]+[0]*N
for n in range(1,N+1):inv8[n]=-sum(p8[k]*inv8[n-k] for k in range(1,n+1))
char8=mul(e4,inv8)
assert char8[:3]==[1,248,4124]
# Fixed H=3: one strong C3 with singleton sequences on both sides.
counts3=[0 if n<3 else n-2 for n in range(15)]
# Fixed H=63 inherited complete two order-8 strong atoms.
counts63=[0 if n<8 else 2*(n-7) for n in range(15)]
out={'precision':'identities through q^60 inclusive, inverse Delta recorded after q shift','Delta':delta[:10],'inverse_Delta_q_shifted':inv[:8],'j_q_shifted':jshift[:8],'Leech_theta':leech[:6],'Golay_Construction_A_theta':golayA[:6],'D24_plus_theta':d24[:6],'E8_character_q_one_third_shifted':char8[:6],'1729_cube_pairs':cubes[1729],'1729_Carmichael_units_checked':sum(gcd(a,1729)==1 for a in range(1,1729)),'E8_roots':roots8,'score_half_coset_roots':comb(8,4),'H3_counts':counts3,'H63_counts':counts63,'status':'ALL ASSERTIONS PASSED'}
print(json.dumps(out,indent=2))
