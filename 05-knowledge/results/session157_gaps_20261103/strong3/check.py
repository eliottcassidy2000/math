from bank import automaton,rows,strong,dp

def need(x,what):
 if not x:raise AssertionError(what)

def step(v,x):
 out=[0]*54
 for o,order in enumerate(((0,1),(1,0))):
  for i in range(27):
   z=v[o*27+i]
   if not z: continue
   s=[i//9,i//3%3,i%3]
   for j in range(3):
    left=order[j-1] if j else None; right=order[j] if j<2 else None
    if s[j]==0 and left is not None and (x>>left&1):continue
    t=s[:];t[j]=2 if right is None or (x>>right&1) else 1
    out[o*27+t[0]*9+t[1]*3+t[2]]+=z
 return out

def P(v):return step(step(v,0),1)
def accept(v):
 v=step(step(v,2),1);out=0
 for o,order in enumerate(((0,1),(1,0))):
  # 1->0, so internal empty accepted only if order == (1,0)
  for i in range(27):
   s=[i//9,i//3%3,i%3]
   if 1 in s:continue
   if s[1]==0 and order!=(1,0):continue
   out+=v[o*27+i]
 return out

def lincomb(vs,c):return [sum(c[j]*vs[j][i] for j in range(len(c))) for i in range(54)]
v=[0]*54;v[0]=v[27]=1
v=P(v)  # family starts n=1
vs=[v]
for _ in range(5):vs.append(P(vs[-1]))
# q(X)=(X-1)^2(X-4)^2(X-9), ascending coefficients.
w=lincomb(vs,[-144,376,-337,123,-19,1])
for k in range(54):
 need(accept(w)==0,('observed annihilator',k,accept(w)))
 w=P(w)
# Cayley-Hamilton says the 54 consecutive zeros prove every later zero.
def g(n):
 z=25+150*n+(668+60*n)*4**n+432*9**n
 need(z%225==0,('integrality',n))
 return z//225
for n in range(1,6):need(automaton((0,1)*n+(2,1),0)==g(n),('initial',n))
for n in range(1,4):
 A=rows((0,1)*n+(2,1),0)
 need(strong(A),('strong',n));need(dp(A)==g(n),('direct DP',n,dp(A),g(n)))
print('Exact recurrence: 54-state Cayley-Hamilton observability certificate passed')
print('Initial values', [g(n) for n in range(1,6)], '; independent graph subset DP n=1..3')
# Modular check against the exact formula, selecting n beyond e so 9^n is negligible mod 3^e.
for e in range(1,8):
 M=3**e
 def gm(n):return ((25+150*n+(668+60*n)*pow(4,n,225*M)+432*pow(9,n,225*M))%(225*M))//225
 vals={gm(n) for n in range(e+2,e+2+M)}
 need(vals==set(range(M)),('3-adic',e,len(vals)))
print('All residues mod 3^e covered for e=1..7')
M=21;D=225*M;period=630
need(pow(4,period,D)==1, '4-period')
need(pow(9,2+period,D)==pow(9,2,D), '9-eventual-period')
need((150*period)%D==0 and (60*period)%D==0,'linear-period')
def g21(n):return ((25+150*n+(668+60*n)*pow(4,n,D)+432*pow(9,n,D))%D)//225
image={g21(n) for n in range(1,2+period)}
need(set(range(M))-image=={2,8,11,14,17,20},('mod21 image',image))
print('Hostile control: eventual period 630 and mod 21 image confirmed')
M=75;D=225*M;period=2250
need(pow(4,period,D)==1,'mod75 4-period')
need(pow(9,2+period,D)==pow(9,2,D),'mod75 9-period')
need((150*period)%D==0 and (60*period)%D==0,'mod75 linear period')
def g75(n):return ((25+150*n+(668+60*n)*pow(4,n,D)+432*pow(9,n,D))%D)//225
image={g75(n) for n in range(1,period+2)}
missing_div3=[a for a in range(0,M,3) if a not in image]
need(missing_div3==[9,12,24,27,54,57,69,72],('mod75 missing',missing_div3))
print('Union obstruction: modulo 75 these divisible-by-3 residues are missed by both families:',missing_div3)
