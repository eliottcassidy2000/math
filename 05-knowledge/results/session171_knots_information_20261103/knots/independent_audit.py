import math,json,itertools

def need(x,m):
    if not x:raise RuntimeError(m)
# A mixed-prime finite nilpotent target D8 x C3, independent of octonion permutations.
G=list(itertools.product(range(4),range(2),range(3)))
e=(0,0,0)
def mul(x,y):return ((x[0]+(-1)**x[1]*y[0])%4,(x[1]+y[1])%2,(x[2]+y[2])%3)
def power(x,n):
    a=e
    for _ in range(n):a=mul(a,x)
    return a
count=0
for p in range(1,9):
    for q in range(1,9):
        if math.gcd(p,q)!=1:continue
        actual={(a,b) for a in G for b in G if power(a,p)==power(b,q)}
        expected={(power(g,q),power(g,p)) for g in G}
        need(actual==expected and len(actual)==24,'mixed-prime nilpotent representation')
        count+=1
# General triangular integer lattices, not just isotropic mZ^2.
lattice_checks=0
for a in range(1,6):
    for c in range(1,6):
        for b in range(a):
            D=a*c
            image={((c*x-b*y)%D,(a*y)%D) for x in range(D) for y in range(D)}
            need(len(image)==D,'quotient order')
            for p in range(-4,5):
                for q in range(-4,5):
                    if math.gcd(p,q)!=1:continue
                    step=((c*p-b*q)%D,a*q%D)
                    v=(0,0);orbit=[]
                    while v not in orbit:
                        orbit.append(v);v=((v[0]+step[0])%D,(v[1]+step[1])%D)
                    d=len(orbit)
                    need(D%d==0,'orbit divides degree')
                    u1=d*(c*p-b*q)//D;u2=d*a*q//D
                    need(d*(c*p-b*q)%D==0 and d*a*q%D==0,'lift in lattice')
                    need(math.gcd(u1,u2)==1,'lift primitive in lattice basis')
                    left=set(image);components=0
                    while left:
                        z=next(iter(left));coset={((z[0]+x)%D,(z[1]+y)%D) for x,y in orbit}
                        need(coset<=left,'orbit partition');left-=coset;components+=1
                    need(components==D//d,'component count')
                    lattice_checks+=1
extrema=[]
for k in range(1,21):
    values=[]
    for a in range(k+1):
        for b in range(k-a+1):
            for c in range(k-a-b+1):
                d=k-a-b-c;values.append(3**a+3**b+3**c+3**d)
    u,r=divmod(k,4)
    need(min(values)==(4-r)*3**u+r*3**(u+1),'minimum occupancy')
    need(max(values)==3**k+3,'maximum occupancy')
    extrema.append({'k':k,'min':min(values),'max':max(values)})
# Explicit isolated defect is an allowed cyclic trefoil representation.
# An involution g gives A=g^3=g, B=g^2=1, hence A^2=B^3=1.
g=(0,1,0)
need(g!=e and power(g,2)==e,'nontrivial involution')
A=power(g,3);B=power(g,2)
need(power(A,2)==power(B,3),'isolated defect can be meridian image')
print(json.dumps({'status':'PASS','mixed_prime_nilpotent_cases':count,'nonisotropic_lattice_slope_cases':lattice_checks,'occupancy_extrema':extrema,'isolated_involution_cyclic_representation':True},indent=2))
