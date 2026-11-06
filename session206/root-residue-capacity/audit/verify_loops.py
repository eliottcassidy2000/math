"""Independent exact checks of the chain's residue-adjusting loops."""
import json

def Q(n):
    q=0
    while n:
        a=n//2
        if n%2:q+=a+a%2
        n=a
    return q

checks=0
for p in range(3,100,2):
    b=1
    while pow(2,b,p)!=1:b+=1
    r=2*((1<<b)-1)
    if r%(2*p):raise RuntimeError(('block residue',p))
    for j in range(20):
        v=2*p*j
        w=(v<<(b+1))+r
        if (Q(w)-Q(v))%p!=(p-1):raise RuntimeError(('Q increment',p,j))
        if w%(2*p):raise RuntimeError(('first coordinate',p,j))
        checks+=1
print(json.dumps({'status':'PASS','odd_moduli':49,'loop_checks':checks,'largest_modulus':99},indent=2))
