from math import isqrt

def check(condition, description):
    if not condition:
        raise AssertionError(description)

for d in range(1, 13):
    q = 1 if d == 1 else 1 + 2**(d-1)
    t = 2**(d+1) + 2 * (3**(d-1) - 1)
    check(q % 7 != 0, ('q7', d))
    check((t % 7 == 0) == (d % 6 in (4,5)), ('t7', d))
    if d % 2 == 0:
        check(q % 3 == t % 3 == 0, ('mod3', d))

def restricted_neighbors(v,n):
    ans=[]
    for s in range(2,isqrt(v+n)+1):
        w=s*s-v
        if 1 <= w <= n and w != v and abs(w-v)%6 not in (4,5):
            ans.append(w)
    return ans
check(restricted_neighbors(10,53)==[], '10 isolated')
check(restricted_neighbors(52,143)==[], '52 isolated')
for v,w in ((52,144),(55,141),(58,138)):
    check(restricted_neighbors(v,166)==[w], ('degree one',v))

# X' has reversed edges 3→1 and 4→1; X has 3→0.
# p maps X' labels 1..4 to X labels 0..3.
p={1:3,2:0,3:1,4:2}
def edge(i,j,back):
    pair=(max(i,j),min(i,j))
    return i>j if pair in back else i<j
for i in range(1,5):
    for j in range(1,5):
        if i != j:
            check(edge(i,j,{(3,1),(4,1)})==edge(p[i],p[j],{(3,0)}), ('isomorphism',i,j))

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
for g,pv in ((1,19),(7,1171)):
    check(prime(pv) and 5+2**(g+6)==7*pv, ('prime instance',g,pv))
print('PASS mod-3 and mod-7 factors, degree certificate, isomorphism, prime cofactors')
