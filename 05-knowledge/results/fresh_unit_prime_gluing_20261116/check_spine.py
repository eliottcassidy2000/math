from itertools import product
def scan(word):
    # Counts for each fixed apex order by 3 increasing runs.
    def step(v,a,b):
        out=[0]*12
        for f,l,m in product(range(2),range(2),range(3)):
            x=v[6*f+3*l+m]
            out[6*a+3*l+m]+=x
            if a or m:out[6*f+3*l+1+b]+=x
            if b or l:out[6*f+3+m]+=x
        return out
    A=B=[1]+[0]*11
    for a,b in word:
        A=step(A,a,b);B=step(B,b,a)
    return A[0]+A[3]+A[1]+A[4]+B[1]+B[4]
def one(w):
    # Single apex: positive odd-cycle packs have at most one cycle.
    s=0
    for j in range(len(w)):
        if w[j]:continue
        for i in range(j):
            if w[i]:s+=1 if j==i+1 else 2**(j-i-2)
    return 1+2*s
checks=0
for n in range(1,8):
    for z in product(range(4),repeat=n):
        w=[divmod(t,2) for t in z]
        p=one([a for a,b in w]);q=one([b for a,b in w])
        if p>1 and q>1 and (p>3 or q>3):
            assert scan(w)>p+q-1
            checks+=1
assert scan([(1,1),(0,0)])==5
print('Transitive shared-spine obstruction checked on',checks,'qualifying tapes; seed equality 3+3−1=5 passes.')
