from itertools import combinations

def schur_count(A):
    S=set(A)
    return sum(1 for a in A for b in A if a+b in S)  # ordered pairs (a,b), c=a+b

# Claim: among n-element sets of positive integers, interval {1..n} maximizes T(A)=#{(a,b,c)in A^3:a+b=c},
# with max = n(n-1)/2, equality iff scaled interval {d,2d,...,nd}.
print("n | T({1..n}) | n(n-1)/2 | brute-max over subsets of [1..N] | maximizers(sample)")
for n in range(2,7):
    N = 3*n+2
    interval = list(range(1,n+1))
    Ti = schur_count(interval)
    best=-1; maximizers=[]
    for A in combinations(range(1,N+1), n):
        t=schur_count(A)
        if t>best:
            best=t; maximizers=[A]
        elif t==best:
            maximizers.append(A)
    # which maximizers are scaled intervals {d,2d,...,nd}?
    def is_scaled_interval(A):
        A=sorted(A); d=A[0]
        return all(A[i]==(i+1)*d for i in range(len(A)))
    si = [A for A in maximizers if is_scaled_interval(A)]
    nonsi = [A for A in maximizers if not is_scaled_interval(A)]
    print(f"{n} | {Ti} | {n*(n-1)//2} | {best} | #max={len(maximizers)}, #scaled-int={len(si)}, #other={len(nonsi)}")
    if nonsi:
        print(f"    NON-scaled-interval maximizers exist: {nonsi[:3]}")

# scaling invariance + a few non-interval APs
print("\nScaling / AP-offset test:")
for A in [(1,2,3),(2,4,6),(3,6,9),(1,3,5),(2,3,4),(1,2,4)]:
    print(f"  T({A})={schur_count(list(A))}  scaled-int={sorted(A)[0]* (len(A)) == max(A) and all(sorted(A)[i]==(i+1)*sorted(A)[0] for i in range(len(A)))}")

# upper bound check T(A) <= n(n-1)/2 over MANY random positive-real sets
import random
random.seed(1)
viol=0; tight=0
for _ in range(200000):
    n=random.randint(2,6)
    A=sorted(random.sample(range(1,40), n))
    t=schur_count(A)
    if t> n*(n-1)//2: viol+=1
    if t==n*(n-1)//2: tight+=1
print(f"\nRandom integer sets: violations of T<=n(n-1)/2: {viol}; tight(=bound): {tight}")
