from fermat_profiles import rooted,canon,extend,profile,strong
T={(0,)}
letters={(1,2)}
C=set()
for n in range(1,7):
    if n>1:T={canon(extend(A,c)) for A in T for c in range(1<<(n-1))}
    for A in T:
        if strong(A):
            p,q=profile(A)
            if q<=33:letters.add((p,q))
            C.update((h,l) for h,l,_,_,_,_ in rooted(A))
D={(1,1)}
for e in range(1,34):
    D.update((h*a,q*b) for h,q in tuple(D) if q==e for a,b in letters if q*b<=33)
if {x for x in D if x[0]>1}!={(3,12),(3,24),(5,32)}:
    raise RuntimeError('low tail catalogue')
if any(h*a+(l*b<<j)==67 for h,l in C for a,b in D for j in range(1,7)):
    raise RuntimeError('67 falsely attained by seam')
S={h+(l<<j) for h,l in C for j in (1,2,4)}
if set(range(1,69,2))-S!={1,7,11,13,21,47,67}:
    raise RuntimeError('wrong low-window exclusions')
print('67 excluded for every one-way near-transitive seam, even with arbitrary source order')
