from replay import DATA, TT, add, profile, strong, substitute

EXPECTED = [
    [(0,1,False),(3,7,True)],
    [(0,1,False),(4,5,True),(32,3,False),(34,7,True),(38,9,True)]
]
for ix,(R,_,a,b) in enumerate(DATA):
    J1,J2=2*a+b,4*a+b
    RR=substitute(R,TT(2))
    actual=[]
    for m in range(1<<(len(R)-1)):
        T=add(R,m)
        L1,L2=profile(T)[0], profile(add(RR,m))[0]
        if L1%J1==0 and L1*J2==L2*J1:
            actual.append((m,L1//J1,strong(T)))
    if actual!=EXPECTED[ix]:
        raise ValueError(actual)
    print("exact multiplier menu",ix+1,actual)
