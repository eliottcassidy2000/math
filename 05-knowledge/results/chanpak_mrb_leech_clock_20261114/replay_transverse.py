from residue_round import eval_count

def kernel(q):
    return (eval_count(q,(1,0),(1,0)),eval_count(q,(1,0),(0,1)),eval_count(q,(0,1),(0,1)))

def state(gaps):
    a,b=1,1<<gaps[0]
    q=0
    for r in gaps[1:]:
        if r==0:
            q+=1
        else:
            e,f,g=kernel(q)
            a,b=e*a+f*b,(f*a+g*b)<<r
            q=0
    assert q==0
    return a,b

def packing(gaps):
    P=[1]
    cum=[0]
    for r in gaps: cum.append(cum[-1]+r)
    for j,x in enumerate(cum):
        P.append(P[-1]+sum((1<<(x-cum[i]))*P[i] if x>cum[i]
                            else 2*P[i]*int((j-i)%2==0) for i in range(j)))
    return P[-1]

words=((1,0,0,3,0,0,1),(2,0,0,0,0,0,0,2))
z=tuple(map(state,words))
assert z==((291,758),(179,852))
assert [47*a+56*b for a,b in z]==[56125]*2
rows=[]
for q in range(18):
    e,f,g=kernel(q)
    rows.append((56*e-47*f,56*f-47*g))
assert [q for q,(a,b) in enumerate(rows) if b==0]==[10,13,14]
assert all(rows[q][1]>0 for q in (15,16,17))
expect={10:((423939,505125),(423941,505125)),14:((4851811,5780875),(4851809,5780875))}
for q in (10,13,14):
    for r in (1,2,3,6,17):
        H=[packing(w+(0,)*q+(r,)) for w in words]
        assert H==[eval_count(q,(1,1<<r),v) for v in z]
        if q in expect:
            assert H==[a+(b<<r) for a,b in expect[q]]
        else:
            assert H[0]==H[1]
print('two physical transverse states; two exact width-independent +/-2 repairs')
