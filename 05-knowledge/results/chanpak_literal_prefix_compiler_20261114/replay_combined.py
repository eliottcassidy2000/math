from bridge import h_f2,bridge,source,rooted_profile
from physical_prefix import prime
from physical_seed import step
from construct_target import actual,strong

C=(8,1,3,6)
assert h_f2(C)==(5,32)
P=bridge(source(C),0,1)
assert P==[62,112,98,102,108,64,1]
assert rooted_profile(P,len(P)-1)==(69,69,-5)
values=[]
for order in ((1,2),(2,1)):
    T=P.copy()
    z=(5,64)
    for q,r in ((14,1),(order[0],1),(order[1],1),(7,2)):
        for j in range(q):T=bridge(T,len(T)-1,0)
        T=bridge(T,len(T)-1,r)
        z=step(q,r,z)
    assert len(T)==40 and actual(T) and strong(T)
    assert sum(z)%7==0
    values.append(sum(z)//7)
assert values==[3539307127,3539307117]
assert prime(values[0]) and values[1]%3==0
print('literal seven-vertex nonterrace source, two 40-vertex switch witnesses, and prime cofactor checked')
