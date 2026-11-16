from itertools import product
def h(w):
    A,p,H=0,0,1
    for c in w:
        e=int(c)
        if not e:H+=A+p
        A,p=2*A+e,e
    return H,A
t1=[]
for n in range(2,10):
    for middle in product('01',repeat=n-2):
        x='1'+''.join(middle)+'1'
        H,A=h(x)
        t=(2*(A+1)-3*H-1)//2
        if t==1:t1.append(x)
        count,Q=h(x+'10'*t);Q+=1
        assert 2*Q-3*count==1
        assert ((85*Q-14)%120,(85*Q-12)%120)==((21,23) if t==1 else (81,83))
assert t1==['101']
assert {(520569*2**t+123357)%120 for t in range(3,7)}=={69,21,45,93}
assert {(867615*2**t+123357)%120 for t in range(3,7)}=={117}
print('Physical bridge residues and transient replayed.')
