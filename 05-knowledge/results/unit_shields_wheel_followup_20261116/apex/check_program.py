from itertools import product
from check_shield import gap, fast, literal

def common(t,k):
    return tuple(3*int(b) for b in f'{t:0{k}b}')

def decode(d,k,u):
    if d%4==0:
        t=(2**k-d)//2
        t=min(t,2**k-1)
        c=(1,2)
    else:
        t=(2**(k+1)-2-d)//4
        c=(0,1,2,3)
    assert 0<=t<2**k
    s=common(t,k)
    p=s+c+s[::-1]
    return u+p+tuple(3-z for z in u[::-1])

def test():
    cases=0
    for k in range(2,8):
        for d in range(-2**k,2**k+1,2):
            for tail in product((0,3), repeat=4):
                u=(3,)+tail
                w=decode(d,k,u)
                a,b=fast(w)
                assert a-b==d and gap(w)==d
                cases+=1
            if k==2:
                w=decode(d,k,(3,))
                assert literal(w,True)-literal(w,False)==d
    print(cases,'exact three-run controls of all signed even residues, plus literal DP at k=2.')

if __name__=='__main__':test()
