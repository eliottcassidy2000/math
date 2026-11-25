def g(u):
    return sum((2 if i-j==1 else 2**(i-j-1)) for i in range(u.bit_length())
               for j in range(i) if u>>i&1 and not u>>j&1)

def m(q,i,j):
    u,v=2*q+i,2*q+j
    return 2*u*v+g(u)+g(v)+2*q+2

def ranges(p):
    x=1<<(p-2)
    E=8*x*x-2*x
    got={(t,i,j):m(x+t,i,j)-E for t in (-1,0)
         for i in (-1,0) for j in (1,2)}
    assert got=={(-1,-1,1):-23*x//2+6,(-1,-1,2):-11*x//2,
                 (-1,0,1):-7*x+4,(-1,0,2):-x,
                 (0,-1,1):5*x,(0,-1,2):19*x//2,
                 (0,0,1):11*x+2,(0,0,2):31*x//2+4}
    assert 16*x-8>11*x+2 and 20*x+10>11*x+2
    assert min(map(abs,got.values()))==x
    if p%2:
        assert E%3==1
        assert [k for k,v in got.items() if v%3==0]==[(0,0,1)]
    return got

if __name__=='__main__':
    for p in range(5,62):ranges(p)
    print('exponent-61 offsets:',ranges(61))
