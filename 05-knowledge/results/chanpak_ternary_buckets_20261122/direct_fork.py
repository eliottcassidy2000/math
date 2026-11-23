"""A direct physical tripling arc with a persistent upper neighbor and both unit neighbors on a dyadic family."""
from itertools import combinations

from physical_navigation import address_height, direct, flipped, star_rows, tagged_data
from e15_pair import join
from fermat_receiver import p as read


PARENT=(0,1,11,19,7,31,63)
CENTER=flipped(PARENT,4,6)
LOWER=flipped(CENTER,2,3)
UPPER=flipped(CENTER,4,5)
TAG=11
EXPECTED=[
    (3,192,21,33,3,3,5),
    (9,336,63,54,9,8,15),
    (5,256,63,50,5,4,17),
    (9,336,63,60,9,8,17),
]


def checks():
    data=[tagged_data(rows,TAG) for rows in (PARENT,CENTER,LOWER,UPPER)]
    assert data==EXPECTED
    for x in (0,1,2,5):
        for y in range(128):
            p,c,l,u=[address_height(d,x,y) for d in data]
            assert c-3*p == -240*x*y-45*x-(x%2)
            assert u-c == 6*x+2
            if x==0:
                assert l-c == 2-4*(y%2)-4*read(y)
                assert (l==c-2) == (y>0 and y&(y+1)==0)
    for y in (0,1,3,5):
        hosts=[star_rows(r,TAG,0,y) for r in (PARENT,CENTER,LOWER,UPPER)]
        assert flipped(hosts[0],4,6)==hosts[1]
        assert flipped(hosts[1],2,3)==hosts[2] and flipped(hosts[1],4,5)==hosts[3]
        assert [direct(t) for t in hosts]==[address_height(d,0,y) for d in data]
    for m in range(8):
        b=tagged_data((2,4,1),m)
        jp,jc,jl,ju=[join(d,b) for d in data]
        for y in range(17):
            p,c,u=[address_height(d,0,y) for d in (jp,jc,ju)]
            assert c==3*p and u-c==2*b[0] and p%3==2*b[0]%3
    for a_marks in range(8):
        a=tagged_data((2,4,1),a_marks)
        for b_marks in range(8):
            b=tagged_data((2,4,1),b_marks)
            p,c=[join(join(a,d),b) for d in data[:2]]
            for y in (0,1,2,7,8):
                delta=address_height(c,0,y)-3*address_height(p,0,y)
                assert delta==-240*a[2]*(b[1]*y+b[3])-b[0]*(45*a[2]+a[4])
                assert (delta==0)==(a_marks==0)
    options=[]
    for i,j in combinations(range(8),2):
        d=tagged_data(flipped(CENTER,i,j),TAG) if j<7 else tagged_data(CENTER,TAG^(1<<i))
        if (d[0],d[2]) == (data[1][0],data[1][2]):
            options.append(((i,j),d[4]-data[1][4],d[6]-data[1][6]))
    assert options==[((0,1),0,0),((4,5),0,2)]
    print("Direct fork profiles",data)
    print("At the seven-vertex center, h,c-preserving core and tag edits:", options)
    # A distinct persistent lower-side module, with every right digit free.
    rows=(0,5,17,7,75,31,47)
    triple=[tagged_data(r,5) for r in (flipped(rows,4,6),rows,flipped(rows,3,4))]
    for y in range(128):
        p,c,l=[address_height(d,0,y) for d in triple]
        assert p==5*read(y)+15*y+5*(y%2)+9 and c==3*p==l+2
    print("Persistent lower-side fork profiles",triple)


if __name__=="__main__":
    checks()
