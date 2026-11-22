"""All-cut obstruction and complete single-homogeneous-run certificates for the E15 pair."""
from e15_pair import A, D, V, EDGE, MARKS, flipped, join, tagged_data

LOWER_61=2**60*(2**61-1)-1


def tape(bit, length):
    q=2**length
    return (1,q,q-1,0,1,0,1) if bit else (1,q,0,q-1,0,1,1)


def fold(blocks):
    result=blocks[0]
    for block in blocks[1:]:
        result=join(result,block)
    return result


def main():
    a,v=tagged_data(A,MARKS[0]),tagged_data(V,MARKS[2])
    sides=[tagged_data(rows,MARKS[1]) for rows in (flipped(D,*EDGE),D)]
    sides=[(a,d,v) for d in sides]
    lower,upper=map(fold,sides)
    assert (a[0],a[4],v[0],v[5])==(5,2,627,376)
    assert lower==(90915,147456000,40180840,50982110,62060,72732,14077887)
    assert upper==lower[:-1]+(lower[-1]+2,)
    h,q,c,t,e,r,k=lower
    assert all(x%5==0 for x in (h,q,c,t,e)) and (r%5,k%5)==(2,2)
    assert LOWER_61%5 not in {k%5,(k+r)%5}
    # A complete check of the quoted half-gaps for insertions of homogeneous runs.
    result={}
    suffix=[tape(int(b),1) for b in '0100']
    for cut in range(8):
        for bit in (0,1):
            heights=[]
            for core in sides:
                blocks=[*core,*suffix]
                vals=[fold(blocks[:cut]+[tape(bit,i)]+blocks[cut:])[-1] for i in (1,2,3)]
                slope=(vals[1]-vals[0])//2
                constant=vals[0]-2*slope
                assert vals==[slope*2**i+constant for i in (1,2,3)]
                heights.append((slope,constant))
            (s,b),(S,B)=heights
            assert s==S
            result[cut,bit]=(B-b)//2
    print('Cuts 0..7, half-gaps for zero / one runs:',[(result[i,0],result[i,1]) for i in range(8)])
    assert [result[i,b] for i in (1,2) for b in (0,1)]==[1255,-1880,-1254,1881]
    print('Unit-preserving exterior lower residues modulo 5:',{k%5,(k+r)%5})


if __name__=='__main__':
    main()
