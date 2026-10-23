"""Independent symmetric congruence over Z/2^precision using 1x1 and 2x2 blocks."""
def valuations(A,precision=40):
    A=[[int(x)%(1<<precision) for x in row] for row in A]
    answer=[];shift=0;mod=1<<precision
    while A:
        n=len(A)
        one=next((i for i in range(n) if A[i][i]%2),None)
        if one is not None:
            idx=[one]+[i for i in range(n) if i!=one]
            A=[[A[i][j] for j in idx] for i in idx];inv=pow(A[0][0],-1,mod)
            A=[[(A[i][j]-A[i][0]*inv*A[0][j])%mod for j in range(1,n)] for i in range(1,n)]
            answer.append(shift)
            continue
        two=next(((i,j) for i in range(n) for j in range(i+1,n) if A[i][j]%2),None)
        if two is not None:
            idx=list(two)+[i for i in range(n) if i not in two]
            A=[[A[i][j] for j in idx] for i in idx];a,b,d=A[0][0],A[0][1],A[1][1]
            inv=pow((a*d-b*b)%mod,-1,mod)
            A=[[(A[i][j]-inv*(d*A[i][0]*A[j][0]-b*(A[i][0]*A[j][1]+A[i][1]*A[j][0])+a*A[i][1]*A[j][1]))%mod for j in range(2,n)] for i in range(2,n)]
            answer.extend([shift,shift])
            continue
        if mod<=2:raise RuntimeError('Precision exhausted, no complete result')
        A=[[x//2 for x in row] for row in A];mod//=2;shift+=1
    return sorted(answer)
