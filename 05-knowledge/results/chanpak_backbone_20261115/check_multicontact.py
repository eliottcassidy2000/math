from replay_rooted_clock import all_h,forward,reverse

def ext(A,m,ports):
    n=len(A); k=len(ports)
    R=forward(n+m+k)
    for i,row in enumerate(A):
        R[i]=(R[i]&~((1<<n)-1))|row
    for i,r in enumerate(ports):reverse(R,r,n+m+i)
    return R
if __name__=='__main__':
    A=[2,4,1]
    for m in range(1,9):
        for ports in ((0,1,2),(0,2,1)):
            b=sum(bool(A[ports[i]]>>ports[j]&1) for i in range(3) for j in range(i))
            H=all_h(ext(A,m,ports))[-1]
            remainder=H-(6*4**m+24*3**m+(12+b)*2**m)
            if remainder != -7:
                raise RuntimeError(m,ports,b,H,remainder)
    print("literal three-contact expansion: all 16 controls agree")
