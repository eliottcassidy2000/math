from functools import lru_cache

def f(z):
    if not z:return 1
    return f(z//2) if z&1 else f(z//2)+z//2+(z//2&1)

@lru_cache(None)
def inverse(N,a,b,c):
    """All z>=0 such that a*z+b*f(z)+c*(z&1) == N."""
    assert a>0 and b>0 and c>=0
    if N<b: return ()
    ret=[0] if N==b else []
    if (N-b)&1==0: # even z; retain u≥1 to avoid z=0 cycling at a fixed target
        if N>=2*a+b:
            ret.extend(2*u for u in inverse(N,2*a+b,b,b) if u)
    if ((N-a-c-b)&1)==0 and N>=a+b+c:
        ret.extend(2*u+1 for u in inverse(N-a-c,2*a,b,0))
    return tuple(sorted(ret))

def rc(s):
    return s[::-1].translate(str.maketrans('01','10'))

def select(p):
    answers=[]
    if p%3==2:
        for z in inverse((p-2)//3,1,1,0):
            if z>=2 and z%2==0:
                answers.append((rc(bin(z)[2:]),0,z))
    for z in inverse(p-4,6,5,2):
        if z>=1:answers.append(('1'+rc(bin(z)[2:]),1,z))
    return answers

if __name__=='__main__':
    from verify import word_rows, join, count
    for a,b,c in ((1,1,0),(6,5,2),(2,9,7),(1,2,1)):
        bank={}
        for z in range(601):
            v=a*z+b*f(z)+c*(z&1)
            if v<=600:bank.setdefault(v,[]).append(z)
        for N in range(601):
            assert inverse(N,a,b,c)==tuple(bank.get(N,[]))
    assert select(30011)==[('1001100111100',0,6246),('1000011001100',0,6558)]
    assert select(30013)==[('10001000110110',1,4727)]
    for p in (30011,30013):
        out=select(p)
        for w,y,z in out: assert count(join((2,4,1),word_rows(w),(0,y)))==p
        print(p,out)
