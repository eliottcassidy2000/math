"""Physical syzygy words preserving both fringes and the dyadic carry."""
from math import lcm
from replay_joint_phase import V, power, initial, direct, solve
from replay_joint_phase import step as original_step
from replay_dyadic import compiler

def order(b,m):
    x=b%m; k=1
    while x!=1%m:
        x=x*b%m;k+=1
    return k

def lift_word(a,m,h2,P2,h_m,P_m,h3):
    if m<2 or m%2==0 or m%3==0:raise ValueError(m)
    base=compiler(a,h2,P2)
    cut=len(base)-2*a-2
    pre,suffix=base[:cut],base[cut:]
    if pre[-a:]!='0'*a or suffix[:a]!='0'*a:raise ValueError('insufficient dyadic guards')
    lam=lcm(4,order(3,1<<a))
    a_ord=lcm(lam,2,order(2,m),order(3*pow(2,-1,m),m))
    guard='0'*(((a+1+a_ord-1)//a_ord)*a_ord)
    tokens=['0','10'+guard]
    def step(v,k):
        z=original_step(v[:6]+(v[6]%2,)+v[7:],m,k)
        return z[:6]+((v[6]+1+k)%lam,)+z[7:]
    def run(v,cmds,out=None):
        for s,d in cmds:
            if d<0:
                t=step(v,s);mag=1
                while t!=v:t=step(t,s);mag+=1
                d=(mag-1)*(-d)
            for _ in range(d):v=step(v,s)
            if out is not None:out.append(tokens[s]*d)
        return v
    u=initial(pre,m)
    start=u[:6]+(len(pre)%lam,)+u[7:]
    changes=[]
    for word in V:
        w=run(start,word)
        if w[:4]!=start[:4] or w[6:8]!=start[6:8]:raise ValueError('fixed state')
        changes.append(tuple((v-b)%m for v,b in zip(w[4:6],start[4:6])))
    h0,P0=direct(base,m)
    target=((h_m-h0)%m,(P_m-P0)*pow(pow(2,len(base),m),-1,m)%m)
    m0,m1=m,1
    while m0%5==0:m0//=5;m1*=5
    r,s=solve(*changes[:2],target,m0)
    t,u=solve(changes[0],changes[2],target,m1)
    x0=(r,s,0);x1=(t,0,u)
    exps=[(r+m0*((s-r)*pow(m0,-1,m1)%m1))%m if m1>1 else r for r,s in zip(x0,x1)]
    dx3=(h3-direct(base,3)[0])%3
    exps=[e+m*((want-e)*pow(m,-1,3)%3) for e,want in zip(exps,(0,0,dx3))]
    out=[pre];state=start
    for w,e in zip(V,exps):state=run(state,power(w,e),out)
    out.append(suffix)
    word=''.join(out)
    if len(word)%lam!=len(base)%lam or not word.startswith('1') or not word.endswith('0'):
        raise ValueError('shape')
    for mod,h,P in ((1<<a,h2,P2),(m,h_m,P_m),(3,h3,h3)):
        if direct(word,mod)!=(h%mod,P%mod):raise ValueError(('target',mod,h,P))
    return word

if __name__=='__main__':
    for a,m,p,q in ((4,5,13,37),(5,7,23,31),(6,5,17,25)):
        h,P=71*p-44*q,21*p-13*q
        hp3=(p%3 or 1)
        word=lift_word(a,m,h,P,h,P,hp3)
        print(a,m,p,q,len(word))
