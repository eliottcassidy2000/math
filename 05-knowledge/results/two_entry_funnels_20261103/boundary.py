"""Old-root growth of the five exceptional families and the singleton-core dual."""
from research import *
from verify import rooted_data
from carrier import field


def obs(T,r):
    p,e,s=bin_values(T,r,2)
    require((s[2]-e[1])%2==0,'observer parity')
    return p[1],(s[2]-e[1])//2,(s[2]+e[1])//2


def main():
    tails=[(),((1,0),),((0,1),),((1,0),(0,1)),((0,1),(1,0))]
    growth=0
    for tail in tails:
        for k in range(1,9):
            T,r=incoming_pair_source(((1,1),)*k+tail)
            U,x=birth(T,r)
            V,v=incoming_pair_source(((1,1),)*(k+1)+tail)
            order=[0,1,x]+list(range(2,r))+[r]
            require(induced(U,order)==V,'old-root growth prepends exactly one 11')
            growth+=1
    ladder=[]
    for k in range(1,9):
        T,r=incoming_pair_source(((1,1),)*k)
        A,a=funnel((0,),{0},(1,)*k);A,b=birth(A,a)
        require(induced(opposite(A),[0,b]+list(range(k,0,-1))+[a])==T,'singleton-core converse transport')
        R=opposite(T)
        va=obs(R,1);vb=obs(T,r);h=3*2**k-1
        require(rooted_data(R,1)[3]-h==1 and rooted_data(T,r)[3]-h==2,'both ready lanes')
        xa=3**k+3*2**(k-1)-2;xb=2**(k+1)-2
        require(va==(h,xa,xa+1) and vb==(h,xb,xb+1),'full observer formulas')
        require((field(va)==field(vb))==(k%6==2),'exact converse phase clock')
        ra,_=reset(R,1);rb,_=reset(T,r)
        require(tables(ra)[0][-1]==9*2**k-1 and tables(rb)[0][-1]==9*2**k+1,'two literal reset values')
        ladder.append(dict(k=k,n=len(T),H=h,lane1_observer=va,lane2_observer=vb,
                           lane1_phase=field(va),lane2_phase=field(vb),branches=[9*2**k-1,9*2**k+1]))
        tables.cache_clear()
    print(json.dumps(dict(status='PASS',literal_exception_growth=growth,singleton_dual_ladder=ladder),indent=2,sort_keys=True))


if __name__=='__main__':main()
