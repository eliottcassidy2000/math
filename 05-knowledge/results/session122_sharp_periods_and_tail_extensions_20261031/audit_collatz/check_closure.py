"""Independent exact phase/endpoint/boundary checks using the audited corridor."""
from fractions import Fraction as Q
from pathlib import Path
import json

def need(ok,msg):
    if not ok: raise RuntimeError(msg)

def content_ones(E,k):return sum(k//(2*5**j) for j in range(E))
def digits(n):
    total=0
    while n:n,d=divmod(n,5);total+=d
    return total

def main():
    controls=0
    for E in range(1,11):
        for ell in range(1,E+1):
            T=2*5**(ell-1);M=5**(E-ell+1)
            bs=set(range(180))|{M-1,M,M+1,3*M+4}|{5**j for j in range(E-ell+2)}
            ss={0,1,T//3,T//2,T-2,T-1}
            for b in sorted(bs):
                for s in ss:
                    if not 0<=s<T:continue
                    k=b*T+s
                    loss=5*(k-s)-8*(content_ones(E,k)-content_ones(E,s))
                    want=2*digits(b%M)+10*(b//M)
                    need(loss==want,'exact phase loss')
                    permitted=b==0 or b in {5**j for j in range(E-ell+1)}
                    need((loss<=2)==permitted,'complete deficit-two phase bank')
                    controls+=1
    root=Path(__file__).resolve().parent
    graph=json.loads((root/'corridor_E10.json').read_text())
    need((graph['nodes'],graph['edges'],graph['retained_entries'])==(3908743,3908754,10),'complete graph binding')
    expected={30:[(342,326)],31:[],32:[(x,x-16) for x in [72,352,392,552,592,1592,6592,31592,156592,781592,3906592]]}
    endpoint_checks=0; exact_multiplier_heights=[]
    for row in graph['records']:
        need([tuple(x) for x in row['pairs']]==expected[row['score']],'independent exact pair set')
        for A,L in row['pairs']:
            need(0<L<A and pow(3,L,5**10)!=pow(2,A,5**10),'actual endpoint exclusion')
            v=0
            while v<10 and pow(3,L,5**(v+1))==pow(2,A,5**(v+1)):v+=1
            need(v in (2,3),'near-tight forgotten multiplier digit')
            exact_multiplier_heights.append([A,L,v])
            endpoint_checks+=1
    need(endpoint_checks==12,'all near-tight pairs excluded')
    # Every omitted letter has slack greater than two, for entries and edges.
    need(-128-(-10)-80+325>2 and -128-80+325>2,'strict omitted-letter gap')
    eta=Q(183105426,48828125);fee=Q(24,25)*Q(33,8)+eta/25;alpha=fee-4
    lam=lambda e:Q(5,8)-Q(1,8*5**(e-1))
    need(Q(24,25)*Q(5,8)+lam(12)/25==lam(14),'same-word interpolation')
    need(fee==Q(5017089801,1220703125) and alpha==Q(134277301,1220703125)>0,'fee and margin')
    for e,value in [(11,Q(31,125)),(12,Q(6,25)),(13,Q(1,5))]:
        t=2*5**(e-1)
        need(lam(14)*t-content_ones(e,t)==value,'all-one primitive price')
    need(fee-2>0,'core invoice coefficient')
    tangent=0
    for x in range(5):
        for y in range(5):
            old=Q((-x)%5,5)
            for u in range(5):
                for v in range(5):
                    X,Y=(x+u)%5,(y+v)%5
                    bill=Q(int(X==0 and Y==0))+Q((-X)%5,5)-old
                    need(bill<=1,'general tangent payment')
                    if(u,v)==(4,4):need(bill<=Q(1,5),'all-one payment')
                    tangent+=1
    q=5**14;rho=3*pow(2,-1,q)%q
    need(pow(rho,2*5**12,q)==1+4*5**13,'native all-one tangent')
    need(Q(61035156,48828125)>alpha,'suffix margin')
    t=2*5**13;s=content_ones(14,t)
    need(Q(s,t)==lam(14)==Q(762939453,1220703125),'sharp all-one endpoint')
    return dict(status='PASS',exact_multiplier_heights=exact_multiplier_heights,phase_controls=controls,near_tight_clock_pairs_excluded=endpoint_checks,
        tangent_controls=tangent,integer_native_identity10_fee='33/8',interpolated_fee=str(fee),
        nonone_S14_margin=str(alpha),sharp_S14_slope=str(lam(14)),sharp_clock=t,sharp_charge=s,
        inherited_rational_certificate_rerun=False,all_layer_and_convergence='OPEN')

if __name__=='__main__':print(json.dumps(main(),sort_keys=True,indent=2))
