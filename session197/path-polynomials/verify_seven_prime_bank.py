#!/usr/bin/env python3
from verify_inverse_bank import *

def run():
    forms=[(3,6),(5,12),(5,18),(11,36),(13,36),(25,72),(19,60),(39,120)]
    candidates=[(a,b,h,u) for a,b in forms for h in range(1,8,2) for u in range(8) if a*h+b*u==7]
    require(not candidates,'seven last-generator obstruction')
    bound=7*127;N=bound.bit_length()-1+5
    require(N==14,'atom cutoff')
    levels=make_profiles(N,bound)
    require(sum(map(len,levels))==1745,'profile census')
    attainable={value(p,1) for level in levels for p in level}
    targets=[7*p for p in primes(127)]
    hits=[h for h in targets if h in attainable]
    missing=[h for h in targets if h not in attainable]
    require(hits==[49,91,133,329,749,763],'exact 7p coverage')
    require(len(targets)==30 and len(missing)==24,'target census')
    x=lambda a:join((0,),a)
    f=lambda a:port_wrapper((4,1,2),(0,1),a)
    g=lambda a:port_wrapper((0,1,3,7),(0,2),a)
    witnesses=[(49,'F(Tr3)',f(transitive(3))),
       (91,'E(D(Tr2))',edge_wrapper(wrapper(transitive(2)))),
       (133,'D(X^2 D(Tr2))',wrapper(x(x(wrapper(transitive(2)))))),
       (329,'D(F(Tr2))',wrapper(f(transitive(2)))),
       (749,'E(G(K1))',edge_wrapper(g((0,)))),
       (763,'E(X^3 D(Tr2))',edge_wrapper(x(x(x(wrapper(transitive(2)))))))]
    output=[]
    for h,recipe,a in witnesses:
        validate(a);require(h_subsets(a)[-1]==h,('physical endpoint DP',h,recipe,a))
        p=tuple(bin_polynomial(a));require(p in levels[len(a)],('full-profile membership',h,p))
        output.append({'H':h,'recipe':recipe,'rows':a,'order':len(a)})
    print(json.dumps({'prime_multiplier_bound':127,'target_bound':bound,'complete_order_bound':N,'profiles_per_order':[len(x) for x in levels[1:]],'total_profiles':sum(map(len,levels)),'targets':targets,'attained':hits,'omitted':missing,'physical_witnesses':output},sort_keys=True,indent=2))

if __name__=='__main__':run()
