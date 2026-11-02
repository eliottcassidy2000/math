#!/usr/bin/env python3
"""Exact endpoint audit for PROOF_N748.md. No numerical optimization."""
from fractions import Fraction as F
from math import factorial
import json
N=748
c=F(249999,1000000)
eta=F(1,55)
delta=F(1,3500)
k=F(8,N-1)
C=4/(c-4*eta-delta)
h=eta-C*eta*eta
checks={}
def pos(name,value):
    if not value>0:
        raise RuntimeError((name,value))
    checks[name]={"positive_exact":str(value)}
def yes(name,value):
    if not value:
        raise RuntimeError(name)
    checks[name]=True
s=sum((F(4)**j)/factorial(j) for j in range(8))
yes('partial_exp4_value',s==F(16319,315))
pos('partial_exp4_exceeds50',s-50)
pos('c0_lower_exceeds_kappa_upper',F(49,200)-k)
pos('V_lower_exceeds4',F(49*(N-1),1600)-4)
yes('exp16_lower_exceeds250000',50**4>250000)
yes('c_value_from_250000',c==F(1,4)*(1-F(1,250000)))
pos('variance_upper_n_over7',F(N,7)-(F(N,8)+F(4*N,N-1)))
yes('variance_bound_uniform_condition',N-1>=224)
eta_count_mass=4*(c*(N-1)/8-F(4*N,N-1)-eta*F(N,2))-F(N,15)
pos('heavy_count_lower_minus_n_over15',eta_count_mass)
yes('heavy_count_endpoint_value',eta_count_mass==F(428091991,1494000000))
pos('heavy_count_derivative_constant',c/2-2*eta-F(1,15))
yes('exponential_endpoint_integer',(2000*N)**15*3**N<4**N)
pos('exponential_monotonicity',F(1,60)-F(1,N))
yes('epsilon_variance_error',4*F(1,2000)*F(1,7)==delta)
pos('light_count_denominator',c-4*eta-delta)
pos('shared_budget_interior',F(1,2)-C*eta)
yes('shared_budget_interior_value',F(1,2)-C*eta==F(2427923,27255846))
pos('shared_budget_contradiction',h-k)
yes('shared_budget_contradiction_value',h-k==F(190787,186634405485))
pos('heavy_weights_exceed_kappa',eta-k)
print(json.dumps({'status':'PASS','N':N,'constants':{'c':str(c),'eta':str(eta),'delta':str(delta),'C':str(C),'h':str(h)},'checks':checks},sort_keys=True,indent=2))
