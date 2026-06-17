import sys
sys.path.insert(0,'/tmp')
from lrc_base import *

# At each safe point of 1..13, WHICH speeds achieve the min 1/14?
S=list(range(1,14))
for t in [F(1,14),F(3,14),F(5,14)]:
    achievers=[v for v in S if nrm(v*t)==F(1,14)]
    print(f"tau={t}: min achieved by speeds {achievers}, sum pairs:",
          [(a,b,a+b) for i,a in enumerate(achievers) for b in achievers[i+1:]])

# Mechanism: tau = k/14 (odd k). nrm(v*k/14)=1/14 iff v*k ≡ ±1 (mod 14).
# For k=1: v ≡ ±1 mod 14 -> v in {1,13}. 
# For k=3: 3v ≡ ±1 mod 14 -> v ≡ ±5 mod 14 (3*5=15≡1) -> v in {5,9}.
# For k=5: 5v ≡ ±1 mod 14 -> v ≡ ±3 mod 14 (5*3=15≡1) -> v in {3,11}.
print("\nMechanism check (mod 14):")
for k in [1,3,5]:
    inv=pow(k,-1,14)
    print(f"  k={k}: v ≡ ±{inv} mod 14 -> achievers in 1..13: {[v for v in range(1,14) if v%14 in (inv, 14-inv)]}")
# So each safe point is propped up by a complementary pair (v, 14-v) both in the set.
# To DESTROY tau=k/14, must remove BOTH speeds of every propping pair -- but removing
# breaks other points. This is the combinatorial obstruction.
print("\nThe pairs {1,13},{5,9},{3,11} prop up tau=1/14,3/14,5/14 respectively.")
print("Also {2,12},{4,10},{6,8},{7} give norms >1/14 at these tau (not binding).")
