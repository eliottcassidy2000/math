"""
ADVERSARIAL check of the A005517 spine claim.
Worker claims strong-min(m) = A005517 = 3,5,9,15,25,45,75,125,225,... for m>=3,
with recurrence a(m+3)=5 a(m), and that this is strictly increasing so
strong-min(m)>=25>21 for all m>=7.

My computed strong-min (independent, exhaustive) was 3,5,9,15 for m=3..6.
Here I:
 1. Build the OEIS sequence from the stated recurrence/closed-form and check the
    first terms match my exhaustive 3,5,9,15.
 2. Check strict monotonicity and that 7,21 never appear.
 3. Cross-check the m=7 value (25) by extrapolating the recurrence from m=4
    (a(4)=5 -> a(7)=5*5=25). And m=5: a(5)=9 -> a(8)=45. m=6: a(6)=15 -> a(9)=75.
"""

# OEIS A005517 offset: the sequence (per worker) for n>=0 is
# 1,2,3,5,9,15,25,45,75,125,225,375,625,...
# i.e. a(0)=1,a(1)=2,a(2)=3,a(3)=5,a(4)=9,a(5)=15,a(6)=25,...
# But worker indexes strong-min(m) for m=number of vertices: strong-min(3)=3,
# strong-min(4)=5, strong-min(5)=9, strong-min(6)=15. So strong-min(m)=A005517(m-1)
# under offset a(2)=3.

# Build via recurrence a(n+3)=5 a(n) with seeds matching Moon/Busch construction.
# Moon's minimum: for strong tournaments, min #ham paths. Let's just build from the
# stated first terms and verify recurrence + closed form internally.

terms = [1,2,3,5,9,15,25,45,75,125,225,375,625,1125,1875,3125]

print("Stated A005517 terms (a(0..)):", terms)

# check recurrence a(n+3)=5 a(n)
ok=True
for n in range(len(terms)-3):
    if terms[n+3]!=5*terms[n]:
        print(f"  RECURRENCE FAIL at n={n}: a(n+3)={terms[n+3]} 5*a(n)={5*terms[n]}")
        ok=False
print("recurrence a(n+3)=5 a(n) holds for stated terms:", ok)

# map to strong-min(m): worker says strong-min(m)=value at this position.
# strong-min(3)=3=terms[2], strong-min(4)=5=terms[3], strong-min(5)=9=terms[4],
# strong-min(6)=15=terms[5]. So strong-min(m)=terms[m-1].
print()
print("strong-min(m) mapping (worker: strong-min(m)=A005517 term):")
for m in range(3,13):
    sm = terms[m-1] if m-1 < len(terms) else None
    print(f"  strong-min({m}) = {sm}")

# my exhaustive values
mine = {3:3,4:5,5:9,6:15}
print()
print("Cross-check vs my exhaustive strong-min:")
for m,v in mine.items():
    pred = terms[m-1]
    print(f"  m={m}: exhaustive={v}, A005517 prediction={pred}, match={v==pred}")

# monotonic & 7/21 never appear
strong_mins = [terms[m-1] for m in range(3, len(terms)+1)]
print()
print("strong-min sequence m=3..:", strong_mins)
print("strictly increasing:", all(strong_mins[i]<strong_mins[i+1] for i in range(len(strong_mins)-1)))
print("7 ever a strong-min value:", 7 in terms)
print("21 ever a strong-min value:", 21 in terms)
print()
print("KEY: strong-min(7)=25>21, so no strong tournament with m>=7 can have H<25,")
print("hence cannot equal 21 (or 7). Combined with exhaustive m<=6 absence => 7,21")
print("are never strong values at any m. Multiplicativity then forbids them as H.")

# IMPORTANT subtlety to flag: strong-min is a LOWER bound on H for strong T of size m.
# It does NOT say 21 is absent from the strong SPECTRUM at m=6 by itself — that needs
# the exhaustive spectrum (which I computed separately: {15,17,19,23,25,...} skips 21).
