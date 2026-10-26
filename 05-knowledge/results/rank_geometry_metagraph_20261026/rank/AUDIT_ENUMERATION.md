# Independent audit of semantic finite-sensor closure

2026-10-26. **ACCEPTED, with the source and domain restrictions retained.** Reviewer: recurrence/enumeration agent. Read `README.md` and `check.py`; replayed normal and optimized modes. Both review receipts match both producer receipts, SHA256 `cc61c221aa4ede11abb611387bb2ab9fbe9f1eed656a852f92c1b814641b596c`.

## The load-bearing quantifiers

- Theorem1 concerns **all positive integer sources in the exact legal cylinder of the chosen balanced block**, not just sources reachable from7. Initial credit can be fixed because balanced prefixes are executable at every positive credit and return that credit unchanged.
- All sensor centers lie in Z2: rational numbers with odd denominator. Their omitted pullback rho also lies in Z2 because its affine denominator is an odd power of3. Thus all finite separation depths are nonnegative integers, the binary CRT congruence is legitimate, and the cylinder test includes rho.
- If a pulled center is missing, choosing exact distance N from it beyond every retained separation depth forces every retained sensor to be constant. Choosing N beyond the fixed binary precision preserves **one compatible** residue class modulo M. The proof does not prescribe an arbitrary independent incompatible binary residue. The odd residue is arbitrary. This yields one observation fiber with arbitrarily large finite successor values, which defeats both exact prediction and every finite upper-bound function on that fiber.
- If every pulled center is present, the displayed affine identity gives the whole successor sensor vector directly. It need not give the successor residue modulo M. The theorem explicitly separates that stronger, generally false closure condition.
- On a finite nonempty set, an injective inverse map permutes the set. Its real contraction factor `(2/3)^L` then forces zero diameter. The only possible singleton is the negative periodic center. The common-center/primitive-word characterization is the inherited periodic-itinerary statement; it is not inferred merely from commuting scalar values. Empty sensor banks and empty block collections have correctly different scopes.

I found no counterexample to these quantified statements. In particular, positive or zero sensors do not defeat necessity; negative centers are allowed, infinity is accounted for in the exact identity, and restriction to positive integer sources does not prevent arbitrarily large CRT lifts.

## Actual-prefix construction

Theorem3 is also sound. Prescribing `P Q^m 0^47` puts the first mismatch from rho at position `3+6m`; its distance is therefore exactly that integer. Since the actual entry has depth12 from rho, for m>=2 the new entry differs from the actual entry at **exactly** depth12. The pullback through the 86-letter actual word adds86, giving depth98 at the original source. The odd congruence `3^86 M0` is sufficient even when M0 contains powers of3: it makes the original pullback an integer and preserves M0. Positivity follows by taking a large lift. The matched parity prefix preserves the native guards and credit history, not the numerical trajectory values.

At credit47 the first46 tail zeros execute and the47th is blocked. The clock is consequently `86+3+6m+46`, excluding the blocked symbol. All these examples are different starting integers and halt; neither they nor the semantic closure theorem decide the original source7. The exact loss of agreement at entry bit13 and source bit99 prevents treating these finite-prefix constructions as a compatible infinite approximation to that source.

## Independent finite controls

The executable block below imports no mathematical producer. It constructs binary cylinders by **literal low-bit lifting**, independently of the producer's affine modular-inverse constructor; builds the 71-word bank by a ballot-word filter instead of ternary-tree recursion; and tests source pullbacks at m2,5,9 with four odd sidecars. Its generic sensor hostile is constructed at the **successor first**, then pulled back through P, rather than copying the missing-center construction. All tests use explicit exceptions and print JSON. This is bounded verification, not the proof of the universal quantifiers above.

```python
from fractions import Fraction as F
from itertools import combinations
import json

def need(ok, why):
    if not ok:
        raise RuntimeError(why)

def val(x):
    x = F(x)
    if x == 0:
        return float('inf')
    a, b, result = abs(x.numerator), x.denominator, 0
    while a % 2 == 0:
        a //= 2
        result += 1
    while b % 2 == 0:
        b //= 2
        result -= 1
    return result

def native(z, q, length):
    word = ''
    for _ in range(length):
        bit = z % 2
        need(z >= 3 and q >= 1 and not (q == 1 and bit == 0), 'native guard')
        word += str(bit)
        z, q = (3*z + bit)//2, q + 3*bit - 1
    return z, q, word

def carry(w):
    return sum(int(c)*2**i*3**(len(w)-1-i) for i, c in enumerate(w))

def center(w):
    return -F(carry(w), 3**len(w)-2**len(w))

def lift_word(w):
    # Literal parity, including the formal last blocked input, determines
    # one new low bit at a time; this function does not use the carry.
    z = 0
    for i, wanted in enumerate(w):
        current = z
        for _ in range(i):
            current = (3*current + current % 2)//2
        if current % 2 != int(wanted):
            z += 2**i
    current = z
    for wanted in w:
        need(current % 2 == int(wanted), 'low-bit cylinder')
        current = (3*current + current % 2)//2
    return z

def crt(a, two, b, odd):
    z = a % two
    return z + two*((b-z)*pow(two, -1, odd) % odd), two*odd

words = []
for e in range(1, 5):
    for places in combinations(range(3*e), e):
        w = ''.join('1' if i in places else '0' for i in range(3*e))
        credit, good = 0, True
        for c in w:
            credit += 3*int(c)-1
            good = good and credit >= 0
        if good:
            words.append(w)
bank = sorted({F(0)} | {center(w) for w in words})
need(len(words) == 71 and len(bank) == 67, 'independent ballot bank')
p, qword = '100', '101000'
rp, rq = center(p), center(qword)
rho = (8*rq-9)/27
actual, credit, prefix = native(7, 1, 86)
need((actual, credit) == (10952205219944621, 47), 'actual state')
signature = tuple(val(actual-s) for s in bank)
need(max(val(rho-s) for s in bank) == 11 and val(actual-rho) == 12, 'separation')

generic = 0
for modulus in (1, 2**18*105):
    e = val(modulus)
    odd = modulus // 2**e
    observations = []
    for depth in (30, 31, 35):
        two = 2**(depth+1)
        binary_y = (rq.numerator*pow(rq.denominator, -1, two) + 2**depth) % two
        odd_y = (27*7+9)*pow(8, -1, 27*odd) % (27*odd)
        y, period = crt(binary_y, two, odd_y, 27*odd)
        y += period
        need((8*y-9) % 27 == 0, 'successor-first integer inverse')
        x = (8*y-9)//27
        need(native(x, 1, 3) == (y, 1, p), 'successor-first legal P')
        need(val(F(y)-rq) == depth, 'unbounded successor sensor')
        observations.append((x % modulus, tuple(val(F(x)-s) for s in bank)))
        generic += 1
    need(len(set(observations)) == 1, 'whole observation fiber')

families = []
for m in (2, 5, 9):
    future = p + qword*m + '0'*47
    two = 2**len(future)
    binary = lift_word(future)
    for odd_sidecar in (1, 9, 35, 5005):
        odd = 3**86*odd_sidecar
        entry, period = crt(binary, two, actual, odd)
        if entry <= actual:
            entry += ((actual-entry)//period+1)*period
        numerator = 2**86*entry-carry(prefix)
        need(numerator % 3**86 == 0, 'ordinary source pullback')
        source = numerator//3**86
        need(source > 7 and (source-7) % odd_sidecar == 0, 'source and odd sidecar')
        need(val(source-7) == 98 and val(entry-actual) == 12, 'sharp fixed binary boundary')
        need(val(F(entry)-rho) == 3+6*m, 'regenerated precision')
        need(tuple(val(F(entry)-s) for s in bank) == signature, 'all 67 sensors')
        need(native(source, 1, 86) == (entry, 47, prefix), 'actual parity and credit prefix')
        end, q, got = native(source, 1, len(prefix+future)-1)
        need(got == (prefix+future)[:-1] and q == 1 and end % 2 == 0, 'exact final blocked input')
        need(len(got) == 135+6*m, 'native execution clock')
        families.append({'m': m, 'odd_sidecar': odd_sidecar, 'executed': len(got)})
print(json.dumps({'status': 'PASS', 'ballot_words': len(words), 'sensors': len(bank),
                  'successor_first_hostiles': generic, 'actual_prefix_families': families}, sort_keys=True))
```

Extract the single Python block and run it normally and with `-O`. Review output is saved in `review-independent-normal.json` and `review-independent-optimized.json`. No production code or previous packet was changed by this review.
