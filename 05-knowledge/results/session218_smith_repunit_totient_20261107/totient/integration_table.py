#!/usr/bin/env python3
"""Exact bounded integration; inverse source enumeration imported from audited core."""
import contextlib, io, json
with contextlib.redirect_stdout(io.StringIO()):
    from verify import inv
from sympy import factorint,isprime,totient
from pathlib import Path

def digits(n,b):
    out=[]
    while n:
        n,r=divmod(n,b);out.append(r)
    return out[::-1] or [0]
rows=[]
for n in range(2,77,2):
    if inv(n):continue
    factors={int(p):int(e) for p,e in factorint(n).items()}
    F=sum(p*e for p,e in factors.items());b=n-F+1
    ds=digits(n,b);w=sum(e*sum(digits(p,b)) for p,e in factors.items())
    assert b>=2 and sum(ds)==w
    repaired=None
    for a in range(1,9):
        pre=inv(n*2**a)
        if pre:
            repaired={'a':a,'target':n*2**a,'preimage':min(pre)}
            assert totient(min(pre))==n*2**a
            break
    rows.append({'n':n,'factorization':factors,'F':F,'base':b,'digits':ds,'factor_digit_weight':w,
                 'prime_half':bool(isprime(n//2)),'phi_n':int(totient(n)),
                 'dyadic_repair':repaired or 'not found for 1<=a<=8'})
assert [r['n'] for r in rows]==[14,26,34,38,50,62,68,74,76]
text=['# Even nontotients through 76: exact integration','',
'Status: FINITE-EXACT. The universe is every even integer 2 through 76, with complete inverse-totient enumeration filtering empty fibers. F counts prime factors with multiplicity. Base digits are integer arrays, not decimal strings. Repair a means phi(x)=2^a n; every smaller a was checked, so a found within 1..8 is the exact first positive repair. A failure through 8 would only be bounded-not-found. No unbounded existence claim is made.','',
'| n | factorization | F(n) | Smith base n−F+1 | digits | factor-digit weight | prime half? | phi(n) | first repair a | preimage x |',
'|---|---|---|---|---|---|---|---|---|---|']
for r in rows:
    fac=' * '.join(str(p)+(f'^{e}' if e>1 else '') for p,e in r['factorization'].items())
    d=r['dyadic_repair']
    a=d['a'] if isinstance(d,dict) else 'none <=8'; x=d['preimage'] if isinstance(d,dict) else 'none <=8'
    text.append(f"| {r['n']} | {fac} | {r['F']} | {r['base']} | {r['digits']} | {r['factor_digit_weight']} | {r['prime_half']} | {r['phi_n']} | {a} | {x} |")
text.extend(['','The moving-base Smith theorem of session218-smith/BASE_COMPLETENESS.md applies to every even nontotient: each is composite (2 is a totient and other even integers are composite), and the even excluded Smith numbers 6,8,12 are themselves totients. This does not characterize nontotients or imply a common-base Smith property. The phi(n) column concerns n as input; the nontotient filter concerns n as output.','',
'Reproduce: python3 session218-scaling/integration_table.py. Uses exact prime-support inverse enumeration from verify.py, independently checked there against a complete totient sieve for targets <=200. Larger repair targets retain the same proven-complete source enumeration, with explicit preimages checked independently by Euler’s totient formula.',''])
root=Path(__file__).parent
(root/'INTEGRATION_TABLE.md').write_text('\n'.join(text))
(root/'integration_receipt.json').write_text(json.dumps({'status':'PASS','universe':'even 2..76','repair_range':'a=1..8','rows':rows},indent=2)+'\n')
print('\n'.join(text))
