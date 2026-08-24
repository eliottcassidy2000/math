#!/usr/bin/env python3
"""Verify and search Epoch's imaginary-quadratic class-rank certificates.

Certificate grammar::

    ell|r|D|a1,b1,c1|...|ar,br,cr

The ``verify`` command is self-contained except for square-free factorization,
for which it uses python-flint (and falls back to SymPy).  Class arithmetic is
implemented directly with Buell's general Gauss-composition formula and exact
positive-definite reduction.  Independence is checked meet-in-the-middle in
O(ell**ceil(r/2)) stored classes rather than by enumerating all ell**r
relations.

The optional PARI audit deliberately uses a second implementation.  The
``search`` command requires a PARI/GP executable and extracts order-ell forms
from PARI's invariant-factor generators whenever a target is reached.

Examples::

    python3 imaginary_quadratic_class_rank_certificate_tool_20260824.py \
      verify '3|3|4447704|390,-96,2857|921,-786,1375|346,-68,3217'

    python3 imaginary_quadratic_class_rank_certificate_tool_20260824.py \
      verify --pari-gp /path/to/gp '3|3|4447704|...'

    python3 imaginary_quadratic_class_rank_certificate_tool_20260824.py \
      search --pari-gp /path/to/gp --ell 3 --target-rank 9 \
      --candidate 4447704 \
      --candidate 217541503961543485618350976479
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
from math import gcd, isqrt
from pathlib import Path
import subprocess
import sys
from typing import Dict, Iterable, List, Optional, Sequence, Tuple


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def xgcd(a: int, b: int) -> Tuple[int, int, int]:
    """Return positive g and x,y with ax+by=g."""
    old_r, r = a, b
    old_s, s = 1, 0
    old_t, t = 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r - q * r
        old_s, s = s, old_s - q * s
        old_t, t = t, old_t - q * t
    if old_r < 0:
        old_r, old_s, old_t = -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def xgcd3(a: int, b: int, c: int) -> Tuple[int, int, int, int]:
    """Return positive g and x,y,z with ax+by+cz=g."""
    g, x, y = xgcd(a, b)
    h, s, z = xgcd(g, c)
    return h, x * s, y * s, z


@dataclass(frozen=True, order=True)
class Form:
    a: int
    b: int
    c: int

    @property
    def discriminant(self) -> int:
        return self.b * self.b - 4 * self.a * self.c

    def validate(self, expected_discriminant: Optional[int] = None) -> None:
        delta = self.discriminant
        require(self.a > 0, f"form {self} is not positive definite: a <= 0")
        require(delta < 0, f"form {self} is not positive definite: disc >= 0")
        require(gcd(self.a, gcd(abs(self.b), self.c)) == 1,
                f"form {self} is not primitive")
        if expected_discriminant is not None:
            require(delta == expected_discriminant,
                    f"form {self} has discriminant {delta}, expected "
                    f"{expected_discriminant}")

    def reduced(self) -> "Form":
        """Return the unique properly reduced representative of this class."""
        a, b, c = self.a, self.b, self.c
        delta = self.discriminant
        require(a > 0 and delta < 0,
                f"reduction needs a positive-definite form, got {self}")
        for _ in range(1_000_000):
            # The unique q makes -a < b+2aq <= a.
            q = (a - b) // (2 * a)
            if q:
                b += 2 * a * q
                numerator = b * b - delta
                require(numerator % (4 * a) == 0,
                        "nonintegral coefficient during reduction")
                c = numerator // (4 * a)
            if a > c:
                # Proper substitution (x,y) -> (-y,x).
                a, b, c = c, -b, a
                continue
            if abs(b) > a:
                continue
            if b < 0 and (abs(b) == a or a == c):
                # Canonicalize the two reduced boundary identifications.
                b = -b
            out = Form(a, b, c)
            require(out.discriminant == delta,
                    "reduction changed the discriminant")
            require(abs(b) <= a <= c,
                    f"internal reduction failure: {out}")
            require(not (b < 0 and (abs(b) == a or a == c)),
                    f"noncanonical reduced boundary: {out}")
            return out
        raise ArithmeticError("positive-definite reduction did not terminate")

    def inverse(self) -> "Form":
        return Form(self.a, -self.b, self.c).reduced()

    def compose(self, other: "Form") -> "Form":
        """General Gauss composition, followed by exact reduction.

        This is Buell's formula.  Unlike the simplified Dirichlet formula, it
        also handles non-coprime leading coefficients and hence squares of the
        sample forms.
        """
        delta = self.discriminant
        require(other.discriminant == delta,
                "cannot compose forms of different discriminants")
        a1, b1, _ = self.a, self.b, self.c
        a2, b2, _ = other.a, other.b, other.c
        require((b1 + b2) % 2 == 0,
                "same-discriminant forms unexpectedly have different parity")
        beta = (b1 + b2) // 2
        n, t, u, v = xgcd3(a1, a2, beta)
        require(n > 0, "zero content in composition")
        require(a1 % n == 0 and a2 % n == 0,
                "composition gcd does not divide leading coefficients")
        require((a1 * a2) % (n * n) == 0,
                "composition leading coefficient is nonintegral")
        A = (a1 * a2) // (n * n)
        half = b1 * b2 + delta
        require(half % 2 == 0, "composition half-coefficient is nonintegral")
        numerator = a1 * b2 * t + a2 * b1 * u + v * (half // 2)
        require(numerator % n == 0,
                "composition middle coefficient is nonintegral")
        B = numerator // n
        tail = B * B - delta
        require(tail % (4 * A) == 0,
                "composition trailing coefficient is nonintegral")
        C = tail // (4 * A)
        out = Form(A, B, C).reduced()
        require(out.discriminant == delta,
                "composition changed the discriminant")
        return out

    def power(self, exponent: int) -> "Form":
        require(exponent >= 0, "negative exponents are not needed here")
        result = principal_form(self.discriminant)
        base = self.reduced()
        e = exponent
        while e:
            if e & 1:
                result = result.compose(base)
            e >>= 1
            if e:
                base = base.compose(base)
        return result


def principal_form(delta: int) -> Form:
    require(delta < 0, "principal form needs a negative discriminant")
    if delta % 4 == 1:
        return Form(1, 1, (1 - delta) // 4)
    require(delta % 4 == 0, f"{delta} is not a quadratic-order discriminant")
    return Form(1, 0, -delta // 4)


def factor_integer(n: int) -> List[Tuple[int, int]]:
    """Factor |n| exactly, preferring python-flint."""
    n = abs(n)
    require(n > 0, "cannot factor zero")
    try:
        from flint import fmpz  # type: ignore
        return [(int(p), int(e)) for p, e in fmpz(n).factor()]
    except ImportError:
        try:
            from sympy import factorint  # type: ignore
        except ImportError as exc:
            raise RuntimeError(
                "fundamental-discriminant checking needs python-flint or SymPy"
            ) from exc
        return sorted((int(p), int(e)) for p, e in factorint(n).items())


def is_squarefree(n: int) -> Tuple[bool, List[Tuple[int, int]]]:
    factors = factor_integer(n)
    return all(e == 1 for _, e in factors), factors


def fundamental_discriminant_data(D: int) -> Tuple[int, List[Tuple[int, int]]]:
    """Validate that -D is fundamental; return square-free core and factors."""
    require(D > 0, "D must be positive")
    delta = -D
    if delta % 4 == 1:
        core = delta
    elif delta % 4 == 0 and (delta // 4) % 4 in (2, 3):
        core = delta // 4
    else:
        raise ValueError(f"{delta} is not a fundamental discriminant")
    squarefree, factors = is_squarefree(core)
    require(squarefree,
            f"{delta} is not fundamental: square-free core {core} factors "
            f"as {factors}")
    return core, factors


def is_prime(n: int) -> bool:
    if n < 2:
        return False
    p = 2
    while p * p <= n:
        if n % p == 0:
            return n == p
        p = 3 if p == 2 else p + 2
    return True


def parse_form(text: str) -> Form:
    pieces = text.split(",")
    require(len(pieces) == 3, f"bad form field {text!r}")
    try:
        return Form(*(int(x.strip()) for x in pieces))
    except ValueError as exc:
        raise ValueError(f"bad integer in form field {text!r}") from exc


def parse_certificate(line: str) -> Tuple[int, int, int, List[Form]]:
    fields = line.strip().split("|")
    require(len(fields) >= 3, "certificate needs at least ell|r|D")
    try:
        ell, rank, D = map(int, fields[:3])
    except ValueError as exc:
        raise ValueError("ell, r, and D must be integers") from exc
    forms = [parse_form(field) for field in fields[3:]]
    require(len(forms) == rank,
            f"declared r={rank}, but found {len(forms)} form fields")
    return ell, rank, D, forms


Exponent = Tuple[int, ...]


def subgroup_table(
    generators: Sequence[Form], ell: int, identity: Form
) -> Tuple[Dict[Form, Exponent], Optional[Tuple[Exponent, Exponent, Form]]]:
    """Enumerate a small elementary-ell subgroup, retaining a collision."""
    table: Dict[Form, Exponent] = {identity: ()}
    for generator in generators:
        powers = [identity]
        for _ in range(1, ell):
            powers.append(powers[-1].compose(generator))
        next_table: Dict[Form, Exponent] = {}
        for value, exponent_prefix in table.items():
            for exponent, power in enumerate(powers):
                product = value.compose(power)
                vector = exponent_prefix + (exponent,)
                if product in next_table:
                    return next_table, (next_table[product], vector, product)
                next_table[product] = vector
        table = next_table
    return table, None


@dataclass(frozen=True)
class Verification:
    ell: int
    rank: int
    D: int
    factors: Tuple[Tuple[int, int], ...]
    reduced_forms: Tuple[Form, ...]
    left_size: int
    right_size: int


def verify_certificate(line: str) -> Verification:
    ell, rank, D, forms = parse_certificate(line)
    require(is_prime(ell), f"ell={ell} is not prime")
    require(rank >= 1, "r must be positive")
    _, factors = fundamental_discriminant_data(D)
    delta = -D
    reduced: List[Form] = []
    for form in forms:
        form.validate(delta)
        reduced.append(form.reduced())
    identity = principal_form(delta)

    for index, form in enumerate(reduced, start=1):
        require(form != identity, f"generator {index} is principal")
        require(form.power(ell) == identity,
                f"generator {index} does not have order ell={ell}")

    split = rank // 2
    left, collision = subgroup_table(reduced[:split], ell, identity)
    if collision is not None:
        raise ValueError(
            f"left-half dependence: exponents {collision[0]} and "
            f"{collision[1]} give {collision[2]}"
        )
    right, collision = subgroup_table(reduced[split:], ell, identity)
    if collision is not None:
        raise ValueError(
            f"right-half dependence: exponents {collision[0]} and "
            f"{collision[1]} give {collision[2]}"
        )
    require(len(left) == ell ** split, "left subgroup has the wrong size")
    require(len(right) == ell ** (rank - split),
            "right subgroup has the wrong size")

    intersections: List[Tuple[Exponent, Exponent, Form]] = []
    for right_value, right_exponents in right.items():
        left_exponents = left.get(right_value.inverse())
        if left_exponents is not None:
            intersections.append((left_exponents, right_exponents, right_value))
    require(len(intersections) == 1,
            "meet-in-the-middle found a nonzero relation: "
            f"{intersections[:2]}")
    only = intersections[0]
    require(not any(only[0]) and not any(only[1]),
            f"the unique half-intersection is not the zero relation: {only}")
    return Verification(
        ell=ell,
        rank=rank,
        D=D,
        factors=tuple(factors),
        reduced_forms=tuple(reduced),
        left_size=len(left),
        right_size=len(right),
    )


def require_gp_success(output: str, prefix: str, context: str) -> None:
    """Reject GP's nonfatal error stream and require a real success row.

    PARI/GP can exit with status zero after an interactive-style error and can
    echo source containing the requested marker.  A substring check is thus
    not a certificate.  Stack-size warnings are harmless; other ``***`` rows
    are treated as errors.
    """
    errors = [
        row for row in output.splitlines()
        if row.lstrip().startswith("***") and "Warning:" not in row
    ]
    successes = [
        row.strip() for row in output.splitlines()
        if row.strip().startswith(prefix)
    ]
    require(not errors, f"PARI/GP {context} reported errors:\n{output}")
    require(len(successes) == 1,
            f"PARI/GP {context} needs exactly one {prefix!r} row:\n{output}")


def pari_audit(line: str, gp_path: str) -> str:
    ell, rank, D, forms = parse_certificate(line)
    literal = ",".join(f"Qfb({f.a},{f.b},{f.c})" for f in forms)
    split = rank // 2
    program = f"""
default(parisizemax,2000000000); allocatemem(128000000);
ell={ell}; D={D}; F=[{literal}]; P=F[1]^0;
span(V,l,P)={{my(S=Set([P]),T,pw);for(i=1,#V,T=List();for(j=1,#S,pw=P;for(e=0,l-1,listput(T,S[j]*pw);pw=pw*V[i]));S=Set(T));return(S)}};
for(i=1,#F,if(F[i]^ell!=P,error("bad order"));if(F[i]==P,error("principal generator")));
k={split}; LF=if(k,F[1..k],[]); RF=if(k<#F,F[k+1..#F],[]);
L=span(LF,ell,P); R=span(RF,ell,P);
RI=Set(vector(#R,i,R[i]^-1)); Inter=setintersect(L,RI);
if(#L!=ell^k,error("left dependence"));
if(#R!=ell^(#F-k),error("right dependence"));
if(#Inter!=1,error("cross dependence"));
q=quadclassunit(-D); rk=sum(i=1,#q.cyc,valuation(q.cyc[i],ell)>0);
if(rk<#F,error("class group rank too small"));
print("PARI_AUDIT=PASS|h=",q.no,"|cyc=",q.cyc,"|ell_rank=",rk,"|halves=",#L,",",#R);
quit
"""
    process = subprocess.run(
        [gp_path, "-fq"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    require(process.returncode == 0,
            f"PARI/GP exited {process.returncode}:\n{process.stdout}")
    require_gp_success(process.stdout, "PARI_AUDIT=PASS|", "audit")
    return process.stdout.strip()


def read_certificate_lines(args: argparse.Namespace) -> List[str]:
    lines: List[str] = []
    if args.line:
        lines.extend(args.line)
    if args.file:
        for raw in Path(args.file).read_text(encoding="utf-8").splitlines():
            raw = raw.strip()
            if raw and not raw.startswith("#"):
                lines.append(raw)
    require(bool(lines), "provide a certificate line or --file")
    return lines


def run_verify(args: argparse.Namespace) -> int:
    for line in read_certificate_lines(args):
        result = verify_certificate(line)
        print(
            "PURE_PYTHON_VERIFY=PASS"
            f"|ell={result.ell}|r={result.rank}|D={result.D}"
            f"|squarefree_core_factors={list(result.factors)}"
            f"|mitm_halves={result.left_size},{result.right_size}"
        )
        print("REDUCED_FORMS=" + "|".join(
            f"{f.a},{f.b},{f.c}" for f in result.reduced_forms
        ))
        if args.pari_gp:
            print(pari_audit(line, args.pari_gp))
    return 0


def candidate_list(args: argparse.Namespace) -> List[int]:
    candidates = list(args.candidate or [])
    if args.candidates_file:
        for raw in Path(args.candidates_file).read_text(encoding="utf-8").splitlines():
            raw = raw.split("#", 1)[0].strip()
            if raw:
                candidates.append(int(raw))
    return candidates


def search_program(args: argparse.Namespace, candidates: Sequence[int]) -> str:
    ell = args.ell
    target = args.target_rank
    require(is_prime(ell), f"ell={ell} is not prime")
    require(target > 0, "target rank must be positive")
    if candidates:
        loop_head = f"C=[{','.join(map(str, candidates))}]; for(ii=1,#C,{{D=C[ii];"
    else:
        require(args.start is not None and args.stop is not None,
                "search needs --candidate(s) or both --start and --stop")
        require(args.step > 0, "step must be positive")
        loop_head = f"forstep(D={args.start},{args.stop},{args.step},{{"
    # q.gen is an invariant-factor basis.  Raising its i-th element by
    # cyc[i]/ell produces an independent element of exact order ell.
    return f"""
default(parisizemax,2000000000); allocatemem(64000000);
ell={ell}; target={target}; checked=0; fundamental=0; best=-1; hist=vector(32,i,0);
{loop_head}
  checked++;
  if(D>0 && isfundamental(-D),
    fundamental++;
    q=quadclassunit(-D);
    rk=sum(i=1,#q.cyc,valuation(q.cyc[i],ell)>0);
    if(rk+1>#hist,hist=concat(hist,vector(rk+1-#hist)));
    hist[rk+1]++;
    if(rk>best,best=rk;print("IMPROVE|D=",D,"|h=",q.no,"|cyc=",q.cyc,"|ell_rank=",rk));
    if(rk>=target,
      line=Str(ell,"|",target,"|",D); taken=0;
      for(i=1,#q.cyc,if(valuation(q.cyc[i],ell)>0 && taken<target,
        g=q.gen[i]^(q.cyc[i]/ell); v=Vec(g);
        line=Str(line,"|",v[1],",",v[2],",",v[3]); taken++));
      print("CERT|",line)
    )
  )
}});
print("SEARCH_SUMMARY|checked=",checked,"|fundamental=",fundamental,"|best=",best,"|hist=",hist);
quit
"""


def run_search(args: argparse.Namespace) -> int:
    candidates = candidate_list(args)
    program = search_program(args, candidates)
    process = subprocess.run(
        [args.pari_gp, "-fq"],
        input=program,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        check=False,
    )
    sys.stdout.write(process.stdout)
    require(process.returncode == 0,
            f"PARI/GP search exited with status {process.returncode}")
    require_gp_success(process.stdout, "SEARCH_SUMMARY|", "search")
    certificates = [
        row.split("CERT|", 1)[1]
        for row in process.stdout.splitlines()
        if row.startswith("CERT|")
    ]
    for certificate in certificates:
        verification = verify_certificate(certificate)
        print(
            "EXTRACTED_CERTIFICATE_PURE_VERIFY=PASS"
            f"|D={verification.D}|r={verification.rank}"
        )
    return 0


def reduced_forms(delta: int) -> List[Form]:
    """Enumerate all primitive reduced positive-definite forms of delta."""
    require(delta < 0, "form enumeration needs negative discriminant")
    forms: List[Form] = []
    # For a reduced form, 3a^2 <= |delta|.
    for a in range(1, isqrt((-delta) // 3) + 2):
        for b in range(-a, a + 1):
            numerator = b * b - delta
            if numerator % (4 * a):
                continue
            c = numerator // (4 * a)
            if a > c:
                continue
            if b < 0 and (abs(b) == a or a == c):
                continue
            form = Form(a, b, c)
            if gcd(a, gcd(abs(b), c)) == 1:
                forms.append(form)
    return forms


SAMPLE = (
    "3|3|4447704|390,-96,2857|921,-786,1375|346,-68,3217"
)


def run_selftest(args: argparse.Namespace) -> int:
    # These include odd/even discriminants, nontrivial 2-torsion, and class
    # groups of orders 1 through 8.  Exhaustive associativity is independent
    # of PARI and exercises the non-coprime composition branch repeatedly.
    deltas = (-3, -4, -7, -8, -15, -20, -23, -24, -31, -39, -40, -47, -84)
    composition_gates = 0
    associativity_gates = 0
    for delta in deltas:
        forms = reduced_forms(delta)
        form_set = set(forms)
        identity = principal_form(delta)
        require(identity in form_set,
                f"principal form missing from enumeration at {delta}")
        for left in forms:
            require(left.reduced() == left, f"noncanonical enumeration at {delta}")
            require(left.compose(identity) == left,
                    f"right identity failed at {delta}")
            require(identity.compose(left) == left,
                    f"left identity failed at {delta}")
            require(left.compose(left.inverse()) == identity,
                    f"inverse failed at {delta}")
            for right in forms:
                product = left.compose(right)
                require(product in form_set,
                        f"closure failed at {delta}: {left}*{right}={product}")
                require(product == right.compose(left),
                        f"commutativity failed at {delta}")
                composition_gates += 1
                for third in forms:
                    require(
                        product.compose(third)
                        == left.compose(right.compose(third)),
                        f"associativity failed at {delta}",
                    )
                    associativity_gates += 1

    sample = verify_certificate(SAMPLE)
    require(sample.left_size == 3 and sample.right_size == 9,
            "sample meet-in-the-middle sizes changed")
    duplicate = (
        "3|3|4447704|390,-96,2857|390,-96,2857|346,-68,3217"
    )
    try:
        verify_certificate(duplicate)
    except ValueError:
        hostile = "PASS"
    else:
        raise ValueError("dependent duplicate-generator hostile was accepted")

    print(
        "PURE_SELFTEST=PASS"
        f"|discriminants={len(deltas)}"
        f"|composition_gates={composition_gates}"
        f"|associativity_gates={associativity_gates}"
        f"|duplicate_hostile={hostile}"
        "|sample=PASS"
    )
    if args.pari_gp:
        checks: List[str] = []
        pari_pairs = 0
        for delta in deltas:
            forms = reduced_forms(delta)
            for left in forms:
                for right in forms:
                    product = left.compose(right)
                    checks.append(
                        "if(Vec(Qfb(%d,%d,%d)*Qfb(%d,%d,%d))"
                        "!=[%d,%d,%d],error(\"composition mismatch\"));"
                        % (
                            left.a, left.b, left.c,
                            right.a, right.b, right.c,
                            product.a, product.b, product.c,
                        )
                    )
                    pari_pairs += 1
        program = "\n".join(checks) + (
            f'\nprint("PARI_COMPOSITION_AUDIT=PASS|pairs={pari_pairs}");\nquit\n'
        )
        process = subprocess.run(
            [args.pari_gp, "-fq"],
            input=program,
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            check=False,
        )
        require(process.returncode == 0,
                f"PARI composition audit exited {process.returncode}")
        require_gp_success(
            process.stdout, "PARI_COMPOSITION_AUDIT=PASS|", "composition audit"
        )
        print(process.stdout.strip())
        print(pari_audit(SAMPLE, args.pari_gp))
    return 0


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    verify = subparsers.add_parser("verify", help="verify certificate lines")
    verify.add_argument("line", nargs="*", help="literal certificate line(s)")
    verify.add_argument("--file", help="file containing one certificate per line")
    verify.add_argument("--pari-gp", help="optional independent PARI/GP executable")
    verify.set_defaults(func=run_verify)

    search = subparsers.add_parser("search", help="search discriminants with PARI/GP")
    search.add_argument("--pari-gp", required=True, help="path to gp executable")
    search.add_argument("--ell", type=int, required=True)
    search.add_argument("--target-rank", type=int, required=True)
    search.add_argument("--candidate", type=int, action="append",
                        help="individual positive D; may be repeated")
    search.add_argument("--candidates-file", help="one positive D per line")
    search.add_argument("--start", type=int)
    search.add_argument("--stop", type=int)
    search.add_argument("--step", type=int, default=1)
    search.set_defaults(func=run_search)

    selftest = subparsers.add_parser(
        "selftest", help="exhaustive small-class arithmetic and hostile tests"
    )
    selftest.add_argument("--pari-gp", help="optional independent PARI/GP executable")
    selftest.set_defaults(func=run_selftest)
    return parser


def main(argv: Optional[Sequence[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    try:
        return args.func(args)
    except (ArithmeticError, RuntimeError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
