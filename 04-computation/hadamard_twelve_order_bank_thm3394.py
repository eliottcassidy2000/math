#!/usr/bin/env python3
"""Inert exact renderer/verifier for the THM-3394 Hadamard bank.

The only external input is the checked-in base85/zlib sign word named below.
This program has no shell, subprocess, eval, network, or dynamic-code path.
It implements a small fixed text-rewrite language because that is the most
literal independent reconstruction of the supplied construction templates.
"""

from __future__ import annotations

import base64
from dataclasses import dataclass
import hashlib
from pathlib import Path
import re
import zlib


ROOT = Path(__file__).resolve().parents[1]
SIGNWORD_PATH = ROOT / "04-computation" / "hadamard_twelve_order_signword_thm3394.b85"
SIGNWORD_SHA256 = "5b5fe8fa42f0d6a8b4e4c9926726d82a6aab8e1070c1ae4d1b430c1277e58db4"
SIGNWORD_ZLIB_SHA256 = "1756297611d2bb403e9c4152ea91146428482983716aad165a3cc21396d5a61c"
SCHEDULE_BITS = (
    "111000101001100011100010110010001100001101111100"
    "110001000110110011000100110111001000010101011000"
    "011001011000000011100110100010001110011011101000"
    "101001110111000011000111100111001100011110101100"
)
SCHEDULE_SHA256 = "c18071311a78311819a20a1ba7823de6581bb907fd9b4b44a3a9b17860a370ff"
KIND_FROM_BITS = {"111": "h", "110": "g", "101": "w", "100": "v", "011": "u"}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise RuntimeError(message)


@dataclass
class Address:
    kind: str
    value: object
    negate: bool = False


@dataclass
class Command:
    address: Address | None
    op: str
    args: tuple = ()


def read_delimited(source: str, start: int, delimiter: str) -> tuple[str, int]:
    out: list[str] = []
    i = start
    while i < len(source):
        char = source[i]
        if char == "\\" and i + 1 < len(source):
            following = source[i + 1]
            if following == "\n":
                out.append("\n")
                i += 2
                continue
            out.extend((char, following))
            i += 2
            continue
        if char == delimiter:
            return "".join(out), i + 1
        out.append(char)
        i += 1
    raise ValueError("unterminated fixed rewrite field")


def parse_rewrite_program(source: str) -> tuple[list[Command], dict[str, int]]:
    commands: list[Command] = []
    labels: dict[str, int] = {}
    i = 0
    while i < len(source):
        while i < len(source):
            if source[i] in " \t\r\n;":
                i += 1
            elif source.startswith("\\\n", i):
                i += 2
            else:
                break
        if i >= len(source):
            break

        address = None
        if source[i].isdigit():
            j = i
            while j < len(source) and source[j].isdigit():
                j += 1
            first = int(source[i:j])
            if j < len(source) and source[j] == ",":
                k = j + 1
                while k < len(source) and source[k].isdigit():
                    k += 1
                address = Address("range", (first, int(source[j + 1:k])))
                i = k
            else:
                address = Address("line", first)
                i = j
        elif source[i] == "/":
            pattern, i = read_delimited(source, i + 1, "/")
            address = Address("regex", pattern)
        if address is not None and i < len(source) and source[i] == "!":
            address.negate = True
            i += 1

        require(i < len(source), "missing fixed rewrite operation")
        op = source[i]
        i += 1
        if op == ":":
            j = i
            while j < len(source) and source[j] not in ";\n":
                j += 1
            name = source[i:j].strip()
            labels[name] = len(commands)
            commands.append(Command(address, op, (name,)))
            i = j
        elif op in "tb":
            j = i
            while j < len(source) and source[j] not in ";\n":
                j += 1
            commands.append(Command(address, op, (source[i:j].strip(),)))
            i = j
        elif op == "s":
            delimiter = source[i]
            pattern, i = read_delimited(source, i + 1, delimiter)
            replacement, i = read_delimited(source, i, delimiter)
            j = i
            while j < len(source) and source[j] not in ";\n":
                if source.startswith("\\\n", j):
                    j += 2
                else:
                    j += 1
            flags = source[i:j].strip().replace("\\", "")
            commands.append(Command(address, op, (pattern, replacement, flags)))
            i = j
        elif op == "y":
            delimiter = source[i]
            left, i = read_delimited(source, i + 1, delimiter)
            right, i = read_delimited(source, i, delimiter)
            commands.append(Command(address, op, (left, right)))
        elif op in "hHgGxPpd":
            commands.append(Command(address, op))
        else:
            raise ValueError(f"operation {op!r} is outside the fixed rewrite language")
    return commands, labels


def basic_regex_to_python(pattern: str) -> str:
    out: list[str] = []
    i = 0
    in_class = False
    while i < len(pattern):
        char = pattern[i]
        if char == "[" and not in_class:
            in_class = True
            out.append(char)
            i += 1
            continue
        if char == "]" and in_class:
            in_class = False
            out.append(char)
            i += 1
            continue
        if char == "\\" and i + 1 < len(pattern):
            following = pattern[i + 1]
            if following in "(){}+?|":
                out.append(following)
            elif following == "n":
                out.append("\n")
            elif following.isdigit():
                out.append("\\" + following)
            elif following == "\\":
                out.append(r"\\")
            else:
                out.append("\\" + following)
            i += 2
            continue
        if not in_class and char in "+?{}|()":
            out.append("\\" + char)
        else:
            out.append(char)
        i += 1
    return "".join(out)


def expand_replacement(replacement: str, match: re.Match[str]) -> str:
    out: list[str] = []
    i = 0
    while i < len(replacement):
        char = replacement[i]
        if char == "&":
            out.append(match.group(0))
            i += 1
            continue
        if char == "\\" and i + 1 < len(replacement):
            following = replacement[i + 1]
            if following.isdigit():
                out.append(match.group(int(following)) or "")
            elif following == "n":
                out.append("\n")
            elif following == "&":
                out.append("&")
            elif following == "\\":
                out.append("\\")
            else:
                out.append(following)
            i += 2
            continue
        out.append(char)
        i += 1
    return "".join(out)


def rewrite_sub(pattern: str, replacement: str, text: str, global_flag: bool) -> tuple[str, bool]:
    expression = re.compile(basic_regex_to_python(pattern), re.DOTALL)
    changed = False

    def replace(match: re.Match[str]) -> str:
        nonlocal changed
        changed = True
        return expand_replacement(replacement, match)

    return expression.sub(replace, text, count=0 if global_flag else 1), changed


def address_matches(address: Address | None, pattern: str, line_number: int) -> bool:
    if address is None:
        return True
    if address.kind == "line":
        result = line_number == address.value
    elif address.kind == "range":
        first, last = address.value
        result = first <= line_number <= last
    else:
        result = re.search(basic_regex_to_python(str(address.value)), pattern, re.DOTALL) is not None
    return not result if address.negate else result


def run_rewrite(source: str, input_lines: list[str], auto_print: bool = True) -> str:
    commands, labels = parse_rewrite_program(source)
    output: list[str] = []
    hold = ""
    steps = 0
    for line_number, initial in enumerate(input_lines, 1):
        pattern = initial
        substitution = False
        program_counter = 0
        deleted = False
        while program_counter < len(commands):
            steps += 1
            require(steps <= 100_000_000, "fixed rewrite exceeded its step bound")
            command = commands[program_counter]
            if not address_matches(command.address, pattern, line_number):
                program_counter += 1
                continue
            op = command.op
            if op == ":":
                program_counter += 1
            elif op == "h":
                hold = pattern
                program_counter += 1
            elif op == "H":
                hold += "\n" + pattern
                program_counter += 1
            elif op == "g":
                pattern = hold
                program_counter += 1
            elif op == "G":
                pattern += "\n" + hold
                program_counter += 1
            elif op == "x":
                pattern, hold = hold, pattern
                program_counter += 1
            elif op == "y":
                left, right = command.args
                pattern = pattern.translate(str.maketrans(left, right))
                program_counter += 1
            elif op == "s":
                old, new, flags = command.args
                pattern, changed = rewrite_sub(old, new, pattern, "g" in flags)
                substitution = substitution or changed
                if changed and "p" in flags:
                    output.append(pattern + "\n")
                program_counter += 1
            elif op == "P":
                output.append(pattern.split("\n", 1)[0] + "\n")
                program_counter += 1
            elif op == "p":
                output.append(pattern + "\n")
                program_counter += 1
            elif op == "d":
                deleted = True
                break
            elif op == "t":
                if substitution:
                    substitution = False
                    program_counter = labels[command.args[0]]
                else:
                    program_counter += 1
            elif op == "b":
                program_counter = labels[command.args[0]]
        if auto_print and not deleted:
            output.append(pattern + "\n")
    return "".join(output)


H_TEMPLATE = "\n".join(
    (
        r"h;y/+/-/;H;y/-/+/;G;s/^=.*\n\(=AA\).*/-++-\2\1\\",
        r"+-+-\3\3\1\3\\",
        r"++--\3\1\3\3\\",
        r"----\1\2/p;g;s/\n.*//",
        r"s/^=\(.*\)/\1\\",
        r"\2\\",
        r"\2/",
        r":4",
        r"s/\n\(.*\)\(.\)\n/\2\\",
        r"\1\\",
        r"/",
        r"t4",
        r"h;s/^====.*/+++-\1\4\3\2/L0rlll",
        r"s/^AAA=ZA.=AAAAZ.*/--+-\1\4\6\5\3\2/L1lrll",
        r"s/^AA=AAZ.=AAZ.*/-+--\1\3\2\4\6\5/L2llrl",
        r"s/^A=AZAA.=AAAZ.*/+---\1\6\5\3\2\4/L3lllr",
        r"d",
        "",
    )
)
# Raw one-line literals cannot end in a lone backslash.  The tuple spells
# those continuations doubled; collapse only a doubled slash before newline.
H_TEMPLATE = H_TEMPLATE.replace("\\\\\n", "\\\n")


def expand_l_macro(line: str) -> str:
    match = re.search(r"L(.)((?:.|\n)*)", line)
    if not match:
        return line
    label, remainder = match.group(1), match.group(2)
    expansion = (
        ";s/.*/&\\\n&/\n"
        + ":" + label + "\n"
        + "P\n"
        + "s/" + remainder + "\\n/\\2\\1\\4\\3\\6\\5\\8\\7\\\n"
        + "/\n"
        + "/^\\(.*\\)\\n\\1/!b" + label + "\n"
        + "g;y/+-/-+/;G"
    )
    return line[:match.start()] + expansion


def make_h_program(length: int, unbordered: bool) -> str:
    lines = H_TEMPLATE.splitlines()
    if unbordered:
        lines = lines[4:]
        lines = [re.sub(r"/[+-]{4}", "/", line, count=1) for line in lines]
    left = length // 8
    right = length // 4 - 1 - left
    percent = rf".\{{{left}\}}.\{{{right}\}}"
    output: list[str] = []
    for line in lines:
        line = line.replace("r", r"\(%\)\(.\)", 1)
        line = line.replace("l", r"\(.\)\(%\)")
        line = expand_l_macro(line)
        line = line.replace("=", r"\(A\)")
        line = line.replace("Z", r"%...\(.\)\(%\)")
        line = line.replace("A", r"%.")
        line = line.replace("%", percent)
        output.append(line)
    return "\n".join(output) + "\n"


def render_h(data: str, unbordered: bool = False) -> list[str]:
    return run_rewrite(make_h_program(len(data), unbordered), [data]).splitlines()


META_W = r''':c
/B............[-+]/s/+\([^-+]*\)$/\\2\1/;/B............[-+]/s/-\([^-+]*\)$/\\1\1/;tc
s/B//;s|L\(.\)\(.\)\(.\)\(.\)\(.\)\(.\)|\
:\1\
/^[-+]*\\nZ\\nZ\\nZ$/s/==/\\2\\1/g;/^[-+]*\\nZ\\nZ$/s/====/\\4\\3\\2\\1/g;/^[-+]*\\nZ$/s/==/\\2\\1/g;s/^\\([-+]*\\)\\n\\(Z\\)/\\2\\1\\\
\\2\\1/\
:\2\
P;s/^\\(Z\\)\3\3\3\3/~/;s/^\\(ZK\\)\4\4\4\4/~/;s/^\\(ZKK\\)\5\5\5\5/~/;s/^\\(ZKKK\\)\6\6\6\6/~/;/^\\(.*\\)\\n\\1/!b\2\
s/^[-+]*\\nZ//;/\\n/b\1\
g;y/+-/-+/;G|;s/J/KKKK\\n/g;s/I/\\(K\\)/g;s/r/\\(%\\)\\(.\\)/g;s/l/\\(.\\)\\(%\\)/g;s/~/\\1\\3\\2\\5\\4\\7\\6\\9\\8/g;s/Z/............/g;s/K/AAAA/g;s/=/\\(A\\)/g;s/A/%./g;s/%/.\\{118\\}/g'''

META_V = r''':c
/S[-+]/s/+\([^-+]*\)$/\\2\1/;/S[-+]/s/-\([^-+]*\)$/\\1\1/;tc
s|X\(.\)\(.\)\(.\)\(.\)\(.\)|\
:\1\
s/^\\([-+]*\\)\\n\\(Z\\)/\\2\\1\\\
\\2\\1/\
:\1\1\
P;s/^\\(Z\\)\2\2\2/~/;s/^\\(ZAAA\\)\2\2\2/~/;s/^\\(ZK\\)\3\3\3/~/;s/^\\(ZKAAA\\)\3\3\3/~/;s/^\\(ZKK\\)\4\4\4/~/;s/^\\(ZKKAAA\\)\4\4\4/~/;s/^\\(ZKKK\\)\5\5\5/~/;s/^\\(ZKKKAAA\\)\5\5\5/~/;/^\\(.*\\)\\n\\1/!b\1\1|;s|Y\(.\)\(K*\)|\
s/^[-+]*\\n//;s/^\\(Z\2\\)QQ/\\1\\3\\2\\5\\4/;s/^Z//;s/Q/\\2\\1/g;/^[-+]*\\nZ\\nZ\\nZ$/s/WW/\\2\\1/g;/\\n/b\1\
g;y/+-/-+/;G|;s/J/KKKK\\n/g;s/I/\\(K\\)/g;s/W/\\(AAA\\)/g;s/Q/=\\(AA\\)/g;s/r/\\(%\\)\\(.\\)/g;s/l/\\(.\\)\\(%\\)/g;s/~/\\1\\3\\2\\5\\4\\7\\6/g;s/Z/..................../g;s/K/AAAAAA/g;s/=/\\(A\\)/g;s/A/%./g;s/%/.\\{56\\}/g'''

META_U = r''':c
/S[-+]/s/+\([^-+]*\)$/\\2\1/;/S[-+]/s/-\([^-+]*\)$/\\1\1/;tc
/S/s/\\.\\./&&&&/g;s|X\(.\)\(.\)\(.\)\(.\)\(.\)|\
:\1\
s/^\\([-+]*\\)\\n\\(Z\\)/\\2\\1\\\
\\2\\1/\
:\1\1\
P;s/^\\(Z\\)\2\2\2\2/~/;s/^\\(ZAAAA\\)\2\2\2\2/~/;s/^\\(ZC\\)\2\2\2\2/~/;s/^\\(ZCAAAA\\)\2\2\2\2/~/\
s/^\\(ZCC\\)\2\2\2\2/~/;s/^\\(ZCCAAAA\\)\2\2\2\2/~/;s/^\\(ZCCC\\)\2\2\2\2/~/;s/^\\(ZCCCAAAA\\)\2\2\2\2/~/\
s/^\\(ZK\\)\3\3\3\3/~/;s/^\\(ZKAAAA\\)\3\3\3\3/~/;s/^\\(ZKC\\)\3\3\3\3/~/;s/^\\(ZKCAAAA\\)\3\3\3\3/~/\
s/^\\(ZKCC\\)\3\3\3\3/~/;s/^\\(ZKCCAAAA\\)\3\3\3\3/~/;s/^\\(ZKCCC\\)\3\3\3\3/~/;s/^\\(ZKCCCAAAA\\)\3\3\3\3/~/\
s/^\\(ZKK\\)\4\4\4\4/~/;s/^\\(ZKKAAAA\\)\4\4\4\4/~/;s/^\\(ZKKC\\)\4\4\4\4/~/;s/^\\(ZKKCAAAA\\)\4\4\4\4/~/\
s/^\\(ZKKCC\\)\4\4\4\4/~/;s/^\\(ZKKCCAAAA\\)\4\4\4\4/~/;s/^\\(ZKKCCC\\)\4\4\4\4/~/;s/^\\(ZKKCCCAAAA\\)\4\4\4\4/~/\
s/^\\(ZKKK\\)\5\5\5\5/~/;s/^\\(ZKKKAAAA\\)\5\5\5\5/~/;s/^\\(ZKKKC\\)\5\5\5\5/~/;s/^\\(ZKKKCAAAA\\)\5\5\5\5/~/\
s/^\\(ZKKKCC\\)\5\5\5\5/~/;s/^\\(ZKKKCCAAAA\\)\5\5\5\5/~/;s/^\\(ZKKKCCC\\)\5\5\5\5/~/;s/^\\(ZKKKCCCAAAA\\)\5\5\5\5/~/;/^\\(.*\\)\\n\\1/!b\1\1|;s|Y\(.\)\(K*\)|\
s/^[-+]*\\n//;s/^\\(Z\2\\)QQQQ/~/;s/^Z//\
s/^WWWW/M/;s/^IWWWW/~/;s/^\\(KK\\)WWWW/~/;s/^\\(KKK\\)WWWW/~/\
/^[-+]*OOO$/s/^VVVVVVVV/M/\
/^[-+]*OOO$/s/^\\(KK\\)VVVVVVVV/~/\
/^[-+]*OO$/s/^VVVVVVVV/M/\
/^[-+]*OO$/s/^\\(KK\\)VVVVVVVV/~/\
/^[-+]*OO$/s/^UUUUUUUU/M/\
/^[-+]*O$/s/^VVVVVVVV/M/\
/^[-+]*O$/s/^\\(KK\\)VVVVVVVV/~/;/\\n/b\1\
g;y/+-/-+/;G|;s/J/KKKK\\n/g;s/I/\\(K\\)/g;s/V/\\(C\\)/g;s/U/\\(CC\\)/g;s/W/\\(A\\)\\(.\\{77\\}\\)/g;s/Q/\\(.\\{66\\}\\)\\(.\\{22\\}\\)/g;s/O/\\nZ\\nZ\\nZ\\nZ\\nZ\\nZ\\nZ\\nZ/g;s/M/\\2\\1\\4\\3\\6\\5\\8\\7/g;s/~/\\1\\3\\2\\5\\4\\7\\6\\9\\8/g;s/r/\\(%\\)\\(.\\)/g;s/l/\\(.\\)\\(%\\)/g;s/Z/............................/g;s/K/.\\{176\\}.\\{176\\}/g;s/C/.\\{88\\}/g;s/=/\\(A\\)/g;s/A/%./g;s/%/.\\{10\\}/g'''

BODIES_W = (
    r'''/p;g;s/\n.*//;h;s/l/\1,\2;/g
:8
s/,\([-+]*\)\([-+]\);/\2,\1;/g;t8
s/,;//g;G;h;y/+-/-+/;G;s/^JJK\(KKK\)\nI.*/\2\1\
''',
    r'''/L01rlll
s/^KIKK\nKKIK\nJIKKI/\1\3\4\2\
''',
    r'''/L23lrll
s/^KKIK\nKKKI\nJ\(KK\).*/\1\2\3\
''',
    r'''/L45llrl
s/^KKKI\nKIKK\nJIKI.*/\1\4\2\3\
''',
    '''/L67lllr
d
''',
)

BODIES_V = (
    r'''/;s/S//g;p;g;s/.*\n//;h;s/l/\1,\2;/g
:8
s/,\([-+]*\)\([-+]\);/\2,\1;/g;t8
s/,;//g;s/===/\1\3\2/g;G;h;y/+-/-+/;G;s/^JJK\(KKK\)\nI.*/\2\1\
''',
    r'''/X0rlllY0
s/^KIKK\nKKIK\nJIKKI/\1\3\4\2\
''',
    r'''/X1lrllY1K
s/^KKIK\nKKKI\nJ\(KK\).*/\1\2\3\
''',
    r'''/X2llrlY2KK
s/^KKKI\nKIKK\nJIKI.*/\1\4\2\3\
''',
    '''/X3lllrY3KKK
d
''',
)

BODIES_U = (
    r'''/;s/S//g;p;g;s/.*\n//;h;s/l/\1,\2;/g
:8
s/,\([-+]*\)\([-+]\);/\2,\1;/g;t8
s/,;//g;s/========/\1\8\7\6\5\4\3\2/g;G;h;y/+-/-+/;G;s/^JJK\(KKK\)\nI.*/\2\1\
''',
    r'''/X0rlllY0
s/^KIKK\nKKIK\nJIKKI/\1\3\4\2\
''',
    r'''/X1lrllY1K
s/^KKIK\nKKKI\nJ\(KK\).*/\1\2\3\
''',
    r'''/X2llrlY2KK
s/^KKKI\nKIKK\nJIKI.*/\1\4\2\3\
''',
    '''/X3lllrY3KKK
d
''',
)

TOP_SUFFIX = {
    "w": r"h;y/-/+/;H;y/+/-/;G;s/^=[-+]*\n[-+]*\n=.*/",
    "v": r"h;y/-/+/;x;H;y/+/-/;G;s/^=[-+]*\n=.*/",
    "u": r"h;y/-/+/;x;H;y/+/-/;G;s/^=[-+]*\n=.*/",
}

COMPOSITE_SPEC = {
    "w": (336, 28, "B", 12, (48, 48, 48, 48), META_W, BODIES_W),
    "v": (880, 44, 20, 20, (120, 120, 120, 120), META_V, BODIES_V),
    "u": (1680, 60, 28, 28, (896, 896, 896, 896), META_U, BODIES_U),
}


def folded_text(chunk: str, width: int, insertion: str | int | None, suffix: str) -> str:
    rows = [chunk[i:i + width] for i in range(0, len(chunk), width)]
    require(rows and all(len(row) == width for row in rows), "non-integral fixed fold")
    if insertion == "B":
        rows = ["B" + row for row in rows]
    elif isinstance(insertion, int):
        rows = [row[:insertion] + "S" + row[insertion:] for row in rows]
    else:
        require(insertion is None, "unknown fixed insertion")
    return suffix + "\\\n".join(rows)


def render_composite(kind: str, primary: str, pieces: list[str]) -> list[str]:
    side_total, top_width, insertion, side_width, side_lengths, meta, bodies = COMPOSITE_SPEC[kind]
    require(len(pieces) == 5, "composite construction requires five side pieces")
    require(len(pieces[0]) == side_total, "wrong composite top-side length")
    require(tuple(map(len, pieces[1:])) == side_lengths, "wrong composite OA-side lengths")
    q_source = folded_text(pieces[0], top_width, insertion, TOP_SUFFIX[kind]) + bodies[0]
    for piece, body in zip(pieces[1:], bodies[1:]):
        q_source += folded_text(piece, side_width, None, "") + body
    generated_program = run_rewrite(meta, q_source.splitlines(), auto_print=True)
    return run_rewrite(generated_program, [primary], auto_print=True).splitlines()


EXPECTED = {
    668: (
        "bdeb5059d77e2703211082627b60441b8c888c928a55cc6f295e011941a387b0",
        "73f1de1539849e1dc7e6085cc69c563fd2965c44970263e8203384bd1a46aa63",
    ),
    716: (
        "3adcb1bb2884467d9e34069a3b32950728adabcdb8b35a4503d20c3312664ee6",
        "79ebae74dbadae11059aeabf77b61144ccf6fd9905c2bf8167ac98854651366c",
    ),
    892: (
        "e77fc79ab287f5f5ba5bbdc10191bdc7593839052fe1015c1fb6a2e974ab54de",
        "73d28c9e58d4c3bac1e2bc0cc0d6f0b504eb7719fa2370e9403692e4dfba7f6c",
    ),
    1132: (
        "7d1c1e892149e90330d58bb0cf9ef2c888078df1b35fb55f8724d580ebf7b743",
        "f0499c89d761d57d330471a1446a0f099723f663b5602285aa9dcc1070b0216b",
    ),
    1244: (
        "4cb747cf511eba1f203582b5121bdf6ab02671133e45579c1d023add8b2da143",
        "cffd1c77f09646b9b17e4c46ed9f6bf091bd56b2c2cc06c83c64054b44d74d33",
    ),
    1388: (
        "a6b92584eb803b87026709d64fe892dec8f7182a120e13de9edd3065cf05bf0b",
        "5ebac7a33165a39cb9e86e5b99eb42cc46e3c7006956b2c97885f8f6b4655db1",
    ),
    1436: (
        "e4d745a4d44f39a5671f9cd86f5c1d0aef93504dcfb2e253451cadc9e4086728",
        "99984bf728bcf37a360ba0af57076a729653070b026271085f96dce35703a00e",
    ),
    1676: (
        "8e919c2bdb4d30c34817eb5650d2dd3d82d7c6504feccd96c5ca22a2191cdb99",
        "787f9ec132cfa6900090e78bbebd8c1132af583793b4c65494ca52f7000a1e68",
    ),
    1772: (
        "1852e951db69c44eb95b37ed741c3ff2e29691267eaf872d6a9da3a977236ba2",
        "020f15819f78aab232da7ec9183caf172cc6ab1cdad191f4fabb9c63ec308af1",
    ),
    1916: (
        "be2073eeaa5399cfe104023829d2c6770b49dd2f07bf6347203f1cbd75577ae9",
        "081197bb0c94e9e1aa628fa3170db50d7df45e2889494bf12ee5abcdfd4266da",
    ),
    1948: (
        "fddc841ebf951f6e17e939551d058ea5df046251ea065b5f6e7ee2fd8d0f62ce",
        "ee885b2ade31910d36fd7bc11321388c20aa64eaef081fecacf05f3b53d5d792",
    ),
    1964: (
        "740b907cd442f1b7fd40dcc31f2b3aae9794842da6dc579f98dac1d0d9e1493d",
        "4a3e301084d016f6b11e32d425b9f8f54e79cd2b52ec8c6e64420635c13aaa53",
    ),
}


def sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def load_signword() -> str:
    encoded = b"".join(SIGNWORD_PATH.read_bytes().split())
    require(
        sha256(encoded) == "7831fb93c91b82e70827345c2b72eee6702a30ffce3f5162b01107e5ce35fbc0",
        "base85 signword hash mismatch",
    )
    compressed = base64.b85decode(encoded)
    require(len(compressed) == 3488, "compressed signword length mismatch")
    require(sha256(compressed) == SIGNWORD_ZLIB_SHA256, "compressed signword hash mismatch")
    signword = zlib.decompress(compressed).decode("ascii")
    require(len(signword) == 23828, "signword length mismatch")
    require(set(signword) == {"+", "-"}, "signword alphabet mismatch")
    require(signword.count("+") == 12188 and signword.count("-") == 11640, "signword count mismatch")
    require(sha256(signword.encode()) == SIGNWORD_SHA256, "signword hash mismatch")
    return signword


def decode_schedule() -> list[tuple[str, int]]:
    require(len(SCHEDULE_BITS) == 192, "schedule is not 192 bits")
    require(set(SCHEDULE_BITS) == {"0", "1"}, "schedule is not binary")
    require(sha256(SCHEDULE_BITS.encode()) == SCHEDULE_SHA256, "schedule hash mismatch")
    records: list[tuple[str, int]] = []
    for start in range(0, len(SCHEDULE_BITS), 16):
        record = SCHEDULE_BITS[start:start + 16]
        require(record[:3] in KIND_FROM_BITS, "unknown schedule opcode")
        records.append((KIND_FROM_BITS[record[:3]], int(record[3:], 2)))
    require(
        records
        == [
            ("h", 664), ("h", 712), ("g", 892), ("g", 1132),
            ("g", 1244), ("v", 1368), ("u", 1408), ("h", 1672),
            ("h", 1768), ("w", 1904), ("g", 1948), ("g", 1964),
        ],
        "decoded schedule mismatch",
    )
    return records


def sign_values(row: str) -> list[int]:
    return [1 if char == "+" else -1 for char in row]


def periodic_contract(data: str, kind: str) -> tuple[list[int], int]:
    require(kind in ("h", "g") and len(data) % 4 == 0, "bad periodic contract input")
    width = len(data) // 4
    sequences = [sign_values(data[i * width:(i + 1) * width]) for i in range(4)]
    expected = -4 if kind == "h" else 0
    for shift in range(1, width):
        total = sum(
            sum(row[i] * row[(i + shift) % width] for i in range(width))
            for row in sequences
        )
        require(total == expected, f"{kind} PAF failure at shift {shift}: {total}")
    row_sums = [sum(row) for row in sequences]
    expected_square_sum = 4 if kind == "h" else len(data)
    require(sum(value * value for value in row_sums) == expected_square_sum, "row-sum square identity failed")
    return row_sums, expected


def oa_contract(rows: list[str]) -> tuple[list[int], int]:
    require(len(rows) == 4 and len({len(row) for row in rows}) == 1, "bad OA sidecar shape")
    values = [sign_values(row) for row in rows]
    length = len(rows[0])
    for row in values:
        require(sum(row) == 0, "OA sidecar has a nonzero first moment")
    for left in range(4):
        for right in range(left):
            require(
                sum(values[left][i] * values[right][i] for i in range(length)) == 0,
                "OA sidecar has a nonzero second moment",
            )
    degree_three = []
    for indices in ((0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)):
        degree_three.append(
            sum(values[indices[0]][j] * values[indices[1]][j] * values[indices[2]][j]
                for j in range(length))
        )
    degree_four = sum(
        values[0][j] * values[1][j] * values[2][j] * values[3][j]
        for j in range(length)
    )
    return degree_three, degree_four


def normalized_rows(rows: list[str]) -> list[str]:
    first = rows[0]
    top_left = first[0] == "-"
    output = []
    for row in rows:
        first_in_row = row[0] == "-"
        output.append(
            "".join(
                "-" if ((char == "-") ^ (first[column] == "-") ^ first_in_row ^ top_left) else "+"
                for column, char in enumerate(row)
            )
        )
    require(set(output[0]) == {"+"}, "normalization failed on first row")
    require(all(row[0] == "+" for row in output), "normalization failed on first column")
    return output


def verify_matrix(rows: list[str], order: int) -> tuple[str, str, int, int, int]:
    require(len(rows) == order, f"row count is not {order}")
    require(all(len(row) == order for row in rows), f"column count is not {order}")
    require(all(set(row) <= {"+", "-"} for row in rows), "matrix has a non-sign entry")
    integers = [int(row.translate(str.maketrans("+-", "10")), 2) for row in rows]
    minimum = order
    maximum = 0
    pairs = 0
    for i, left in enumerate(integers):
        for right in integers[:i]:
            distance = (left ^ right).bit_count()
            minimum = min(minimum, distance)
            maximum = max(maximum, distance)
            pairs += 1
            require(distance == order // 2, f"orthogonality failure in order {order}")
    require(pairs == order * (order - 1) // 2, "pair count mismatch")
    raw_hash = sha256(("\n".join(rows) + "\n").encode())
    normal_hash = sha256(("\n".join(normalized_rows(rows)) + "\n").encode())
    require((raw_hash, normal_hash) == EXPECTED[order], f"matrix hash mismatch at order {order}")
    return raw_hash, normal_hash, pairs, minimum, maximum


def single_entry_hostile(rows: list[str]) -> tuple[int, int, int]:
    altered = rows.copy()
    altered[0] = ("-" if rows[0][0] == "+" else "+") + rows[0][1:]
    left = int(altered[0].translate(str.maketrans("+-", "10")), 2)
    for index, row in enumerate(altered[1:], 1):
        right = int(row.translate(str.maketrans("+-", "10")), 2)
        distance = (left ^ right).bit_count()
        if distance != len(rows) // 2:
            return 0, index, distance
    raise RuntimeError("single-entry hostile was not rejected")


def main() -> None:
    signword = load_signword()
    schedule = decode_schedule()
    print("THM-3394 INERT TWELVE-ORDER HADAMARD BANK")
    print(
        f"signword length={len(signword)} plus={signword.count('+')} minus={signword.count('-')} "
        f"sha256={SIGNWORD_SHA256}"
    )
    print(f"signword_zlib length=3488 sha256={SIGNWORD_ZLIB_SHA256}")
    print(f"schedule bits=192 records={len(schedule)} sha256={SCHEDULE_SHA256}")

    cursor = 0
    first_matrix: list[str] | None = None
    for index, (kind, primary_length) in enumerate(schedule, 1):
        primary_start = cursor
        primary = signword[cursor:cursor + primary_length]
        cursor += primary_length
        require(len(primary) == primary_length, "truncated primary payload")
        side_ranges: list[tuple[int, int]] = []
        if kind in COMPOSITE_SPEC:
            side_total, _, _, _, side_lengths, _, _ = COMPOSITE_SPEC[kind]
            lengths = (side_total,) + side_lengths
            pieces = []
            for length in lengths:
                start = cursor
                pieces.append(signword[cursor:cursor + length])
                cursor += length
                require(len(pieces[-1]) == length, "truncated composite side payload")
                side_ranges.append((start, cursor))
            rows = render_composite(kind, primary, pieces)
            degree_three, degree_four = oa_contract(pieces[1:])
            contract = (
                f"OA={len(pieces[1])} d1=d2=0 "
                f"d3={','.join(map(str, degree_three))} d4={degree_four}"
            )
        else:
            rows = render_h(primary, unbordered=(kind == "g"))
            row_sums, paf = periodic_contract(primary, kind)
            contract = (
                f"period={primary_length // 4} rowsums={','.join(map(str, row_sums))} "
                f"offpeak_PAF_sum={paf}"
            )
        order = len(rows)
        raw_hash, normal_hash, pairs, minimum, maximum = verify_matrix(rows, order)
        require(order in EXPECTED, "unexpected output order")
        if first_matrix is None:
            first_matrix = rows
        print(
            f"{index:02d} kind={kind} primary=[{primary_start},{primary_start + primary_length}) "
            f"sides={side_ranges or '-'} consumed={cursor} shape={order}x{order}"
        )
        print(f"   contract {contract}")
        print(
            f"   raw_sha256={raw_hash} normalized_sha256={normal_hash} "
            f"pairs={pairs} hamming_min={minimum} hamming_max={maximum}"
        )

    require(cursor == len(signword) == 23828, "schedule did not consume exactly the signword")
    require(first_matrix is not None, "empty matrix bank")
    hostile = single_entry_hostile(first_matrix)
    print(
        f"hostile single_entry_flip=REJECTED pair={hostile[0]},{hostile[1]} "
        f"hamming={hostile[2]} expected=334"
    )
    print("consumption exact=23828/23828")
    print("verdict all_twelve_HH_transpose_equals_order_I=YES")
    print("scope finite_orders_through_2000_only; full_Hadamard_conjecture=NOT_CLAIMED")


if __name__ == "__main__":
    main()
