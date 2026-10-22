#!/usr/bin/env python3
"""Audit helper (audit_ssj): convert an (untrusted) downloaded HTML file to plain text
and report keyword hits.  The HTML is only parsed as data (html.parser), never executed.
Usage: python3 -I html_to_text.py <in.html> <out.txt> [keyword ...]
"""
import sys
import re
from html.parser import HTMLParser


class Strip(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out = []
        self.skip = 0

    def handle_starttag(self, tag, attrs):
        if tag in ("script", "style"):
            self.skip += 1
        if tag in ("p", "div", "br", "li", "h1", "h2", "h3", "h4", "h5", "h6", "tr", "section", "math"):
            self.out.append("\n")
        # keep the alttext of MathML, which carries the LaTeX source
        if tag == "math":
            for k, v in attrs:
                if k == "alttext" and v:
                    self.out.append(" $" + v + "$ ")
                    self.skip += 1  # skip the presentation MathML inside

    def handle_endtag(self, tag):
        if tag in ("script", "style"):
            self.skip -= 1
        if tag == "math":
            self.skip = max(self.skip - 1, 0)

    def handle_data(self, data):
        if self.skip == 0:
            self.out.append(data)


def main():
    src, dst = sys.argv[1], sys.argv[2]
    kws = sys.argv[3:]
    with open(src, encoding="utf-8", errors="replace") as fh:
        html = fh.read()
    p = Strip()
    p.feed(html)
    text = "".join(p.out)
    text = re.sub(r"[ \t]+", " ", text)
    text = re.sub(r"\n\s*\n+", "\n\n", text)
    with open(dst, "w", encoding="utf-8") as fh:
        fh.write(text)
    print(f"text length {len(text)} chars written to {dst}")
    low = text.lower()
    for kw in kws:
        hits = [m.start() for m in re.finditer(re.escape(kw.lower()), low)]
        print(f"keyword {kw!r}: {len(hits)} hits")
        for h in hits[:6]:
            ctx = text[max(0, h - 90): h + 90].replace("\n", " ")
            print(f"    ...{ctx}...")


if __name__ == "__main__":
    main()
