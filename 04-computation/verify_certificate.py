#!/usr/bin/env python3
"""Compatibility wrapper for exact LRC certificate verification."""

from __future__ import annotations

import json
import sys
from pathlib import Path

from lrc_exact_tools import verify_certificate


def main() -> int:
    if len(sys.argv) != 2:
        print(json.dumps({"ok": False, "reason": "usage", "usage": "python3 verify_certificate.py certificate.json"}))
        return 2
    print(json.dumps(verify_certificate(Path(sys.argv[1])), sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
