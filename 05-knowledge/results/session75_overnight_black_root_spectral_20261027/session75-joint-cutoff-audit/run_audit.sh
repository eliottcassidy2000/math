#!/bin/sh
set -eu
PACKET=${1:?Usage: run_audit.sh PRODUCER_PACKET REPOSITORY_CHECKOUT}
REPO=${2:?Usage: run_audit.sh PRODUCER_PACKET REPOSITORY_CHECKOUT}
HERE=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
export PYTHONDONTWRITEBYTECODE=1
python "$HERE/independent_audit.py" --packet "$PACKET" --repo "$REPO" --out "$HERE/independent_receipt.json"
python -O "$HERE/independent_audit.py" --packet "$PACKET" --repo "$REPO" --out "$HERE/independent_receipt.optimized.json"
cmp "$HERE/independent_receipt.json" "$HERE/independent_receipt.optimized.json"
python "$HERE/supplementary_checks.py" --packet "$PACKET" --out "$HERE/supplementary_receipt.json"
python -O "$HERE/supplementary_checks.py" --packet "$PACKET" --out "$HERE/supplementary_receipt.optimized.json"
cmp "$HERE/supplementary_receipt.json" "$HERE/supplementary_receipt.optimized.json"
