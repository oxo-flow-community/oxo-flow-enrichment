#!/usr/bin/env bash
# Acceptance test for the oxo-flow-enrichment port.
# Usage: ./test/run.sh            (uses ./main.oxoflow)
set -euo pipefail
cd "$(dirname "$0")/.."
OXO=${OXO:-oxo-flow}

echo "==> validate"
"$OXO" validate main.oxoflow

echo "==> lint (warnings are acceptable, errors are not)"
"$OXO" lint main.oxoflow

echo "==> dry-run with default config"
"$OXO" dry-run main.oxoflow > /tmp/oxo-dryrun-$$.txt 2>&1
grep -q "would execute" /tmp/oxo-dryrun-$$.txt

echo "==> debug: expanded commands contain no literal placeholders"
# WARN log lines are excluded: the ORA_GSEApy/RcisTarget aggregates carry an
# intentional txt_gene_sets fan that contributes zero inputs on the default
# config, and the engine's plan-time note for that quotes the raw pattern.
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1
grep -v WARN /tmp/oxo-debug-$$.txt | grep -Eq '\{config\.|\{region_set\}|\{gene_set\}' && { echo "unexpanded placeholders in debug output"; exit 1; } || true

echo "PASS"
