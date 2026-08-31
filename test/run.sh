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
# Lines carrying the intentional {gene_set} fan are excluded: the txt_gene_sets
# list is default-empty, so that expand_inputs entry contributes zero inputs.
# oxo-flow <= 0.16.0 injects the raw pattern (with literal {gene_set}) into the
# aggregates' {input}; >= 0.16.1 (#254) returns an empty product and logs a WARN
# quoting the raw pattern instead. Both flavors are filtered here; genuinely
# unexpanded {config.*}/{region_set} placeholders still fail the check.
"$OXO" debug main.oxoflow > /tmp/oxo-debug-$$.txt 2>&1
grep -v WARN /tmp/oxo-debug-$$.txt | grep -v '{gene_set}' | grep -Eq '\{config\.|\{region_set\}' && { echo "unexpanded placeholders in debug output"; exit 1; } || true

echo "PASS"
