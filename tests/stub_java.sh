#!/usr/bin/env bash
# Fake "java" for the smoke test: stands in for `java -jar juicer_tools.jar pre ...`
# without needing a JRE or the juicer_tools jar. It just finds the .hic output
# argument and writes a non-empty placeholder so 03_Generate_hic_and_TADs.r can
# complete its own PRE-format validation (the part we actually want to test).
set -euo pipefail
for a in "$@"; do
  case "$a" in
    *.hic) printf 'stub-hic\n' > "$a" ;;
  esac
done
exit 0
