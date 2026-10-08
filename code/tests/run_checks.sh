#!/bin/sh
# Run from code/: sh tests/run_checks.sh (on macOS or Linux).
set -eu
CC=${CC:-gcc}
check_binary=$(mktemp "${TMPDIR:-/tmp}/dpd-check.XXXXXX")
trap 'rm -f "$check_binary"' EXIT HUP INT TERM
sources=""
for source in ./*.c; do
    [ "$source" = "./main.c" ] || sources="$sources $source"
done
# Source filenames are fixed by this project and contain no spaces.
$CC ${CFLAGS:--O2 -Wall -Wextra -Wpedantic} -I. $sources tests/test_dpd.c -o "$check_binary" -lm
"$check_binary"
