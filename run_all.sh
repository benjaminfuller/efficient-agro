#!/usr/bin/env bash
# efficient-agro sweep. The structure is always built in a dry-run SAM. With
# --crypto, the snapshot is installed into an encrypted Path OSAM+ tree and
# the queries run there.
#
#   ./run_all.sh            # dry-run read counts: every dataset, every tail,
#                           # 1000 queries, one run at a time
#   ./run_all.sh --crypto   # encrypted timing and stash: small datasets,
#                           # group tail + SPARQ, 100 queries, 4 runs in parallel
#
# Read counts are identical in the two modes, so the encrypted sweep only
# measures wall-clock time and stash size, and is kept small. Every query is
# padded to its worst-case budget, so the semi and median tails, and the 3D
# datasets, dominate the encrypted cost.
#
# Environment overrides (defaults: dry-run / crypto):
#   QUERIES        1000 / 100
#   KINDS          layered tails for the small datasets: group,semi,quantile / group
#   KINDS_LARGE    tails for the 1M-point datasets: group,semi / group
#   EXTRA_SEMI     crypto only: datasets that also run the semi tail
#                  (default "spitz-1024x1024 nh_64"; "" for none)
#   LARGE          include the 1M-point datasets: 1 / 0
#   DATASETS       explicit dataset list (overrides the small/large lists)
#   JOBS           runs in parallel: 1 / 4 (each run holds a whole tree in memory)
#   CAPS           fan-in caps to sweep: "raw"
#
# Logs go to results/<mode>/<dataset>.<cap>.log. Per-run CSVs are merged
# into results/<mode>.csv at the end.
set -euo pipefail
cd "$(dirname "$0")"
MODE_FLAG="${1:-}"
if [ "$MODE_FLAG" = "--crypto" ]; then
  MODE=crypto
  QUERIES="${QUERIES:-100}"
  KINDS="${KINDS:-group}"
  KINDS_LARGE="${KINDS_LARGE:-group}"
  EXTRA_SEMI="${EXTRA_SEMI-spitz-1024x1024 nh_64}"
  LARGE="${LARGE:-0}"
  JOBS="${JOBS:-192}"
else
  MODE=dryrun
  MODE_FLAG=""
  QUERIES="${QUERIES:-1000}"
  KINDS="${KINDS:-group,semi,quantile}"
  KINDS_LARGE="${KINDS_LARGE:-group,semi}"
  EXTRA_SEMI=""
  LARGE="${LARGE:-1}"
  JOBS="${JOBS:-1}"
fi
CAPS="${CAPS:-raw}"
BIN=./target/release/sparq_bench
cargo build --release
OUT="results/$MODE"
mkdir -p "$OUT"

SMALL_LIST="amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k"
LARGE_LIST="gowalla-1d synthetic_2d_1m synthetic_2d_1m_sparse synthetic_2d_1m-1024x1024 synthetic_3d_1m_128 synthetic_3d_1m_256"

# One line per run: dataset, tails, cap.
jobs_file=$(mktemp)
trap 'rm -f "$jobs_file"' EXIT
add() { for cap in $CAPS; do echo "$1 $2 $cap" >> "$jobs_file"; done; }
kinds_for() {
  local ds=$1 kinds=$2
  for s in $EXTRA_SEMI; do
    if [ "$s" = "$ds" ] && [[ ",$kinds," != *",semi,"* ]]; then kinds="$kinds,semi"; fi
  done
  echo "$kinds"
}
if [ -n "${DATASETS:-}" ]; then
  for ds in $DATASETS; do add "$ds" "$(kinds_for "$ds" "$KINDS")"; done
else
  for ds in $SMALL_LIST; do add "$ds" "$(kinds_for "$ds" "$KINDS")"; done
  if [ "$LARGE" = 1 ]; then
    for ds in $LARGE_LIST; do add "$ds" "$(kinds_for "$ds" "$KINDS_LARGE")"; done
  fi
fi

run() {
  local ds=$1 kinds=$2 cap=$3
  local extra=()
  [ "$cap" != raw ] && extra=(--fanin-cap "$cap")
  local log="$OUT/${ds}.cap${cap}.log"
  local start=$(date +%s)
  echo "== start $ds kinds=$kinds cap=$cap mode=$MODE queries=$QUERIES"
  if $BIN --dataset "datasets/$ds.pts" --kinds "$kinds" --queries "$QUERIES" --check \
       ${MODE_FLAG:+$MODE_FLAG} ${extra[@]+"${extra[@]}"} --csv "$OUT/${ds}.cap${cap}.csv" \
       > "$log" 2>&1; then
    echo "== done  $ds ($(( $(date +%s) - start )) s)"
  else
    echo "== FAILED $ds after $(( $(date +%s) - start )) s (see $log)"
  fi
  grep '^result' "$log" | cut -c1-160 || true
}
export -f run
export BIN OUT MODE MODE_FLAG QUERIES

echo "$(wc -l < "$jobs_file") runs, $JOBS at a time"
xargs -P "$JOBS" -L 1 bash -c 'run "$0" "$1" "$2"' < "$jobs_file"

# Merge the per-run CSVs (one header).
merged="results/$MODE.csv"
first=1
: > "$merged"
for f in "$OUT"/*.csv; do
  [ -e "$f" ] || continue
  if [ $first = 1 ]; then cat "$f" >> "$merged"; first=0; else tail -n +2 "$f" >> "$merged"; fi
done
echo "merged CSV: $merged"
