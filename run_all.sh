#!/usr/bin/env bash
# efficient-agro sweep. The structure is always built in a dry-run SAM. With
# --crypto, the snapshot is bulk-loaded into an encrypted Path OSAM+ tree
# and the queries run there.
#
#   ./run_all.sh            # dry-run read counts: every dataset, every tail,
#                           # 1000 queries, one run at a time
#   ./run_all.sh --crypto   # encrypted timing: every dataset, every tail and
#                           # SPARQ, 100 queries
#
# In crypto mode every (dataset, scheme, tail) is its own run with its own
# structure, so each encrypted tree holds only what its queries need (a group
# run carries no segment or wavelet trees). The encrypted tree is held in
# memory: 4 * 256 B per leaf, i.e. 1.07 GB at 2^20 leaves and 68.7 GB at 2^26.
# Small datasets run JOBS at a time, then the 1M-point datasets JOBS_LARGE at
# a time. MAX_SERVER_GB makes a run skip (with a "skip" line in its log)
# instead of allocating a larger tree.
#
# Environment overrides (defaults: dry-run / crypto):
#   QUERIES        1000 / 100
#   KINDS          layered tails for the small datasets: group,semi,quantile
#   KINDS_LARGE    tails for the 1M-point datasets: group,semi
#   LARGE          include the 1M-point datasets: 1 / 1
#   DATASETS       explicit dataset list (overrides the small/large lists;
#                  runs with the small-dataset settings)
#   JOBS           small-dataset runs in parallel: 1 / 192
#   JOBS_LARGE     1M-point runs in parallel: 1 / 1
#   MAX_SERVER_GB  crypto: skip runs whose encrypted tree exceeds this (unset = no limit)
#   SKIP_DONE      crypto: skip runs whose log already has a result: 1
#   CAPS           fan-in caps to sweep: "raw"
#
# Logs: dry-run  results/dryrun/<dataset>.cap<cap>.log (both schemes, all tails)
#       crypto   results/crypto/<dataset>.<scheme>-<tail>.cap<cap>.log
#                (<scheme>-<tail> is layered-group, layered-semi,
#                layered-quantile or sparq)
# Per-run CSVs are merged into results/<mode>.csv at the end.
set -euo pipefail
cd "$(dirname "$0")"
MODE_FLAG="${1:-}"
KINDS="${KINDS:-group,semi,quantile}"
KINDS_LARGE="${KINDS_LARGE:-group,semi}"
LARGE="${LARGE:-1}"
JOBS_LARGE="${JOBS_LARGE:-1}"
if [ "$MODE_FLAG" = "--crypto" ]; then
  MODE=crypto
  QUERIES="${QUERIES:-100}"
  JOBS="${JOBS:-192}"
  SKIP_DONE="${SKIP_DONE:-1}"
else
  MODE=dryrun
  MODE_FLAG=""
  QUERIES="${QUERIES:-1000}"
  JOBS="${JOBS:-1}"
  SKIP_DONE=0
fi
MAX_SERVER_GB="${MAX_SERVER_GB:-}"
CAPS="${CAPS:-raw}"
BIN=./target/release/sparq_bench
cargo build --release
OUT="results/$MODE"
mkdir -p "$OUT"

SMALL_LIST="amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k"
LARGE_LIST="gowalla-1d synthetic_2d_1m synthetic_2d_1m_sparse synthetic_2d_1m-1024x1024 synthetic_3d_1m_128 synthetic_3d_1m_256"

# One line per run: dataset, scheme, tails, cap, label.
small_jobs=$(mktemp)
large_jobs=$(mktemp)
trap 'rm -f "$small_jobs" "$large_jobs"' EXIT
add() { # file dataset kinds
  local file=$1 ds=$2 kinds=$3 cap k
  for cap in $CAPS; do
    if [ "$MODE" = crypto ]; then
      for k in ${kinds//,/ }; do echo "$ds layered $k $cap layered-$k" >> "$file"; done
      echo "$ds sparq group $cap sparq" >> "$file"
    else
      echo "$ds both $kinds $cap all" >> "$file"
    fi
  done
}
if [ -n "${DATASETS:-}" ]; then
  for ds in $DATASETS; do add "$small_jobs" "$ds" "$KINDS"; done
else
  for ds in $SMALL_LIST; do add "$small_jobs" "$ds" "$KINDS"; done
  if [ "$LARGE" = 1 ]; then
    for ds in $LARGE_LIST; do add "$large_jobs" "$ds" "$KINDS_LARGE"; done
  fi
fi

run() {
  local ds=$1 scheme=$2 kinds=$3 cap=$4 label=$5
  local extra=()
  [ "$cap" != raw ] && extra=(--fanin-cap "$cap")
  [ -n "$MAX_SERVER_GB" ] && [ -n "$MODE_FLAG" ] && extra+=(--max-server-gb "$MAX_SERVER_GB")
  local log csv
  if [ "$MODE" = crypto ]; then
    log="$OUT/${ds}.${label}.cap${cap}.log"
    csv="$OUT/${ds}.${label}.cap${cap}.csv"
  else
    log="$OUT/${ds}.cap${cap}.log"
    csv="$OUT/${ds}.cap${cap}.csv"
  fi
  if [ "$SKIP_DONE" = 1 ] && [ -f "$log" ] && grep -q '^result' "$log"; then
    echo "== skip  $ds $label (already done: $log)"
    return 0
  fi
  rm -f "$csv"
  local start=$(date +%s)
  echo "== start $ds $label cap=$cap mode=$MODE queries=$QUERIES"
  if $BIN --dataset "datasets/$ds.pts" --scheme "$scheme" --kinds "$kinds" --queries "$QUERIES" --check \
       ${MODE_FLAG:+$MODE_FLAG} ${extra[@]+"${extra[@]}"} --csv "$csv" \
       > "$log" 2>&1; then
    echo "== done  $ds $label ($(( $(date +%s) - start )) s)"
  else
    echo "== FAILED $ds $label after $(( $(date +%s) - start )) s (see $log)"
  fi
  grep -E '^(result|skip)' "$log" | cut -c1-160 || true
}
export -f run
export BIN OUT MODE MODE_FLAG QUERIES MAX_SERVER_GB SKIP_DONE

echo "$(wc -l < "$small_jobs") small runs, $JOBS at a time"
xargs -P "$JOBS" -L 1 bash -c 'run "$0" "$1" "$2" "$3" "$4"' < "$small_jobs"
if [ -s "$large_jobs" ]; then
  echo "$(wc -l < "$large_jobs") large runs, $JOBS_LARGE at a time"
  xargs -P "$JOBS_LARGE" -L 1 bash -c 'run "$0" "$1" "$2" "$3" "$4"' < "$large_jobs"
fi

# Merge the per-run CSVs (one header).
merged="results/$MODE.csv"
first=1
: > "$merged"
if [ "$MODE" = crypto ]; then
  csvs=("$OUT"/*.layered-*.cap*.csv "$OUT"/*.sparq.cap*.csv)   # this script's runs only
else
  csvs=("$OUT"/*.csv)
fi
for f in "${csvs[@]}"; do
  [ -e "$f" ] || continue
  if [ $first = 1 ]; then cat "$f" >> "$merged"; first=0; else tail -n +2 "$f" >> "$merged"; fi
done
echo "merged CSV: $merged"
