#!/usr/bin/env bash
# Full efficient-agro sweep. Build happens in a dry-run SAM (no crypto); with
# --crypto the snapshot is installed into an encrypted Path OSAM+ tree and the
# queries run there. Results are appended to results/<mode>.csv and one log
# per run is kept in results/.
#
#   ./run_all.sh            # dry-run counts for every dataset
#   ./run_all.sh --crypto   # encrypted Path OSAM+ (needs a large-memory server)
#
# Environment overrides: QUERIES (default 1000), CAPS ("raw"),
# KINDS_SMALL / KINDS_LARGE (the median tail is very large on 1M-point sets).
set -euo pipefail
cd "$(dirname "$0")"
MODE_FLAG="${1:-}"
MODE=$([ "$MODE_FLAG" = "--crypto" ] && echo crypto || echo dryrun)
QUERIES="${QUERIES:-1000}"
CAPS="${CAPS:-raw}"
KINDS_SMALL="${KINDS_SMALL:-group,semi,quantile}"
KINDS_LARGE="${KINDS_LARGE:-group,semi}"
BIN=./target/release/sparq_bench
cargo build --release
mkdir -p results

SMALL="amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k"
LARGE="gowalla-1d synthetic_2d_1m synthetic_2d_1m_sparse synthetic_2d_1m-1024x1024 synthetic_3d_1m_128 synthetic_3d_1m_256"

run() {
  local ds=$1 kinds=$2 cap=$3
  local extra=()
  [ "$cap" != raw ] && extra=(--fanin-cap "$cap")
  echo "== $ds cap=$cap kinds=$kinds mode=$MODE"
  $BIN --dataset "datasets/$ds.pts" --kinds "$kinds" --queries "$QUERIES" --check \
       ${MODE_FLAG:+$MODE_FLAG} "${extra[@]}" --csv "results/$MODE.csv" \
       > "results/${ds}.${MODE}.cap${cap}.log" 2>&1 || echo "   FAILED (see log)"
  grep '^result' "results/${ds}.${MODE}.cap${cap}.log" | cut -c1-160 || true
}

for ds in $SMALL; do for cap in $CAPS; do run "$ds" "$KINDS_SMALL" "$cap"; done; done
for ds in $LARGE; do for cap in $CAPS; do run "$ds" "$KINDS_LARGE" "$cap"; done; done
