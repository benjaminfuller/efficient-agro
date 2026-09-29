for ds in cali-1024x1024 spitz-1024x1024 nh_64 gowalla_2d_100k; do
 for cfg in "2 raw" "2 2" "4 raw" "4 4" "8 8" "30 8" "30 16" "30 raw"; do
  set -- $cfg; b=$1; c=$2; extra=""; [ "$c" != raw ] && extra="--fanin-cap $c"
  ./target/release/sparq_bench --dataset datasets/$ds.pts --kinds group --queries 500 --branching $b $extra > results/ops/$ds.b$b.cap$c.log 2>&1
 done
done
echo done > results/ops/SWEEPDONE
