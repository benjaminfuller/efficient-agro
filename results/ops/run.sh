for ds in amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k; do
  ./target/release/sparq_bench --dataset datasets/$ds.pts --kinds group,semi,quantile --queries 500 --fanin-cap 30 > results/ops/$ds.cap30.log 2>&1
done
echo done > results/ops/DONE
