for ds in amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k; do
  kinds=group,semi,quantile; [ $ds = gowalla_2d_100k ] && kinds=group,semi
  ./target/release/sparq_bench --dataset datasets/$ds.pts --kinds $kinds --queries 500 --check > results/balanced/$ds.raw.log 2>&1
done
echo done > results/balanced/DONE
