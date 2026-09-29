for ds in amazon-books spitz-1024x1024 cali-1024x1024 gowalla_2d_50k gowalla_2d_100k nh_64 gowalla_3d_23k; do
  kinds=group,semi,quantile; case $ds in gowalla_2d_100k|gowalla_3d_23k) kinds=group,semi;; esac
  ./target/release/sparq_bench --dataset datasets/$ds.pts --kinds $kinds --queries 500 --check > results/b64/$ds.log 2>&1
done
./target/release/sparq_bench --dataset datasets/amazon-books.pts --kinds group,semi --queries 20 --check --crypto > results/b64/amazon-books.crypto.log 2>&1
echo done > results/b64/DONE
