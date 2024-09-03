ls -l ./ | awk '{print $NF}' | while read frqcase
do
  if [ -d ${frqcase} ]; then
    echo -n "make ${frqcase} ... "
    start=$(date +%s)
    cd ${frqcase}
    wmake -j4
    #mkdir test
    cd ..
  fi
done
