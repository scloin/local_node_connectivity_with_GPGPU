#!/bin/bash
N=1500
for i in {1..50}
do
#pp=$((i * 0.0002))
pp=$(bc <<< "scale=10; $i/500")
./timecheck.sh -n $N -p $pp >> probcheck3.txt

done