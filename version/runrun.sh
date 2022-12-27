cd cudagraphtest/test
make
cd ..
cd ..
start1=`date +%s.%N`
exec/run.x
end1=`date +%s.%N`
runtime1=$( echo "$end1 - $start1" | bc -l )
echo  "\033[1;36m(cugraphver execute)\033[0m    done!($(printf %.4f $runtime1))   \n"

DD1=$(diff -y --suppress-common-lines result/test2.txt result/gtest.txt | wc -l)
echo "diff: ${DD1}"