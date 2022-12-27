  
home=/home/deepl/sooho/test/K-com/test
 echo -ne "\033[1;36m(pyver execute)\033[0m                                    \r"
 start1=`date +%s.%N`
 python ${home}/pyver.py
 end1=`date +%s.%N`
 runtime1=$( echo "$end1 - $start1" | bc -l )
 echo -ne "\033[1;36m(pyver execute)\033[0m    done!($(printf %.4f $runtime1))   \n"
