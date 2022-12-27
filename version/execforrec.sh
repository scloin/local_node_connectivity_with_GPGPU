#init
home=/home/deepl/sooho/test/K-com/test
#default value
N=20
P=0.4
S=20
C=0

target=revise5



while getopts n:s:p:c: opt # if first argument str length is not 0
do
  #echo $1
  case ${opt} in
    n) N="$OPTARG";;
    s) S="$OPTARG";;
    p) P="$OPTARG";;
    c) C="$OPTARG";;
    esac
done
  # echo -ne "\033[1;36m(gengraph)\033[0m         N:"$N" S:"$S" P:"$P"\033[0m\n"
  # python ${home}/gengraph.py ${N} ${P} ${S}

#echo -ne "\033[1;36m(compile..)\033[0m             \r"

  #echo -ne "\033[1;36m(gengraph)\033[0m         N:"$N" S:"$S" P:"$P"\033[0m\n"
  python ${home}/gengraph.py ${N} ${P} ${S}


#echo -ne "\033[1;36m(compile..)\033[0m        done!                   \n" 

nvcc -o ${home}/exec/revise4 ${home}/revise4.cu -I ${home}/.. -rdc=true -g
  #echo -ne "\033[1;36m(cuver execute) \033[0m\r"
  # start0=`date +%s.%N`
  # ${home}/exec/${target}
  # end0=`date +%s.%N`
  # runtime0=$( echo "($end0 - $start0)" | bc -l )
  # echo -ne "\033[1;36m(${target} execute)\033[0m  done!($(printf %.4f $runtime0))   \n"
  start0=`date +%s.%N`
  ${home}/exec/revise4
  end0=`date +%s.%N`
  runtime0=$( echo "($end0 - $start0)" | bc -l )
  echo "(revise4 execute)done!($(printf %.4f $runtime0))"
  # start2=`date +%s.%N`
  # ${home}/exec/test1
  # end2=`date +%s.%N`
  # runtime2=$( echo "($end2 - $start2)" | bc -l )
  # echo -ne "\033[1;36m(test1 execute)\033[0m    done!($(printf %.4f $runtime2))   \n"
  # start2=`date +%s.%N`
  # ${home}/exec/test1
  # end2=`date +%s.%N`
  # runtime2=$( echo "($end2 - $start2)" | bc -l )
  # echo -ne "\033[1;36m(test1 execute)\033[0m     done!($(printf %.4f $runtime2))   \n"
  
  #echo -ne "\033[1;36m(pyver execute)\033[0m                                    \r"
  start1=`date +%s.%N`
  python ${home}/pyver.py
  end1=`date +%s.%N`
  runtime1=$( echo "$end1 - $start1" | bc -l )
  echo "(pyver execute)  done!($(printf %.4f $runtime1))"
  
  #echo -ne "\033[1;36m(diff)                              \033[0m\n"
  DD=$(diff -y --suppress-common-lines result/test1.txt result/test2.txt | wc -l)
  DD1=$(diff -y --suppress-common-lines result/test2.txt result/test3.txt | wc -l)
  #echo -ne "\033[1;36m(done)                              \033[0m\n"
  echo "diff: ${DD}"
  #echo "diff: ${DD1}"

