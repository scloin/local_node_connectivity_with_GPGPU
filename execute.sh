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

echo -ne "\033[1;36m(compile..)\033[0m             \r"

  echo -ne "\033[1;36m(gengraph)\033[0m         N:"$N" S:"$S" P:"$P"\033[0m\n"
  python ${home}/gengraph.py ${N} ${P} ${S}


for entry in `ls $home | grep .cu$`;
do  
  fn="${entry:0:-3}"
  nvcc -o ${home}/exec/${fn} ${home}/${entry} -I ${home}/.. -rdc=true -g
  echo -ne "\033[1;36m(compile..)\033[0m        done!($fn)              \r"
done
echo -ne "\033[1;36m(compile..)\033[0m        done!                   \n" 

if [ ${C} == 0 ]
then

  start0=`date +%s.%N`
  ${home}/exec/revise4
  end0=`date +%s.%N`
  runtime0=$( echo "($end0 - $start0)" | bc -l )
  echo -ne "\033[1;36m(revise4  execute)\033[0m  done!($(printf %.4f $runtime0))   \n"

  sh runrun.sh


fi