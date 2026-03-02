# Dont forget to change the gurobi lib paths
clear
mkdir -p log
mkdir -p Results
mkdir -p code_formal

BLUE='\033[0;34m'
NC='\033[0m'
echo -e "${BLUE}Begin task: - $(date '+%Y-%m-%d %H:%M:%S')${NC}"
INDIM=32

FILENAME=$1
HH=$2

echo ${FILENAME}
echo ${HH}


runname=main_${FILENAME}_${HH}

LOG_FILE=log/log_${runname}_$(date '+%Y-%m-%d-%H:%M:%S').log

rm ${runname}
g++ -march=native -O3 -fopenmp -std=c++17 main.cpp -DHH=${HH} -DINDIM=${INDIM} -o ${runname} -I ${GUROBI_HOME}/include -L ${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi120 -lm
./${runname} ${FILENAME} |tee ${LOG_FILE} 



echo -e "${BLUE}Finish task: - $(date '+%Y-%m-%d %H:%M:%S')${NC}"