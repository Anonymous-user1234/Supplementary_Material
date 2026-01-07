#!/bin/bash

clear
mkdir -p log


OPS_FILE=$1

DEPTH_LIMIT=$2


TIE_BREAKING=$3 #1:MN 2:ASM 3:SOM

TAR_STRATEGY=$4 #1:AT 2:ST


RUN=bp66_depth${DEPTH_LIMIT}_${OPS_FILE}_tie${TIE_BREAKING}_tar${TAR_STRATEGY}



if [ $OPS_FILE = "aes10_ops" ]
then
	FILENAME=\"AES10\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=32
fi

if [ $OPS_FILE = "aes10v2_ops" ]
then
	FILENAME=\"AES10v2\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=32
fi
if [ $OPS_FILE = "aes10v3_ops" ]
then
	FILENAME=\"AES10v3\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=33
fi
if [ $OPS_FILE = "aes12_ops" ]
then
	FILENAME=\"AES12\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=34
fi
if [ $OPS_FILE = "aes12v2_ops" ]
then
	FILENAME=\"AES12v2\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=32
fi
if [ $OPS_FILE = "aes12v3_ops" ]
then
	FILENAME=\"AES12v3\"
	SBOX_FILE=aes
	DIM=8
	NL_NUM=33
fi
if [ $OPS_FILE = "saturnin_ops" ]
then
	FILENAME=\"Saturnin\"
	SBOX_FILE=saturn
	DIM=16
	NL_NUM=48
fi
if [ $OPS_FILE = "snow_ops" ]
then
	FILENAME=\"SNOW\"
	SBOX_FILE=snow3g
	DIM=8
	NL_NUM=90
fi


THREADS=32

echo "OPS_FILE: ${OPS_FILE}, DIM: ${DIM}, SBOX_FILE: ${SBOX_FILE}, TIE_BREAKING: ${TIE_BREAKING}, Target Condition: ${TAR_STRATEGY}"



LOG_FILE=log/${RUN}_${FILENAME}_dim${DIM}_tie_${TIE_BREAKING}_tc${TAR_STRATEGY}.log

rm -f ${RUN}.run

g++ -O3 -std=c++17 -m64 -g -o ${RUN}.run src/main.cpp  -DFILENAME=${OPS_FILE} -DWORKERS=${THREADS} -DDIMENSION=${DIM} -DNL_NUM=${NL_NUM} -DDEPTH_LIMIT=${DEPTH_LIMIT} -DSTRATEGY=${TIE_BREAKING} -DTARGET_COND=${TAR_STRATEGY} -I/${GUROBI_HOME}/include/ -L/${GUROBI_HOME}/lib -lgurobi_c++ -lgurobi120 -lm -lpthread


./${RUN}.run s-box/${SBOX_FILE}.txt data/${OPS_FILE}.txt | tee  ${LOG_FILE}


