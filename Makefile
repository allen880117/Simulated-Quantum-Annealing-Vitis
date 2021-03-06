HLS_INCLUDE=/home/edci/Xilinx/Vitis_HLS/2021.1/include/
PYTHON_PATH=/usr/include/python3.8

CFLAG=-g -std=c++14 -I$(HLS_INCLUDE) -I$(XILINX_XRT)/include -I$(PYTHON_PATH) -L$(XILINX_XRT)/lib
CFLAG_POST=-lxrt_coreutil -pthread -lpython3.8 -DWITHOUT_NUMPY 

SQA_FLAG=-DAM=0 -DREPLAY=0 -DU50=0 -DCUR_SIM=2

SRC=src/csim_host/main.cpp src/csim_helper/helper.cpp src/kernel/qmc_basic.cpp 

all: naive.exe

naive.exe: $(SRC)
	g++ $(CFLAG) -o $@ $(SRC) $(CFLAG_POST) $(SQA_FLAG)

clean:
	rm -f naive.exe
	rm -f time_log.basic.txt out_log.basic.txt