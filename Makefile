.PHONY : document build cleantmp cleanexe cleanobs cleanres cleanall

document:
	@echo
	@echo "#################### Package Document ####################"
	@echo
	@echo "This is the file to complete the structure of the project,"
	@echo "and complier the code into executive file."
	@echo
	@echo "First of all, use 'make build' to build the storage structure;"
	@echo
	@echo "Type 'make cleanall' to avoid impact of previous work;"
	@echo
	@echo "TYPE 'make' + purpose to generate executive file."
	@echo "For example, 'make forward' to complier code for forward modeling"
	@echo "Here we offer three options"
	@echo "forward -- forward modeling"
	@echo "migrate -- reverse time migration"
	@echo "inverse -- full waveform inversion"
	@echo
	@echo "This package is developed by Shuo Zhang in UTDallas"
	@echo "Contact email: shuo.zhang@utdallas.edu"
	@echo "Copyright belongs to Shuo Zhang and Dr. Hejun Zhu's group"
	@echo
	@echo "##########################################################"
	@echo

# subroutine
DIR_CUR = $(shell pwd)
DIR_SRC = $(DIR_CUR)/src
DIR_BIN = $(DIR_CUR)/bin
DIR_INC = $(DIR_CUR)/include
DIR_LIB = $(DIR_CUR)/library
DIR_TMP = $(DIR_CUR)/temp
DIR_OBS = $(DIR_CUR)/model
DIR_RES = $(DIR_CUR)/result

OBS_SUB = obs sig
RES_SUB = ite ker sig wvf syn
OBS_DET = $(addprefix $(DIR_OBS)/,$(OBS_SUB))
RES_DET = $(addprefix $(DIR_RES)/,$(RES_SUB))

CODE_DIR = $(DIR_INC) $(DIR_LIB) $(DIR_SRC) $(DIR_TMP) $(DIR_BIN)
DATA_DIR = $(RES_DET) $(OBS_DET)
ALL_DIR = $(CODE_DIR) $(DATA_DIR)

# executive file
mod_aco = mod.exe
mig_aco = mig.exe
cfwi_aco = cfwi.exe
rfwi_aco = rfwi.exe
cfwi_aco_amf = cfwi_amf.exe

# complier parameter
CC = icc
FC = gfortran
MPICC = mpiicc
MPIFC = mpif90
CFLAG = -g -I $(DIR_INC) -lfftw3 -lm -O3

# object depandency
object_mod = gene.o tool.o model.o mod_aco.o
object_mig = gene.o tool.o fft.o mig.o mig_aco.o
object_cfwi = gene.o tool.o fft.o amf.o cfwi.o cfwi_aco.o
object_rfwi = gene.o tool.o fft.o rfwi.o rfwi_aco.o
object_cfwi_amf = gene.o tool.o fft.o cfwi.o amf.o cfwi_aco_amf.o

forward : $(mod_aco)

migrate : $(mig_aco)

inverse : $(cfwi_aco)

#---------complier acoustic modeling operater---------#

$(DIR_TMP)/gene.o : $(DIR_LIB)/gene.c $(DIR_INC)/gene.h
	$(MPICC) $(CFLAG) -c $< -o $@
	
$(DIR_TMP)/tool.o : $(DIR_LIB)/tool.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h
	$(MPICC) $(CFLAG) -c $< -o $@
	
$(DIR_TMP)/fft.o : $(DIR_LIB)/fft.c $(DIR_INC)/gene.h $(DIR_INC)/fft.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/amf.o : $(DIR_LIB)/amf.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/amf.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/xcorr.o : $(DIR_LIB)/xcorr.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/xcorr.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier forward modeling  operater----------#

$(DIR_TMP)/model.o : $(DIR_LIB)/model.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h  $(DIR_INC)/model.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/mod_aco.o : $(DIR_SRC)/mod_aco.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/model.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier reverse time migration operator----------#

$(DIR_TMP)/mig.o : $(DIR_LIB)/mig.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/mig.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/mig_aco.o : $(DIR_SRC)/mig_aco.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/mig.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier conventional full waveform inversion operater----------#

$(DIR_TMP)/cfwi.o : $(DIR_LIB)/cfwi.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/cfwi.h $(DIR_INC)/amf.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/cfwi_aco.o : $(DIR_SRC)/cfwi_aco.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/cfwi.h $(DIR_INC)/amf.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier reflection full waveform inversion operater----------#
$(DIR_TMP)/rfwi.o : $(DIR_LIB)/rfwi.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/rfwi.h
	$(MPICC) $(CFLAG) -c $< -o $@

$(DIR_TMP)/rfwi_aco.o : $(DIR_SRC)/rfwi_aco.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/rfwi.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier conventional full waveform inversion operater with amf----------#

$(DIR_TMP)/cfwi_aco_amf.o : $(DIR_SRC)/cfwi_aco_amf.c $(DIR_INC)/gene.h $(DIR_INC)/tool.h $(DIR_INC)/fft.h $(DIR_INC)/cfwi.h $(DIR_INC)/amf.h
	$(MPICC) $(CFLAG) -c $< -o $@

#----------complier executive file----------#

$(mod_aco) : $(patsubst %, $(DIR_TMP)/%, $(object_mod))
	$(MPICC) $^ $(CFLAG) -o $@

$(mig_aco) : $(patsubst %, $(DIR_TMP)/%, $(object_mig))
	$(MPICC) $^ $(CFLAG) -o $@

$(cfwi_aco) : $(patsubst %, $(DIR_TMP)/%, $(object_cfwi))
	$(MPICC) $^ $(CFLAG) -o $@

$(rfwi_aco) : $(patsubst %, $(DIR_TMP)/%, $(object_rfwi))
	$(MPICC) $^ $(CFLAG) -o $@

$(cfwi_aco_amf) : $(patsubst %, $(DIR_TMP)/%, $(object_cfwi_amf))
	$(MPICC) $^ $(CFLAG) -o $@

# complete document structure
build :
	mkdir -p $(ALL_DIR)

# clean files
cleantmp :
	-rm $(DIR_TMP)/*.o -rf
cleanexe :
	-rm $(DIR_CUR)/*.exe -rf
cleanobs :
	-rm $(patsubst %, %/*, $(OBS_DET)) -rf
cleanres :
	-rm $(patsubst %, %/*, $(RES_DET)) -rf
cleanall :
	-rm $(DIR_RES)/*.bin
	-rm $(DIR_CUR)/*.exe $(DIR_TMP)/*.o
	-rm $(patsubst %, %/*, $(OBS_DET)) -rf
	-rm $(patsubst %, %/*, $(RES_DET)) -rf
