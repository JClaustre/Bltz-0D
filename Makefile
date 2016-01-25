# **** Executable
exe	= run_BOD
all 	: $(exe)
# ******************************************************************** #

# **** Compilateur
#FC = gfortran
FC = gcc
# ******************************************************************** #


# **** Basic directory
ROOTDIR		?=.
F90_SRCDIR	?=$(ROOTDIR)/SRC_F90
OBJ_DIR		?=$(ROOTDIR)/obj
F90_OBJDIR	?=$(OBJ_DIR)/o.f90
MOD_DIR		?=$(F90_OBJDIR)
# ******************************************************************** #

# **** Options de compilation
USEDEBUGG       := 0
USEOPENMP	:= 1
# ******************************************************************** #

# ******************************************************************** #
ifeq ($(USEOPENMP),0)
	FLAGS += -fopenmp
endif
# ******************************************************************** #
ifeq ($(USEDEBUGG),1)
	FLAGS += -g -O0 -Wall -fbounds-check -fbacktrace -fbackslash
else
	FLAGS += -O2 -fbackslash -march=native -funroll-loops -pipe
endif
# ******************************************************************** #

LIBS	= -lgfortran -lm -fopenmp -fbackslash
# ******************************************************************** #
# ******************************************************************** #

# **** Rules
.SUFFIXES:
.SUFFIXES: .f90 .o 
$(F90_OBJDIR)/%.o	: $(F90_SRCDIR)/%.f90
	@$(FC) $(FLAGS) -o $@ -c $< -J $(MOD_DIR)
	@echo "--> Compiled: " $^

SRCF90	= f90_kind.f90 param.f90 excit.f90 l-xchange.f90 penn_asso.f90\
	recomb.f90 ioniz.f90 radiff.f90 heat.f90 read_input.f90 \
	evolution.f90 main.f90

F90_SRC := $(addprefix $(F90_SRCDIR)/, $(SRCF90))
F90_OBJ := $(F90_SRC:$(F90_SRCDIR)/%.f90=$(F90_OBJDIR)/%.o)

$(exe)	: $(F90_OBJ)
	@$(FC) $(F90_OBJ) $(LIBS) -o $@
	@echo "~~> Make Exec : " [$(exe)]
	@echo "*** Include   : " [$(FLAGS)]
	@echo "*** Library   : " [$(LIBS)]
# ******************************************************************** #
makedirectories:
	$(VERBOSE)mkdir -p $(OBJ_DIR)
	$(VERBOSE)mkdir -p $(F90_OBJDIR)
	$(VERBOSE)mkdir -p $(MOD_DIR)

.PHONY	: clean, superClean, hyperClean

clean	:
	@rm $(F90_OBJDIR)/*.o
	@rm $(MOD_DIR)/*.mod

superClean	: clean
	@rm -rf $(exe)

hyperClean	: superClean
	@rm $(F90_SRCDIR)/*.f90~
	@rm $(ROOTDIR)/datFile/*.dat
env :
	echo $(OBJ)
