##################################################
FC = ifort #gfortran
FFLAGS=-O3 -assume byterecl -heap-arrays -nogen-interfaces -mkl=parallel -mcmodel=large -fPIC -shared-intel -g -traceback -lpthread -parallel
INC = -I/opt/intel/include -I/opt/local/include
LIB = -L/opt/intel/lib -L/opt/intel/mkl/lib -L/opt/local/lib

#=========
TARGETS=\
intpt_EGM96\
intpt_EGM08\
allocate_VS\
regional_biftag\
txt_vector\
print_rivvec

# Suffix
.SUFFIXES : .f90
.f90:
	$(FC) ${FFLAGS} ${INC} ${LIB} $^ -o $@

# commands
all: $(TARGETS)
# clean
clean:
	rm -rf  $(TARGETS)


# # .SUFFIXES : .o .f90 .F
# # .f90.o:
# # 	$(FC) -c ${FFLAGS} ${INC} $(MODS) $*.f90 -o $*.o
# # src: $(OBJECTS)
# # 	ar -rv srcs.a $(OBJECTS)
# # ##################################################
# # # simulation  
# # ##################################################
# # TARGET1 = txt_vector
# # TARGET2 = print_rivvec
# # TARGET3 = allocate_VS

# # ##################################################
# # $(TARGET1) : $(TARGET1).o $(LIBS)
# # 	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

# # $(TARGET2) : $(TARGET2).o $(LIBS)
# # 	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)

# # $(TARGET3) : $(TARGET3).o $(LIBS)
# # 	$(FC) $(LFLAGS) $@.o $(LIBS) $(MOD) -o $@ $(LIB)
# # ##################################################
# # # Commands
# # ##################################################
# # all : ${TARGET1} ${TARGET2} ${TARGET3} ##${TARGET4}
# # clean: 
# # 	$(RM) -f *.o *.s core *~ *trace *.mod

