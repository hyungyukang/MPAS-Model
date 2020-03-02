.SUFFIXES: .F .o

OBJS = mpas_ocn_forward_mode.o \
       mpas_ocn_fortrilinos_imp_mod.o \
       mpas_ocn_fortrilinos_model_evaluator_mod.o \
       mpas_ocn_time_integration.o \
       mpas_ocn_time_integration_rk4.o \
       mpas_ocn_time_integration_si.o \
       mpas_ocn_time_integration_split.o

all: forward_mode

forward_mode: $(OBJS) 

mpas_ocn_time_integration.o: mpas_ocn_time_integration_rk4.o mpas_ocn_time_integration_si.o  mpas_ocn_time_integration_split.o

mpas_ocn_time_integration_rk4.o: 

mpas_ocn_time_integration_si.o: mpas_ocn_fortrilinos_imp_mod.o

mpas_ocn_time_integration_split.o: 

mpas_ocn_fortrilinos_imp_mod.o: mpas_ocn_fortrilinos_model_evaluator_mod.o

mpas_ocn_fortrilinos_model_evaluator_mod.o:

mpas_ocn_forward_mode.o: mpas_ocn_time_integration.o \
                         mpas_ocn_time_integration_rk4.o \
                         mpas_ocn_time_integration_si.o \
                         mpas_ocn_time_integration_split.o

clean:
	$(RM) *.o *.mod *.f90

.F.o:
	$(RM) $@ $*.mod


ifeq "$(FORTRILINOS)" "true"
	$(FC) $(CPPFLAGS) $(FFLAGS) -c mpas_ocn_fortrilinos_model_evaluator_mod.F90 $(CPPINCLUDES) $(FCINCLUDES) $(FORTRILINOS_INCLUDES)
	$(FC) $(CPPFLAGS) $(FFLAGS) -c mpas_ocn_fortrilinos_imp_mod.F90 $(CPPINCLUDES) $(FCINCLUDES) $(FORTRILINOS_INCLUDES)
endif


ifeq "$(GEN_F90)" "true"
	$(CPP) $(CPPFLAGS) $(CPPINCLUDES) $< > $*.f90
	$(FC) $(FFLAGS) -c $*.f90 $(FCINCLUDES)
else
	$(FC) $(CPPFLAGS) $(FFLAGS) -c $*.F $(CPPINCLUDES) $(FCINCLUDES) $(FORTRILINOS_INCLUDES)
endif
