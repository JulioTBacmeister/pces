##################################
#  $Id: GNUmakefile,v 1.24 2008/11/30 02:42:01 bacmj Exp $
##################################

#EXEDIR=.
#EXEN=clstr.x

RNTG2=cluster

VPATH=share

COMMON_FILES =                  \
	numerical_utilities.F90 \
	fish.F90                \
	comf.F90                \
	pois3d.F90              \
	besselfuncs.F90         \
	ce1_types.F90           \
	ce1_consts.F90          \
	ppm.F90                 \
	reader.F90              \
	get_driver_data.F90     \
	ce1_diags.F90           \
	ce1_inform.F90          \
	ce1_utils.F90           \
	ce1_gcm_coupler.F90     \
	mp_morr_two_moment.F90 \
	ce1_uphys.F90           \
	ce1_create_destroy.F90  \
	ce1_precipitation.F90   \
	ce1_interactions.F90    \
	ce1_prsolvph00.F90        \
	ce1_prsolvph2.F90        \
	ce1_pr3d.F90        \
	ce1_vdiff.F90        \
	ce1_dynamics.F90        \
	ce1_adjustments.F90  


 SRCS= $(COMMON_FILES) \
	ce1.F90 \
	pce_comp.F90 \
	single_clstr_drv.F90 \

 TEST_SRCS= $(COMMON_FILES) \
	test.F90 \

 SRCS_F77= \
	fftpack.f \
	xerbla.f \
	sgtsv.f \

# Make this arch dependent
#COMPILER=f90
#COMPILER=ifort
#COMPILER=lf95
#COMPILER=pgf95
COMPILER=gfortran


OBJS= ${SRCS:.F90=.o}
#OBJSUFF= ${SRCS:.F90=.o}
#OBJS=$(patsubst share/%.F90,%.o,$(SRCS))
OBJS_F77= ${SRCS_F77:.f=.o}
TEST_OBJS= ${TEST_SRCS:.F90=.o}
DRYTEST_OBJS= ${DRYTEST_SRCS:.F90=.o}

F77_FLAGS=

STD_FLAGS=-c

ifeq ($(COMPILER),gfortran)
  #STD_FLAGS := $(STD_FLAGS) -O3 -fdollar-ok -DGFORTRAN
  STD_FLAGS := $(STD_FLAGS) -g -fdollar-ok -fbounds-check -fdefault-real-8 -DGFORTRAN
endif
ifeq ($(COMPILER),pgf95)
  STD_FLAGS := $(STD_FLAGS) -g
endif

F90_FLAGS=$(PFLAGS) $(STD_FLAGS) $(USER_FDEFS)

CPP_FLAGS=-P -DALPHA_MACH


$(RNTG2) : $(OBJS) $(OBJS_F77)
	$(COMPILER) -O1 $(OBJS) $(OBJS_F77) -o clstr.x

#-----------------------------------------------------
# These files (RHS) are used by almost everyone.
# So, safer just to rebuild everything if they change.
#-----------------------------------------------------
$(OBJS) : ce1_consts.F90 ce1_types.F90
#$(CLSTROBJS) : ce1_consts.F90 ce1_types.F90

#$(OBJS) : %.o: %.F90
#XXgoldyXX: Changed this into an implicit rule to allow VPATH to work
%.o: %.F90
	$(COMPILER) $(F90_FLAGS) $<
#	$(COMPILER) $(F90_FLAGS) $*.F90



#$(OBJS_F77) : %.o: %.f
#XXgoldyXX: Changed this into an implicit rule to allow VPATH to work
%.o: %.f
	$(COMPILER) $(F77_FLAGS) $(F90_FLAGS) $<
#	$(COMPILER) $(F77_FLAGS) $(F90_FLAGS) $*.f


objtest: GNUmakefile
	echo $(OBJS)


testx: $
	$(COMPILER) -O1 $(TEST_SRCS) $(SRCS_F77) -o tst.x

toyp:
	$(COMPILER) -O3 -r8 -c share/besselfuncs.F90
	$(COMPILER) -O3 -r8 -c fish.F
	$(COMPILER) -O3 -r8 -c fftpack.f
	$(COMPILER) -O3 -r8 -c comf.f
	$(COMPILER) -O3 -r8 -c pois3d.f
	$(COMPILER) -O3 -r8 -c toyp3.F90
	$(COMPILER) -O3 -r8 besselfuncs.o fish.o comf.o fftpack.o pois3d.o toyp3.o -o toyp.x

wojp:
	ln -sf /project/convection/juliob/wojtek/fort.17 wojtek.dat
	$(COMPILER) -O3 -c fish.F90
	$(COMPILER) -O3 -c fftpack.f
	$(COMPILER) -O3 -c comf.F90
	$(COMPILER) -O3 -c pois3d.F90
	$(COMPILER) -O3 -c solve_wojtek_pr.F90
	$(COMPILER) -O3 fish.o comf.o fftpack.o pois3d.o solve_wojtek_pr.o -o wojp.x

testdiff:
	$(COMPILER) -c -g xerbla.f
	$(COMPILER) -c -g sgtsv.f
	$(COMPILER) -c -g share/numerical_utilities.F90
	$(COMPILER) -c -g share/ce1_types.F90
	$(COMPILER) -c -g share/ce1_consts.F90
	$(COMPILER) -c -g share/ce1_inform.F90
	$(COMPILER) -c -g share/ce1_utils.F90
	$(COMPILER) -c -g share/ce1_vdiff.F90
	$(COMPILER) -c -g test_diff.F90
	$(COMPILER) xerbla.o sgtsv.o numerical_utilities.o ce1_types.o ce1_inform.o ce1_utils.o ce1_vdiff.o test_diff.o -o test_diff.x

cleancode : 
	/bin/rm -f *.o *.a *.x *.mod *.dvi *.ps *.pdf *.ofl *.aux *.log *~

clean : 
	/bin/rm -f *.o *.a *.x *.mod *.dvi *.ps *.pdf *.ofl *.aux *.log *~ fort.* share/*~ idlpros/*~

