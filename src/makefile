include ../makefile.preample

.PHONY: all clear

%.o: %.f90
	$(F95) $(FFLAGS) -c $*.f90 $(FLINK)
%.o: %.f
	$(F95) $(FFLAGS) -c $*.f $(FLINK)
	
# Files to compile:
commonfiles = strings.o dpsort.o timeseries.o
findmaxfiles = d1mach.o fmin.o FindSpectrumMax.o 

all: spec clean lclean sclean window lowpass highpass bandpass bandstop
ifeq ($(F95),g95)
	@echo "========================================"
	@echo " Please note:"
	@echo " g95 does not support OpenMP directives"
	@echo "========================================"
endif

d1mach.o: d1mach.f makefile
	$(FC) $(FFLAGS) -c d1mach.f $(LBLAS)

fmin.o: fmin.f makefile
	$(FC) $(FFLAGS) -c fmin.f

spec: $(commonfiles) spec.o makefile
	$(F95) $(FFLAGS) -o spec $(commonfiles) spec.o $(FLINK)

clean: $(commonfiles) $(findmaxfiles) clean.o makefile
	$(F95) $(FFLAGS) -o clean $(commonfiles) $(findmaxfiles) clean.o $(FLINK)

sclean: $(commonfiles) $(findmaxfiles) sclean.o makefile
	$(F95) $(FFLAGS) -o sclean $(commonfiles) $(findmaxfiles) sclean.o $(FLINK)

lclean: $(commonfiles) $(findmaxfiles) lclean.o makefile
	$(F95) $(FFLAGS) -o lclean $(commonfiles) $(findmaxfiles) lclean.o $(FLINK)

window: $(commonfiles) integrate.o window.o makefile
	$(F95) $(FFLAGS) -o window $(commonfiles) integrate.o window.o $(FLINK)

lowpass: $(commonfiles) lowpass.o makefile
	$(F95) $(FFLAGS) -o lowpass $(commonfiles) lowpass.o $(FLINK)

highpass: $(commonfiles) highpass.o makefile
	$(F95) $(FFLAGS) -o highpass $(commonfiles) highpass.o $(FLINK)

bandpass: $(commonfiles) bandpass.o makefile
	$(F95) $(FFLAGS) -o bandpass $(commonfiles) bandpass.o $(FLINK)

bandstop: $(commonfiles) bandstop.o makefile
	$(F95) $(FFLAGS) -o bandstop $(commonfiles) bandstop.o $(FLINK)

# For cleaning up:
clear:
	$(DELETE) *.mod *.o *.obj spec$(EXT) clean$(EXT) sclean$(EXT) lclean$(EXT) window$(EXT) lowpass$(EXT) highpass$(EXT) bandpass$(EXT) bandstop$(EXT) converttime$(EXT) core* *~
