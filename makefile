include makefile.preample

all:
	make -C src

install:
	$(MOVE) "src/spec$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/clean$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/sclean$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/lclean$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/window$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/lowpass$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/highpass$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/bandpass$(EXT)" "bin/$(ARCH)/"
	$(MOVE) "src/bandstop$(EXT)" "bin/$(ARCH)/"
	
clean:
	make --directory=src clear
