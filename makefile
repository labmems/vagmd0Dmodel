# Project name
PROJNAME=vagmd0Dmodel

# Definition of the compiler to be used
CC=gcc

# Binary path
BINPATH=./bin/$(PROJNAME)

# Definition of C files
CFILES=$(wildcard ./src/*.c) $(wildcard ./src/*/*.c)

# Definition of the object files
OBJ=$(subst .c,.o,$(CFILES))

# Get help on how to run the binary
help:
	@ ./bin/$(PROJNAME) -help | less

# Run the binary
run:
	@ time ./bin/$(PROJNAME) \
	-membrane_area 25.92 \
	-vacuum_pressure -81325.0 \
	-number_channels 4 \
	-entry_temperature_feed 60.0 \
	-entry_temperature_cool 25.0 \
	-entry_salinity_feed 0.035 \
	-entry_salinity_cool 0.035 \
	-feed_mass_flow_rate 0.08333333333333333333333333 \
	-cool_mass_flow_rate 0.08333333333333333333333333
	@ cat ./results/report.csv

# Build the binary
build: binfolder $(BINPATH)
	@ mkdir -p results

$(BINPATH): $(OBJ)
	@ $(LINK.C) -s -O3 -o $@ $^ $(LDLIBS)
	@ rm -rf ./src/*.o
	@ rm -rf ./src/*/*.o

# Create bin folder
binfolder:
	@ mkdir -p bin

# Clean all files in the bin and results folders
cleanall:
	@ rm -rf ./bin ./graphs
	@ rm -rf ./*.o ./src/*.o ./src/*/*.o

# Inclusion of PETSc config information
include ${PETSC_DIR}/lib/petsc/conf/variables
include ${PETSC_DIR}/lib/petsc/conf/rules
include ${PETSC_DIR}/lib/petsc/conf/test