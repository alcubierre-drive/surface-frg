# Surface FRG for a cubic lattice

This program is based on the
[divERGe](https://doi.org/10.21468/SciPostPhysCodeb.26) library and implements
the renormalization of four-point interactions on the surface of a Hubbard-like
model on the cubic lattice.

## Installation
To build the program, you have to install the
[divERGe](https://doi.org/10.21468/SciPostPhysCodeb.26) library (e.g., from
[git](https://git.rwth-aachen.de/frg/diverge)) and point the Makefile to the
right location. Practically, this means creating a ``Makefile.local`` file with
the following contents:

    INCLUDES += -I$(DIVERGE_LOCATION)/src
    LDFLAGS += -L$(DIVERGE_LOCATION) -Wl,-rpath=$(DIVERGE_LOCATION)

Then type ``make`` and you should see a program called ``surf`` in the repo's
root directory.

## Usage
The options can be printed using ``./surf -h``. Each invocation of the program
corresponds to a single FRG simulation. Typically, you'd wrap the binary with a
shell script for phase diagram scans that you can submit on a cluster (or run
locally for TUFRG).

## Theory
The surface Green's function is calculated by the method presented in
[Lopez Sancho's 1985 paper](https://doi.org/10.1088/0305-4608/15/4/009).
Everything else follows the usual FRG derivation, as the theory (and the divERGe
implementation) is fully Green's function based.