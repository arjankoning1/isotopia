# ISOTOPIA Makefile
#
# Examples:
#   make
#   make FC=gfortran FFLAGS="-w -O3 -ffp-contract=off"
#   make FC=ifx FFLAGS="-O3"
#   make clean
#   make check

ifeq ($(origin FC), default)
  FC = gfortran
endif

ifeq ($(origin FFLAGS), undefined)
  ifneq ($(findstring gfortran,$(FC)),)
    FFLAGS = -w -O3 -ffp-contract=off
  else
    FFLAGS =
  endif
endif

ISOTOPIA_DIR := $(abspath ..)
BIN_DIR      := $(ISOTOPIA_DIR)/bin
TARGET       := $(BIN_DIR)/isotopia

ISOTOPIA_LIBS_CHECK := $(abspath $(ISOTOPIA_DIR)/../isotopia.libs)

SOURCES := $(sort $(wildcard *.f) $(wildcard *.f90))

.PHONY: all install clean check

all install: $(TARGET)

$(TARGET): $(SOURCES)
	@echo "Compiler: $(FC)"
	@echo "Flags:    $(FFLAGS)"
	@mkdir -p "$(BIN_DIR)"
	$(FC) $(FFLAGS) -o "$(TARGET)" $(SOURCES)
	@rm -f *.o *.mod
	@echo "Installed: $(TARGET)"

clean:
	rm -f *.o *.mod isotopia "$(TARGET)"

check: all
	@test -d "$(ISOTOPIA_LIBS_CHECK)" || { \
	  echo "ISOTOPIA check error: sibling isotopia.libs directory not found:" >&2; \
	  echo "  $(ISOTOPIA_LIBS_CHECK)" >&2; \
	  exit 1; \
	}
	cd "$(ISOTOPIA_DIR)/samples" && \
	  ISOTOPIA_DIR="$(ISOTOPIA_DIR)" \
	  ISOTOPIA_USER="Arjan Koning" \
	  ./verify
