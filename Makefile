LDLIBS += -lm
CFLAGS += -Wall -Wextra -O3
PYTHON := python2

all: tests-plonka tests-fdct

#
# generate test-fdct-{4,8,16,...} rules
#
fdct%: fdct%.o

FDCT_N = 2 3 4
define DEFINE_C_FDCT_TEST
$(eval DIM = $(shell echo $$((1 << $(1)))))
fdct$(DIM).c: gen_c.py template.c
	@echo generate fdct$(DIM).c
	@$(PYTHON) gen_c.py $(1)
FDCT_SOURCES += fdct$(DIM).c
test-fdct-$(DIM): fdct$(DIM)
	@echo test-fdct-$(DIM)
	@./fdct$(DIM)
FDCT_TESTS += test-fdct-$(DIM)
endef
$(foreach N,$(FDCT_N),$(eval $(call DEFINE_C_FDCT_TEST,$(N))))
tests-fdct: $(FDCT_TESTS)

#
# generate test-plonka-{cosI,cosII,...}-{2,4,8,...} rules
#
TFMS = cosI cosII cosIII cosIV sinI
TFMS_BITS = 1 2 3 4 5 6 7 8
define DEFINE_PLONKA_TEST
test-plonka-$(1)-$(2):
	@echo test-plonka-$(1)-$(shell echo $$((1 << $(2))))
	@$(PYTHON) plonka.py $(1) $(2)
PLONKA_TESTS += test-plonka-$(1)-$(2)
endef
$(foreach BITS,$(TFMS_BITS),\
    $(foreach T,$(TFMS),\
        $(eval $(call DEFINE_PLONKA_TEST,$(T),$(BITS)))))
tests-plonka: $(PLONKA_TESTS)

clean:
	$(RM) $(FDCT_SOURCES)
distclean: clean
	$(RM) $(FDCT_TESTS)
