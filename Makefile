LDLIBS += -lm
CFLAGS += -Wall -Wextra -O3
PYTHON := python2

tests: tests-plonka tests-dct

#
# generate test-dct-{4,8,16,...} rules
#
dct-%: dct-%.o

DCT_N = 2 3 4 5
define DEFINE_C_DCT_TEST
$(eval DIM = $(shell echo $$((1 << $(1)))))
dct$(DIM).c: gen_c.py template.c
	@echo generate dct$(DIM).c
	@$(PYTHON) gen_c.py $(1)
DCT_SOURCES += dct$(DIM).c
DCT_BINS += dct$(DIM)
test-dct-$(DIM): dct$(DIM)
	@echo test-dct-$(DIM)
	@./dct$(DIM)
DCT_TESTS += test-dct-$(DIM)
endef
$(foreach N,$(DCT_N),$(eval $(call DEFINE_C_DCT_TEST,$(N))))
tests-dct: $(DCT_TESTS)

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
	$(RM) $(DCT_SOURCES)
distclean: clean
	$(RM) $(DCT_BINS)
