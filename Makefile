LDLIBS += -lm
CFLAGS += -Wall -Wextra -O3
PYTHON := python2

tests: tests-plonka tests-dct

#
# generate test-dct-{4,8,16,...} rules
#
DCT_N = 2 3 4 5
define DEFINE_C_DCT_TEST
$(eval DIM = $(shell echo $$((1 << $(1)))))
$(eval TEST_TOOL = dct$(DIM))
$(eval C_FILE = $(TEST_TOOL).c)
$(C_FILE): gen_c.py template.c
	@echo generate $(C_FILE)
	@$(PYTHON) gen_c.py $(1)
$(TEST_TOOL): $(TEST_TOOL).o
test-$(TEST_TOOL): $(TEST_TOOL)
	@echo $$@
	@./$(TEST_TOOL)
DCT_BINS += $(TEST_TOOL)
DCT_TESTS += test-$(TEST_TOOL)
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
	$(RM) $(addsuffix .c,$(DCT_BINS))
	$(RM) $(addsuffix .o,$(DCT_BINS))
	$(RM) $(DCT_BINS)

.PHONY: tests tests-plonka tests-dct $(PLONKA_TESTS) $(DCT_TESTS) clean
