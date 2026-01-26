# CC  ?= clang
CC=/usr/bin/cc
# CXX ?= clang++
LD  = $(CC)

# --- Cross-platform OpenSSL detection ---
UNAME := $(shell uname -s)
ifeq ($(UNAME),Darwin)
  # macOS - try to detect OpenSSL via brew
  OPENSSL_PREFIX := $(shell brew --prefix openssl@3 2>/dev/null)
  ifeq ($(OPENSSL_PREFIX),)
    OPENSSL_PREFIX := $(shell brew --prefix openssl 2>/dev/null)
  endif
  ifeq ($(OPENSSL_PREFIX),)
    $(error "OpenSSL not found on macOS. Please install via 'brew install openssl@3' and try again.")
  endif
else
  # Linux and other systems - use system OpenSSL
  OPENSSL_PREFIX := $(shell pkg-config --variable=prefix openssl 2>/dev/null)
  ifeq ($(OPENSSL_PREFIX),)
    # Fallback to common system paths
    ifneq ($(wildcard /usr/include/openssl),)
      OPENSSL_PREFIX := /usr
    else ifneq ($(wildcard /usr/local/include/openssl),)
      OPENSSL_PREFIX := /usr/local
    else
      $(error "OpenSSL not found on Linux. Please install libssl-dev (Ubuntu/Debian) or openssl-devel (RedHat/CentOS)")
    endif
  endif
endif

ifndef PROJ
PROJ = ref

endif


SRC_DIR  = ./src
SRC_EXT_DIRS = ./src/ref
UTIL_DIR = ./utils
NIST_DIR = ./utils/nistkat
MASKING_DIR = ./masking
TEST_DIR = ./test

CFLAGS +=  -w -Wall -Wextra -Wpedantic -Wredundant-decls \
  -Wshadow -Wpointer-arith -O3 -Werror -fomit-frame-pointer


ifdef PARAM
ifeq ($(PARAM),3)
CFLAGS    += -D_OV256_112_44
else ifeq ($(PARAM),4)
CFLAGS    += -D_OV256_184_72
else ifeq ($(PARAM),5)
CFLAGS    += -D_OV256_244_96
else
CFLAGS    += -D_OV16_160_64
endif
else
PARAM=3
endif


ifdef VARIANT
ifeq ($(VARIANT),2)
CFLAGS += -D_OV_PKC
else ifeq ($(VARIANT),3)
CFLAGS += -D_OV_PKC_SKC
else ifeq ($(VARIANT),4)
CFLAGS += -D_OV_CLASSIC -D_4ROUND_AES_
else ifeq ($(VARIANT),5)
CFLAGS += -D_OV_PKC -D_4ROUND_AES_
else ifeq ($(VARIANT),6)
CFLAGS += -D_OV_PKC_SKC -D_4ROUND_AES_
else
CFLAGS += -D_OV_CLASSIC
endif
else
VARIANT=1
endif


INCPATH := \
	-I.\
	-I$(OPENSSL_PREFIX)/include \
	-I/usr/local/include -I/opt/local/include -I/usr/include \
	-I$(SRC_DIR) -I$(SRC_EXT_DIRS) -I$(UTIL_DIR) -I$(MASKING_DIR)

# Link flags and libraries
ifeq ($(UNAME),Darwin)
  # macOS - explicitly link to brew-installed OpenSSL
  LDFLAGS += -L$(OPENSSL_PREFIX)/lib
else
  # Linux - use pkg-config if available, otherwise default paths
  OPENSSL_LIBS := $(shell pkg-config --libs openssl 2>/dev/null)
  ifneq ($(OPENSSL_LIBS),)
    LIBS = $(OPENSSL_LIBS)
  else
    LDFLAGS += -L$(OPENSSL_PREFIX)/lib
    LIBS = -lssl -lcrypto
  endif
endif
# Fallback if LIBS not set above
ifeq ($(LIBS),)
  LIBS = -lssl -lcrypto
endif

EXE = bench_test sign_test

# Macro
ORDER=1
MACRO = -D MASKING_ORDER=$(ORDER)

# Source
SRC := $(wildcard $(SRC_EXT_DIRS)/*.c) \
	$(wildcard $(SRC_DIR)/*.c) \
	$(wildcard $(UTIL_DIR)/*.c) \
	$(wildcard $(NIST_DIR)/*.c) \
	$(wildcard $(MASKING_DIR)/*.c) \

TEST_SRC := $(wildcard $(TEST_DIR)/*.c)


.PHONY: all clean
all: $(EXE)

bench_test: $(SRC) $(TEST_DIR)/bench_test.c
	$(CC) $(CFLAGS) $(MACRO) $(INCPATH) $(LDFLAGS) -o $@ $^ $(LIBS)

sign_test: $(SRC) $(TEST_DIR)/sign_test.c
	$(CC) $(CFLAGS) $(MACRO) $(INCPATH) $(LDFLAGS) -o $@ $^ $(LIBS)


clean:
	rm -f $(EXE)
