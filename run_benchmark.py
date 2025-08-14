from os import system, popen
import sys
PATH = "./bench_res.txt"
open(PATH, 'w').close() # clear file
cur = ""
res = ""

VERBOSE_COMPILE = True
REDIRECT = ""
MAX_ORDER = 6

if not VERBOSE_COMPILE:
	REDIRECT=">/dev/null"
with open(PATH,'a') as f:
  for PARAM in range(3, 4):
    if PARAM == 3:
        variant = "UOV-I"
    elif PARAM == 4:
        variant = "UOV-III"
    else:  # PARAM == 5
        variant = "UOV-V"
    var = f"Compiling for {variant}"
    print(var)
    f.write(var + "\n")
    for i in range(1, MAX_ORDER+1):
      print("Compiling for masking of order", i)
      system("make clean > /dev/null && make bench_test" + " PARAM=" + str(PARAM) + " ORDER=" + str(i)+REDIRECT)
      print("Running tests...", end ='')
      sys.stdout.flush()
      cur = popen("./bench_test").read()
      print("Writing to "+PATH+" ...")
      f.write(cur)
      res += cur
      print(cur)
      print("Done.")