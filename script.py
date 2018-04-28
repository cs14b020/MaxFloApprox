import os
os.system("make clean")
os.system("make")
nos = [100, 1000]
maxvars = [0.001, 0.02, 0.07, 0.1, 0.25, 0.5, 1.0]
for no in nos :
	for maxvar in maxvars :
		for i in range(10):
			s = "./main " + str(no) + " " + str(i) + " " + str(maxvar)
			os.system(s)