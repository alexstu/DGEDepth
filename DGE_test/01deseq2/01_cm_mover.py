import sys
import subprocess
import string

def qq(p):
	s = subprocess.check_output(['ls'])
	l = string.split(s, '\n')
	for f in l:
		if (('.txt' in f) == False):
			continue
		with open(f, 'r') as fin:
	 	   data = fin.read().splitlines(True)
		with open(f, 'w') as fout:
	  	  fout.writelines(data[1:])

qq('.')



