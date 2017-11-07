####### Change this #########
fname = "polymer_10_10_origin"
num_add = 4
#########################
import numpy as np
import random
fdir = "./"+fname+".msh"
# fout = "./"+fname+"_rand_"+str(num_add)+".msh"
fout = "test4.msh"

print("start")
# read node, element number:
num_nodes = 0
num_elements = 0
ct = 0
found_element = False
with open(fout,'w') as new_file:
	with open(fdir, 'r') as origin_file:
		lines = origin_file.readlines()

		for i in range(len(lines)):
			if (found_element):
				line = str(int(num_elements+num_add))+"\n"
				found_element = False
			else:
				line = lines[i]
				if ("$EndElements" in line):
					break
			if (str(num_elements) == line.split()[0]):
				ct += 1
				if (ct == 2):
					break;
			if ("$Nodes" in line):
				num_nodes = int(lines[i+1])
			if ("$Elements" in line):
				num_elements = int(lines[i+1])
				found_element = True
			new_file.write(line)

	for i in range(num_add):
		a = 0
		b = 0
		while (a == b):
			a = int(random.randrange(1, num_nodes, 1))
			b = int(random.randrange(1, num_nodes, 1))
		add_line = str(num_elements+i+1)+" "+str(1)+" 2 1 1 "
		add_line+= str(a)+" "+str(b)
		add_line+= '\n'
		new_file.write(add_line)
	new_file.write("$EndElements")

	

new_file.close()
origin_file.close()
