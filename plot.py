import matplotlib.pyplot as plt
import numpy as np

############# change #################
L_origin = 50
TIME_STEP = 1e-2
vel = 10.0
long_link_switch = True
#######################################

fname = "forces"
fname2 = "remain_chains"
fdir = "./"+fname+".txt"
fdir2 = "./"+fname2+".txt"
if (long_link_switch):
	fname3 = "add_long_link_info"
	fdir3 = "./"+fname3+".txt"
fout = fname+".png"
fout_peak = fname+"peak.png"
stress = []
stretch1 = []
remain = []
stretch2 = []

# force and stretch
deltaY = TIME_STEP*vel
L = L_origin
i = 0
maxF = 0
maxI = 0
for line in open(fdir, 'r'):
	if (i>13):
		list_of_words = line.split()
		sigma = float(list_of_words[1])
		if (-sigma < 9.99*10**(-5)):
				break
		else:
			stress.append(-sigma)
			L += deltaY
			stretch1.append((L-L_origin)/L_origin)
			if (-sigma > maxF and i >= 100):
				maxF = -sigma
				maxI = i-13
	i += 1

F_initial = stress[0]
stretch_u = round(stretch1[-1],2)
stretch_max_f = round(stretch1[stress.index(maxF)],2)

#zoom in at peak
critical_index = stress.index(maxF)
start_index = int(0.9*critical_index)
end_index = critical_index+int(0.3*(len(stretch2)-critical_index))
stress_d = []
remain_d = []
stretch_d = []


# chain remains
L = L_origin
for line in open(fdir2, 'r'):
	num = int(line.split()[0])
	if (num > 0):
		remain.append(num)
		L += deltaY
		stretch2.append((L-L_origin)/L_origin)




for i in range(start_index,end_index):
	stress_d.append(stress[i])
	remain_d.append(remain[i])
	stretch_d.append(stretch2[i])


if (long_link_switch):
# for long link status:
	force_collector = []
	stretch_collector = []
	i = 0
	for line in open(fdir3, 'r'):
		list_of_words = line.split()
		if (i==0):
			num = int(len(list_of_words)/3)
			for j in range(num):
				force_collector.append([])
				stretch_collector.append([])
		elif (i>len(stretch1)-1):
			break
		else:
			for j in range(num):
				force_val = list_of_words[3*j]
				if (force_val == 0):
					break
				else: 
					force_collector[j].append(force_val)
					stretch_collector[j].append(stretch1[i])
		i += 1


# plot force vs stretch
plt.figure()
plt.xlabel("stretch")
plt.ylabel("Force", color="b")
plt.tick_params(axis="y", labelcolor="b")
plt.plot(stretch1,stress,'b',label='max force is '+str(maxF))
plt.legend(loc=4)

plt.twinx()
plt.ylabel("Remaining Bonds", color="r")
plt.tick_params(axis="y", labelcolor="r")
plt.plot(stretch2,remain,'r',label='initial force is '+str(F_initial))
plt.legend(loc=2)

if (long_link_switch):
	for j in range(num):
		link_name = "long link #"+str(j+1)
		plt.plot(stretch_collector[j],force_collector[j],label=link_name)

ant = 'stretch@Fmax '+str(stretch_max_f)+' ; stretch final: '+str(stretch_u)
plt.title(ant)
plt.savefig(fout)


# plot zoom in
plt.figure()
plt.xlabel("stretch")
plt.ylabel("Force", color="b")
plt.tick_params(axis="y", labelcolor="b")
plt.plot(stretch_d,stress_d,'b')
plt.twinx()
plt.ylabel("Remaining Bonds", color="r")
plt.tick_params(axis="y", labelcolor="r")
plt.plot(stretch_d,remain_d,'r')
plt.savefig(fout_peak)

print(str(F_initial)+'\t'+str(maxF)+'\t'+str(stretch_u)+'\t'+str(stretch_max_f))
