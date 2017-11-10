import matplotlib.pyplot as plt
import numpy as np

############# change #################
fname = "0.166667_false"
fname2 = "0.166667_remain_chains"
L_origin = 50
TIME_STEP = 1e-2
vel = 10.0
#######################################


fdir = "./"+fname+".txt"
fdir2 = "./"+fname2+".txt"
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
			if (-sigma > maxF):
				maxF = -sigma
				maxI = i-13
	i += 1

# chain remains
L = L_origin
for line in open(fdir2, 'r'):
	num = int(line.split()[0])
	if (num > 0):
		remain.append(num)
		L += deltaY
		stretch2.append((L-L_origin)/L_origin)


# if (len(stress)<len(stretch)):
# 	for i in range(len(stretch)-len(stress)):
# 		stress.append(0)

# print((stress))
F_initial = stress[0]
stretch_u = round(stretch1[-1],2)
stretch_max_f = round(stretch1[stress.index(maxF)],2)

#zoom in at peak
critical_index = stress.index(maxF)
# start_index = critical_index-int(0.5/TIME_STEP)
# end_index = critical_index+10
start_index = int(0.9*critical_index)
# end_index = int((critical_index+len(stretch1))/3)
end_index = critical_index+int(0.3*(len(stretch2)-critical_index))
stress_d = []
remain_d = []
stretch_d = []

# print(len(stress))
# print(start_index)
# print(end_index)
for i in range(start_index,end_index):
	stress_d.append(stress[i])
	remain_d.append(remain[i])
	stretch_d.append(stretch2[i])

plt.figure()
plt.xlabel("stretch")
plt.ylabel("Force", color="b")
plt.tick_params(axis="y", labelcolor="b")
plt.plot(stretch1,stress,'b',label='max force is '+str(maxF))
plt.legend(loc=1)

plt.twinx()
plt.ylabel("Remaining Bonds", color="r")
plt.tick_params(axis="y", labelcolor="r")
plt.plot(stretch2,remain,'r',label='initial force is '+str(F_initial))
plt.legend(loc=2)

# ant = 'max top plate force is '+str(maxF)
# plt.annotate(ant, xy=(maxI, maxF), xytext=(maxI-100, maxF-0.1),
#             arrowprops=dict(facecolor='black'),
#             )
ant = 'stretch@Fmax '+str(stretch_max_f)+' ; stretch final: '+str(stretch_u)
# plt.annotate(ant, xy=(maxI, maxF), xytext=(maxI, maxF),
#             arrowprops=dict(facecolor='black'),
#             )
plt.title(ant)
plt.savefig(fout)
# plt.show()

#zoom in
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
