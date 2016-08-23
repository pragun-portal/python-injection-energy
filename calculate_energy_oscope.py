import sys
import matplotlib.pyplot as plt
import numpy as np
from numpy import convolve
import os

V_start = 50					#m/s
V_pierce_max = 238				#m/s
V_pierce_min = 200 				#m/s
T_pierce_start_max = 13 		#ms  (Tstart+13.5ms) time until pierce is started
T_pierce_end_max = 37	 	#ms  (Time from Tpierce start till Tpierce end)
V_follow_max = 130 				#m/s
V_follow_nom = 107.45			#m/s
V_follow_min = 90				#m/s
T_follow_start = 90 			#ms
V_followed = 80					#m/s
V_end = 50						#m/s
T_end = 10						#ms T_followed_end +T_end

#Standard Parameters
pi=3.14156					#pi
Dn = 167e-6					#Nozzel Diameter
An = pi*pow(Dn,2)/4			#Nozzel Area
ro = 999.97					#Water Density
conversion_factor = 0.3		#Kistler CF Setting
moving_average_window = 10


def moving_average (values, window):
    weights = np.repeat(1.0, window)/window
    sma = np.convolve(values, weights, 'valid')
    return sma
 

def raw2velocity(filename):	
	#Read from CSV
	time, value = np.loadtxt(filename, delimiter=',',unpack=True)
	#Data Processing
	#Determine Sensor offset 
	offset = np.mean(value[1:100])
	#print offset
	#Convert from Voltage to Stream Velocity and take out offset. ABS relevant for noise near zero force condition as sqrt doesn't like negitive values. Does not effect critical injection measurements
	value = moving_average(value, moving_average_window)
	value = np.sqrt(np.abs((conversion_factor * (value - offset))/An/ro))
	value = np.append(np.zeros(moving_average_window-1),value)
	return (time, value)


def interpolation(x1,y1,x2,y2,ydes): #returns xdes
	dx = x2-x1
	dy = y2-y1
	xdes = x1 + dx/dy * (ydes-y1)
	return xdes

def find_crossing_time_up(time,value,desired_crossing,start_time,window):
	correspondences=0
	time_sum=0
	for i in range(start_time,len(value)):
		if value[i]>desired_crossing and value[i-1]<desired_crossing:
			if np.mean(value[i:i+window])>desired_crossing and np.mean(value[i-window:i])<desired_crossing:
				correspondences = correspondences + 1
				Crossing_i1 = [i]
				Crossing_i0 = [i-1]
				time_sum = interpolation(time[Crossing_i1],value[Crossing_i1],time[Crossing_i0],value[Crossing_i0],desired_crossing)
				return (time_sum, Crossing_i0[0]) 

def find_crossing_time_down(time,value,desired_crossing,start_time,window):
	correspondences=0
	time_sum=0
	for i in range(start_time,len(value)):
		if value[i]<desired_crossing and value[i-1]>desired_crossing:
			if np.mean(value[i:i+window])<desired_crossing and np.mean(value[i-window:i])>desired_crossing:
				correspondences = correspondences + 1
				Crossing_i1 = [i]
				Crossing_i0 = [i-1]
				time_sum = interpolation(time[Crossing_i1],value[Crossing_i1],time[Crossing_i0],value[Crossing_i0],desired_crossing)
                                return (time_sum, Crossing_i0[0])

w = 10
def processinjection(time, value):
	time50U, index50U = find_crossing_time_up(time,value,V_start,0,w)
	time200U, index200U = find_crossing_time_up(time,value,V_pierce_min,0,w)
	time200D, index200D = find_crossing_time_down(time,value,V_pierce_min,0,w)
	time130D, index130D = find_crossing_time_down(time,value,V_follow_max,0,w)
	time90D, index90D = find_crossing_time_down(time,value,V_follow_min,0,w)
	time50D, index50D = find_crossing_time_down(time,value,50,0,w)

	peaki = np.where(value==np.amax(value))[0]
	peak_mean = value[peaki[0]]# peak_mean = np.mean(value[peaki-10:peaki+10])
	points_x = [time50U,time200U,time200D,time130D,time90D,time50D,time[peaki]]
	indices = [index50U, index200U, index200D, index130D, index90D, index50D, peaki]
	points_y = [V_start,V_pierce_min,V_pierce_min,V_follow_max,V_follow_min,V_end,peak_mean]
	return (points_x,points_y, indices)

def graph(time, value, points_x, points_y):	
	#Plotting
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1, axisbg='white')
	plt.plot(time,value)
	plt.plot([points_x[3],points_x[3],points_x[4],points_x[4],points_x[3]],[V_follow_min,V_follow_max,V_follow_max,V_follow_min,V_follow_min],'g.-')	#Plot Follow Through Window
	plt.plot([points_x[1],points_x[1],points_x[2],points_x[2],points_x[1]],[V_pierce_min,V_pierce_max,V_pierce_max,V_pierce_min,V_pierce_min],'r.-')		#Plot Pierec Window
	#plt.plot(points_x,points_y,'ro')
	plt.title('Injection Crtical Values')
	plt.ylabel('Stream Velocity (m/s)')
	plt.xlabel('Time (s)')
	plt.grid()
	plt.show()

def report(time,value,points_x,points_y,indices, reportfile):
	TStartPierce = 1000*(points_x[1]-points_x[0])
	pierce_duration = 1000*(points_x[2]-points_x[1])
	TStartFollow = 1000*(points_x[3]-points_x[0])
	follow_fall_time = 1000*(points_x[3]-points_x[2])
	follow_duration = 1000*(points_x[4]-points_x[3])
	end_fall_time = 1000*(points_x[5]-points_x[4])

	energy = calculate_energy(indices[0],indices[5],time,value)
	pierce_energy = calculate_energy(indices[0],indices[3],time,value)
        ft_energy = calculate_energy(indices[3],indices[5],time,value)
        

	Pierce_max = np.amax(value)
	Pierce_avg = np.mean(value[np.where(value>V_pierce_min)])
	Follow_avg = np.mean(value[np.logical_and(value<130,value>80)])
	print '\n\n\n\n'
	print '################### Phase Timing #####################\n'
	print 'TStart_Pierce = %0.2f ms'%TStartPierce
	print 'Expecting Less than %0.2f ms\n'%T_pierce_start_max
	print 'Pierce Duration  = %0.2f ms'%pierce_duration
	print 'Maximum %0.2f ms\n'%T_pierce_end_max
	print 'TStart_Follow = %0.2f ms'%TStartFollow
	print "Expected Less than %0.2f ms \n"%T_follow_start
	print 'Follow Transition Time = %0.2f ms'%follow_fall_time
	print 'No Spec Specified\n'

	print 'Follow Time  = %0.2f ms'%follow_duration
	print 'No Spec Specified\n'
	print 'Stop Time = %0.2f ms'%end_fall_time
	print 'Expecting Less than %0.2f ms\n'%T_end
	print '############## Velocity Measurements ################\n'
	print 'Average Pierce Velocity = %0.2f m/s'%Pierce_avg
	print 'Maximum Pierce Velocity = %0.2f m/s'%Pierce_max
	print 'Expected value range between %0.2f m/s and %0.2f m/s\n'%(V_pierce_min,V_pierce_max)
	print 'Average Follow Velocity = %0.2f m/s'%Follow_avg
	print 'Expected value range between %0.2f m/s and %0.2f m/s\n'%(V_follow_min,V_follow_max)
	print "Energy in Jet: %0.2fJ\n"%(energy,)
	print "Energy in pierce: %0.2fJ\n"%(pierce_energy,)
	print "Energy in follow-through: %0.2fJ\n"%(ft_energy,)

	reportfile.write("%0.2f\t%0.2f\t%0.2f"%(energy,pierce_energy,ft_energy))
	print '#################### Complete #######################\n'


def calculate_energy(starting_index,ending_index,time,velocity):
    energy = 0
    for i in range(starting_index,ending_index):
        delta_time = time[i+1] - time[i]
        dvelocity = (velocity[i+1]+velocity[i])/2.0
        delta_energy = An*ro*dvelocity*dvelocity*dvelocity*delta_time*(0.5)
        energy += delta_energy
    return energy



reportfile = open("energy_report.txt",'w')

filenames = os.listdir(".")
for filename in filenames:
    if filename.endswith(".csv"):
        print "Processing File: %s"%(filename,)
        time_stamps, stream_velocity = raw2velocity(filename)
        points_x, points_y, indices = processinjection(time_stamps,stream_velocity)
        print points_x,points_y
        graph(time_stamps,stream_velocity,points_x,points_y)

        reportfile.write("%s\t"%(filename,))
        report(time_stamps,stream_velocity,points_x,points_y,indices, reportfile)
        reportfile.write("\n")

reportfile.close()


