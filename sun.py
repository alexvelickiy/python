#!/usr/bin/python -tt

import math
import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import sys

def normalize(vec) :
	return vec / la.norm(vec)

def get_person_pos(person_angle_pos_rad) :

	z = earth_r * math.sin(person_angle_pos_rad[1])

	pr_len = math.cos(person_angle_pos_rad[1])

	x = earth_r * math.cos(person_angle_pos_rad[0]) * pr_len
	y = earth_r * math.sin(person_angle_pos_rad[0]) * pr_len

	return np.array([x, y, z])

def get_longitude_tang(person_angle_pos_rad) :
	z = earth_r * math.cos(person_angle_pos_rad[1])

	pr_len = -math.sin(person_angle_pos_rad[1])

	x = earth_r * math.cos(person_angle_pos_rad[0]) * pr_len
	y = earth_r * math.sin(person_angle_pos_rad[0]) * pr_len

	return normalize(np.array([x, y, z]))

def get_latitude_tang(person_angle_pos_rad) :
	z = 0

	pr_len = math.cos(person_angle_pos_rad[1])

	x = - earth_r * math.sin(person_angle_pos_rad[0]) * pr_len
	y = earth_r * math.cos(person_angle_pos_rad[0]) * pr_len

	return normalize(np.array([x, y, z]))

def get_pr_coords(sky_pos_pr, long_tang, latid_tang) :
	mat = np.array([latid_tang, long_tang]).T

	sol = la.lstsq(mat, sky_pos_pr)
	
	return sol[0]

def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis/math.sqrt(np.dot(axis, axis))
    a = math.cos(theta/2.0)
    b, c, d = -axis*math.sin(theta/2.0)
    aa, bb, cc, dd = a*a, b*b, c*c, d*d
    bc, ad, ac, ab, bd, cd = b*c, a*d, a*c, a*b, b*d, c*d
    return np.array([[aa+bb-cc-dd, 2*(bc+ad), 2*(bd-ac)],
                     [2*(bc-ad), aa+cc-bb-dd, 2*(cd+ab)],
                     [2*(bd+ac), 2*(cd-ab), aa+dd-bb-cc]])


def rotate_vec(vec, theta, axis) :
	return np.dot(rotation_matrix(axis,theta), vec) 

def plot_sun(person_lat = 50):
	circle_coords = np.array([]).reshape(0,2)
	# drawing vision circle
	for ang in range(360) :
		circle_coords = np.vstack((circle_coords, np.array([math.cos(ang), math.sin(ang)])))
	plt.scatter(circle_coords[:, 0], circle_coords[:, 1], color='b', s=1)



	sun_durations = []



	for eartg_ang in range(0, 180, 30) :

		sun_pos = np.array([sun_distance, 0, sun_h_max * math.cos(math.radians(eartg_ang))])
		sun_pr_coords = np.array([]).reshape(0,3)

		for ang in range(0, 360, 1) :

			person_angle_pos = np.array([ang, person_lat])

			person_angle_pos_rad = np.radians(person_angle_pos)


			person_pos = get_person_pos(person_angle_pos_rad)
			sky_pos = normalize(sun_pos - person_pos) * sky_r
			person_pos_normal = person_pos / earth_r

			sky_pos_normal_dot_prod = np.dot(sky_pos, person_pos_normal)
			if sky_pos_normal_dot_prod < 0 :
				pass
			else :
				sky_pos_pr = sky_pos - sky_pos_normal_dot_prod * person_pos_normal


				long_tang = get_longitude_tang(person_angle_pos_rad)
				latid_tang = get_latitude_tang(person_angle_pos_rad)

				pr_coords = get_pr_coords(sky_pos_pr, long_tang, latid_tang)
				# pr_coords2 = rotate_vec(sky_pos_pr, math.radians(180 + 90 - person_angle_pos[0]), [0, 0, 1])
				# pr_coords2 = rotate_vec(pr_coords2, math.radians(90 - person_angle_pos[1]), [1, 0, 0])
				# print(pr_coords2)

				sun_pr_coords = np.vstack((sun_pr_coords, np.append(pr_coords, ang)))

		sun_pr_coords = sun_pr_coords.T
		plt.scatter(sun_pr_coords[0], sun_pr_coords[1])
		if len(sun_pr_coords[2]) :
			min_ang = np.min(sun_pr_coords[2])
			max_ang = np.max(sun_pr_coords[2])
			sun_duration = (max_ang - min_ang) / 360 * 24
			sun_durations.append(str(round(sun_duration, 1)) + " hours of sun")
		else :
			sun_durations.append('-')
	plt.legend(['-'] + sun_durations)
	plt.show()



sun_distance = -149.6 * 10**6 / 6371

# 23.5 is the Earth tilt angle
sun_h_max = abs(math.sin(math.radians(23.5)) * sun_distance)


earth_r = 1
sky_r = 1


