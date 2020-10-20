import os
from pathlib import Path
import csv
import re

def load_steps(file_name):
	steps = []
	data = open(file_name, "r")

	for line in data:
		data_array = re.split(r'\s+', line)
		if (len(data_array) >= 10):
			try:
				#print(float(data_array[7]))
				steps.append(float(data_array[7]))
			except ValueError:
				pass
				#print("Not a float")
	return steps

steps_filename = "tensile_dog_bone_1cpu.sta"
load_steps = load_steps(steps_filename)
print("Load steps:", load_steps)