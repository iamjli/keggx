#!/usr/bin/env python3

import numpy as np

import matplotlib
import matplotlib.pyplot as plt


class Node: 

	def __init__(self, attribs): 

		for key,val in attribs.items():
			setattr(self, key, val)

	@property
	def anchor(self): 
		# Bottom left corner
		if self.shape == 'circle': 
			return np.array([self.x, self.y])
		else: 
			return np.array([self.x-self.width/2, self.y-self.height/2])

	@property
	def center(self): 
		x,y = self.anchor
		return np.array([x+self.width/2, y+self.height/2])

	@property
	def left(self):
		x,y = self.anchor
		return np.array([x, y+self.height/2])

	@property
	def right(self):
		x,y = self.anchor
		return np.array([x+self.width, y+self.height/2])

	@property
	def top(self):
		x,y = self.anchor
		return np.array([x+self.width/2, y+self.height])

	@property
	def bottom(self):
		x,y = self.anchor
		return np.array([x+self.width/2, y])
	



def shortest_arrow(source, target): 

	source_pos, target_pos = None, None
	min_distance = np.inf

	for pos1 in [source.left, source.right, source.top, source.bottom]: 
		for pos2 in [target.left, target.right, target.top, target.bottom]: 
			distance = ((pos1 - pos2)**2).sum()
			if distance < min_distance: 
				source_pos, target_pos = pos1, pos2
				min_distance = distance

	return source_pos, target_pos


def set_grid(xlim, ylim, scale=1): 

	x_min, x_max = xlim
	y_min, y_max = ylim
	x_pad = (x_max - x_min) * 0.1
	y_pad = (y_max - y_min) * 0.1

	fig, ax = plt.subplots(1, 1, figsize=((x_max - x_min)/100*scale, (y_max - y_min)/100*scale))

	ax.set(xlim=[x_min-x_pad, x_max+x_pad], ylim=[y_min-y_pad, y_max+y_pad])
	ax.invert_yaxis()
	ax.axis('off')

	return fig, ax