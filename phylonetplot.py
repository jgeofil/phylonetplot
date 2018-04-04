import networkx as nx
import matplotlib.pyplot as plt
import random, math, string
import numpy as np


class InvalidSplit(Exception):
	pass


class Split:

	def __init__(self, start, end, weight):

		if weight < 0:
			raise InvalidSplit('Split weight cannot be negative.')

		if start > end:
			raise InvalidSplit('Split is invalid (start > end).')

		self.start = start
		self.end = end
		self.weight = weight

	def is_trivial(self):
		return self.start == self.end


class SplitSystem:

	def __init__(self, samples):
		self._splits = []
		self._samples = samples

	def add_split(self, start, end, weight):
		if end-start > self.num_samples - 2:
			raise InvalidSplit
		self._splits.append(Split(start, end, weight))

	@property
	def num_splits(self):
		return len(self._splits)

	@property
	def num_samples(self):
		return len(self._samples)

	@property
	def trivial_splits(self):
		return [(split, self._samples[split.start]) for split in self._splits if split.is_trival]


class Drawnet:

	def __init__(self, samples, weights, splits):
		self._samples = samples
		self._splits = SplitSystem(self._samples)
		self._process_splits(splits, weights)

		self._graph = nx.Graph()

		self._pos = {}
		self._v_labels = {}
		self._e_labels = {}
		self._mass = {}

	def _process_splits(self, splits, weights):
		for s, w in zip(splits, weights):
			self._splits.add_split(s[0], s[1], w)

	@property
	def num_samples(self):
		return len(self._samples)

	def _set_populations(self, populations):
		self._pops = list(populations)

	@staticmethod
	def _dist(a, b):
		return np.linalg.norm([x - y for x, y in zip(a, b)])

	def _add_center_node(self):
		self._graph.add_node(0)
		self._pos[0] = (0, 0)
		#self._v_labels[0] = ''

	def _add_leaves(self):

		angle = 0  # first node at 3 o'clock (0 radians)

		for split, sample in self._splits.trivial_splits:

			self._graph.add_node(sample)
			self._pos[sample] = (math.cos(angle) * split.weight, math.sin(angle) * split.weight)
			self._mass[sample] = (math.cos(angle), math.sin(angle))
			self._v_labels[sample] = sample

			angle = angle + 2.0 * math.pi / len(self._samples)  # orient leafs in circular order clock-wise
			self._graph.add_edge(0, s)

	def _add_interior_nodes(self):

		for ii, sp in enumerate(self._splits):
			if sp[1] - sp[0] > 0 and sp[1] - sp[0] < len(self._samples)-2 and self._weights[ii] > 0:

				splitpart = [self._samples[i] for i in range(sp[0], sp[1]+1)]

				pathnodes = []
				for i in range(sp[0]+1, sp[1]+1):
					path = nx.shortest_path(self._graph, source=self._samples[i-1], target=self._samples[i])[1:-1]
					pathnodes = pathnodes+[pn for pn in path if pn not in pathnodes]

				# duplicate nodes
				prev = None
				for n in pathnodes:
					newnode = self._graph.number_of_nodes()
					self._pos[newnode] = self._pos[n]
					self._graph.add_node(newnode)

					# update links (some edges will remain at the old node, some
					# are changed to the newly created one

					for nbr in list(self._graph.neighbors(n)):
						if nbr not in splitpart and nbr not in pathnodes:
							self._graph.remove_edge(n, nbr)
							self._graph.add_edge(newnode, nbr)

					# add edge between old and new node
					self._graph.add_edge(n, newnode)

					# add edge between new nodes along the duplicated path
					if prev:
						self._graph.add_edge(newnode, prev)
					prev = newnode

				# move nodes included in the path, as well as nodes included
				# in the split, away from the rest of the nodes
				move_nodes = pathnodes+splitpart

				# compute vector giving the angle to which the nodes are moved
				# based on the angle to which the nodes in each split part were
				# pointing in the star graph in the beginning
				mass_cent_l_x = np.mean([self._mass[n][0] for n in splitpart])
				mass_cent_l_y = np.mean([self._mass[n][1] for n in splitpart])
				mass_cent_r_x = np.mean([self._mass[n][0] for n in set(self._samples) - set(splitpart)])
				mass_cent_r_y = np.mean([self._mass[n][1] for n in set(self._samples) - set(splitpart)])
				vector = [mass_cent_r_x-mass_cent_l_x, mass_cent_r_y-mass_cent_l_y]
				vector = [x/np.linalg.norm(vector) for x in vector]

				# update positions
				for n in move_nodes:
					self._pos[n] = [self._pos[n][i]-vector[i]*self._weights[ii] for i in range(2)]

	def draw(self):

		self._add_leaves()
		self._add_interior_nodes()

		# draw graph
		nx.draw(self._graph,
				pos=self._pos,
				labels=self._v_labels,
				node_size=0,
				width=1,
				with_labels=False)

		nx.draw_networkx_labels(self._graph, self._pos, self._v_labels, font_size=16, font_color='r')

		fn = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
		self.netsp = 1
		plt.show()
		#plt.savefig(str(self.netsp)+"nn"+fn+".pdf")


o = ['a','b','c','d','e']
w = [1,1,3,2,1,2,3,4]
s = [
	(0,0),
	(1,1),
	(2,2),
	(3,3),
	(4,4),
(1,3),
	(0,1),

	(3,4)
]
p = ['r','r','b','b','g']

dn = Drawnet(o, w, s)
dn._set_populations(p)
dn.draw()