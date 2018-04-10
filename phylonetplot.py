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

	@property
	def is_trivial(self):
		return self.start == self.end

	def __str__(self):
		return 'Split {0}-{1}: {2}'.format(self.start, self.end, self.weight)


class SplitSystem:

	def __init__(self, samples):
		self._splits = []
		self._samples = samples

	def add_split(self, start, end, weight):
		print(start, end, weight)
		if end-start > self.num_samples - 1:
			raise InvalidSplit
		self._splits.append(Split(start, end, weight))

	def sort(self):
		self._splits = sorted(self._splits, key=lambda x: x.start-x.end)

	@property
	def num_splits(self):
		return len(self._splits)

	@property
	def num_samples(self):
		return len(self._samples)

	@property
	def trivial_splits(self):
		return [(split, self._samples[split.start]) for split in self._splits if split.is_trivial]

	@property
	def non_trivial_splits(self):
		return [split for split in self._splits if not split.is_trivial]

	def get_split_sample_range(self, split):
		return self._samples[split.start:split.end+1]


class Drawnet:

	def __init__(self, samples):
		self._samples = samples
		self._splits = SplitSystem(self._samples)

		plt.figure(figsize=(10, 10))
		self._graph = nx.Graph()

		self._pos = {}
		self._v_labels = {}
		self._e_labels = {}
		self._mass = {}




	def add_split_matrix(self, splits, weights):

		melt = []

		for s in splits:
			s = list(s)
			if s.count(1) == 1:
				start = s.index(1)
				end = start
			elif s.count(0) == 1:
				start = s.index(0)
				end = start
			else:
				start = s.index(1) if s[0] == 0 else s.index(0)
				end = len(s) - (s[::-1].index(1) if s[-1] == 0 else s[::-1].index(0)) -1

				if end < start:
					start = 0
			print(start,end)
			melt.append((start,end))

		self._process_splits_tuples(melt, weights)

	def _process_splits_tuples(self, splits, weights):
		for s, w in zip(splits, weights):
			self._splits.add_split(s[0], s[1], w)
		self._splits.sort()

	@property
	def num_samples(self):
		return len(self._samples)

	def _set_populations(self, populations):
		self._pops = list(populations)

		NUM_COLORS = len(np.unique(self._pops))

		cm = plt.get_cmap('gist_rainbow')

		cols = [cm(1. * i / NUM_COLORS) for i in range(NUM_COLORS)]

		counter = 0
		categories = {}
		self._colormap = []
		for x in self._pops:
			if x in categories:
				self._colormap.append(categories[x])
			else:
				categories[x] = cols[counter]
				self._colormap.append(categories[x])
				counter += 1

	@staticmethod
	def _dist(a, b):
		return np.linalg.norm([x - y for x, y in zip(a, b)])

	def _add_center_node(self):
		self._graph.add_node(0)
		self._pos[0] = (0, 0)
		self._v_labels[0] = ''

	def _add_leaves(self):

		angle = 0  # first node at 3 o'clock (0 radians)

		done = []

		for split, sample in self._splits.trivial_splits:
			done.append(sample)

			self._graph.add_node(sample)
			self._pos[sample] = (math.cos(angle) * split.weight, math.sin(angle) * split.weight)
			self._mass[sample] = (math.cos(angle), math.sin(angle))
			self._v_labels[sample] = sample

			angle = angle + 2.0 * math.pi / len(self._samples)  # orient leafs in circular order clock-wise
			self._graph.add_edge(0, sample)

		for sample in self._samples:
			if sample not in done:
				self._graph.add_node(sample)
				self._pos[sample] = (math.cos(angle) * split.weight, math.sin(angle) * split.weight)
				self._mass[sample] = (math.cos(angle), math.sin(angle))
				self._v_labels[sample] = sample

				angle = angle + 2.0 * math.pi / len(self._samples)  # orient leafs in circular order clock-wise
				self._graph.add_edge(0, sample)


	def _add_interior_nodes(self):

		for ii, split in enumerate(self._splits.non_trivial_splits):

			splitpart = self._splits.get_split_sample_range(split)

			pathnodes = []
			for i in range(split.start+1, split.end+1):
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
			print(set(self._samples) - set(splitpart))
			vector = [x/np.linalg.norm(vector) for x in vector]

			# update positions
			for n in move_nodes:
				self._pos[n] = [self._pos[n][i]-vector[i]*split.weight for i in range(2)]

	def draw(self):
		self._add_center_node()
		self._add_leaves()
		self._add_interior_nodes()

		# draw graph
		nx.draw(self._graph,
				pos=self._pos,
				labels=self._v_labels,
				node_size=0,
				width=1,
				with_labels=False)

		for s, c in zip(self._samples, self._colormap):
			nx.draw_networkx_labels(self._graph, {s: self._pos[s]}, {s: self._v_labels[s]}, font_size=10, font_color=c)

		fn = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(3))
		self.netsp = 1
		plt.show()
		#plt.savefig(str(self.netsp)+"nn"+fn+".pdf")

'''
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
p = ['C0','C1', 'C2','C1','C3']

dn = Drawnet(o, w, s)
dn._set_populations(p)
dn.draw()
'''