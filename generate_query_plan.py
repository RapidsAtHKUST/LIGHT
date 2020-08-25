import copy

p1_vertex_num = 4
p1_edge_list = [(0, 1), (0, 3), (1, 2), (2, 3)]
p1_rules = [(0, 1), (0, 2), (0, 3), (1, 3)]

p2_vertex_num = 4
p2_edge_list = [(0, 1), (0, 2), (0, 3), (1, 2), (2, 3)]
p2_rules = [(0, 2), (1, 3)]

p3_vertex_num = 4
p3_edge_list = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
p3_rules = [(0, 1, 2, 3)]

p4_vertex_num = 5
p4_edge_list = [(0, 1), (0, 2), (0, 4), (1, 3), (1, 4), (2, 3)]
p4_rules = [(0, 1)]

p5_vertex_num = 6
p5_edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (0, 5), (1, 2), (2, 3), (3, 4), (4, 5)]
p5_rules = [(2, 4)]

p6_vertex_num = 5
p6_edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3)]
p6_rules = [(0, 1), (2, 3)]

p7_vertex_num = 5
p7_edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 4)]
p7_rules = [(0, 1, 2, 3, 4)]


p8_vertex_num = 5
p8_edge_list = [(0, 1), (0, 2), (0, 3), (0, 4), (1, 2), (2, 3), (3, 4)]
p8_rules = [(2, 3)]

g_min_cost = -1.0
g_optimal_order = []
g_optimal_operation_order = []
g_bfs_level = []
g_bn = []
g_si_nbr = []
g_si_c = []
g_enable_local_cache = []


def create_pattern(vertex_num, edge_list, rules):
	# convert the edge list to adjacent list.
	adj = [[] for i in range(vertex_num)]
	for edge in edge_list:
		adj[edge[0]].append(edge[1])
		adj[edge[1]].append(edge[0])

	for nbrs in adj:
		nbrs.sort()

	# create rules for each vertex.
	dag = [[] for i in range(vertex_num)]
	for rule in rules:
		for i in range(len(rule)):
			for j in range(i + 1, len(rule)):
				dag[rule[j]].append(rule[i])

	for nbrs in dag:
		nbrs.sort()

	return adj, dag


def is_set_cover(status):
	for key, value in status.iteritems():
		if value == 0:
			return False
	return True


def bfs_level(root, adj, level):
	n = len(adj)
	visited = [False] * n
	bfs_queue = []

	bfs_queue.append(root)
	visited[root] = True
	level[root] = 0

	while len(bfs_queue) != 0:
		cur_vertex = bfs_queue.pop(0)
		cur_level = level[cur_vertex]
		for nbr in adj[cur_vertex]:
			if not visited[nbr]:
				bfs_queue.append(nbr)
				visited[nbr] = True
				level[nbr] = cur_level + 1


def estimate_cost(operation_order, bn, si_nbr, si_c):
	# Expand factor and reduction factor. We use fixed values on the two factors instead of the complex estimation
	# SEED, because we find that it is hard to give an accurate estimation on the two values just based on the distribution
	# of degrees. We make the expand factor to be large so that the cost on computation generally dominate the total cost.
	a = 1000.0
	b = 0.1

	mat_cost = 0
	comp_cost = 0
	estimated_size = 1

	# Determine whether the pattern is a star that we can use local cache (set intersection cache) to reduce the number
	# of set intersections.
	enable_local_cache = False
	for item in operation_order[2:]:
		vertex = item[0]
		operation = item[1]
		if operation == 'COMP':
			if len(bn[vertex]) == 2 and operation_order[0][0] == bn[vertex][0] and bn[vertex][1] in adj[bn[vertex][0]]:
				enable_local_cache = True
			elif len(bn[vertex]) == 1:
				continue
			else:
				enable_local_cache = False
				break

	if enable_local_cache:
		comp_cost += estimated_size * a * a
	for item in operation_order[2:]:
		vertex = item[0]
		operation = item[1]
		if operation == 'MAT':
			# We set the minimum value (10) to be greater than 1, because we find that for the small pattern,
			# the number of matches increases with the number of vertices in the pattern.
			estimated_size = estimated_size * (max(a * (b ** (len(bn[vertex]) - 1)), 10.0))
			mat_cost += estimated_size
		elif operation == 'COMP':
			if not enable_local_cache:
				comp_cost += estimated_size * a * (max(len(si_nbr[vertex]) + len(si_c[vertex]) - 1, 0))
	total_cost = mat_cost + comp_cost

	return total_cost, enable_local_cache


def recursive_msc(collection, msc, status, current_solution, depth):
	if len(current_solution) + 1 >= len(msc) or depth >= len(collection):
		return

	current_value = collection[depth]
	# Add current value.
	current_solution.append(current_value)
	for vertex in current_value[1]:
		status[vertex] += 1
	if is_set_cover(status):
		if len(current_solution) < len(msc):
			msc[:] = []
			msc.extend(current_solution)
			current_solution.pop()
			for vertex in current_value[1]:
				status[vertex] -= 1
			return
	recursive_msc(collection, msc, status, current_solution, depth + 1)
	current_solution.pop()
	for vertex in current_value[1]:
		status[vertex] -= 1

	# Do not add current value.
	recursive_msc(collection, msc, status, current_solution, depth + 1)


def find_msc(universe, collection):
	current_solution = []
	msc = [(vertex, {vertex}) for vertex in universe]
	status = {vertex: 0 for vertex in universe}
	recursive_msc(collection, msc, status, current_solution, 0)
	return msc


def generate_operation_order(vertex_number, order, adj):
	materialized = [False] * vertex_number
	bn = [[] for i in range(vertex_number)]

	operation_order = [(order[0], 'COMP')]

	for i in range(1, vertex_number):
		vertex = order[i]
		for j in range(i):
			visited_vertex = order[j]
			if visited_vertex in adj[vertex]:
				bn[vertex].append(visited_vertex)
				if not materialized[visited_vertex]:
					materialized[visited_vertex] = True
					operation_order.append((visited_vertex, 'MAT'))

		operation_order.append((vertex, 'COMP'))

	si_nbr = [[] for i in range(vertex_number)]
	si_c = [[] for i in range(vertex_number)]
	for i in range(vertex_number):
		vertex = order[i]
		if not materialized[vertex]:
			materialized[vertex] = True
			operation_order.append((vertex, 'MAT'))

	for i in range(1, vertex_number):
		vertex = order[i]
		universe = set(bn[vertex])
		collection = []
		for u in bn[vertex]:
			collection.append((u, {u}))

		for j in range(0, i):
			u = order[j]
			if len(bn[u]) < 2:
				continue
			u_bn = set(bn[u])
			if u_bn <= universe:
				valid = True
				for value in collection:
					if u_bn == value[1]:
						valid = False
						break
				if valid:
					collection.insert(0, (u, u_bn))

		msc = find_msc(universe, collection)
		for value in msc:
			if len(value[1]) == 1:
				si_nbr[vertex].append(value[0])
			else:
				si_c[vertex].append(value[0])
	return operation_order, bn, si_nbr, si_c


def recursive_permutation(order, visited, vertex_number, adj, dag, cur_bfs_level):
	if len(order) == vertex_number:
		operation_order, bn, si_nbr, si_c = generate_operation_order(vertex_number, order, adj)
		cur_cost, enable_local_cache = estimate_cost(operation_order, bn, si_nbr, si_c)

		global g_min_cost
		global g_optimal_order
		global g_optimal_operation_order
		global g_bn
		global g_si_nbr
		global g_si_c
		global g_bfs_level
		global g_enable_local_cache
		if g_min_cost < 0 or cur_cost < g_min_cost:
			g_min_cost = cur_cost
			g_optimal_order[:] = []
			g_optimal_operation_order[:] = []
			g_bn[:] = []
			g_si_nbr[:] = []
			g_si_c[:] = []
			g_bfs_level[:] = []
			g_enable_local_cache[:] = []
			g_optimal_order.append(copy.deepcopy(order))
			g_optimal_operation_order.append(operation_order)
			g_bn.append(bn)
			g_si_nbr.append(si_nbr)
			g_si_c.append(si_c)
			g_bfs_level.append(cur_bfs_level)
			g_enable_local_cache.append(g_enable_local_cache)
		elif cur_cost == g_min_cost:
			# Break ties by prioritizing the order that puts the vertices constrained by partial orders before other vertices.
			is_tie = 0

			for i in range(vertex_number):
				min_vertex = g_optimal_order[0][i]
				cur_vertex = order[i]

				if len(dag[cur_vertex]) > len(dag[min_vertex]):
					is_tie = -1
					break
				elif len(dag[cur_vertex]) < len(dag[min_vertex]):
					is_tie = 1
					break

			# Break ties by prioritizing the order that brings more backward neighbors at an early stage.
			if is_tie == 0:
				for i in range(vertex_number):
					min_vertex = g_optimal_order[0][i]
					cur_vertex = order[i]

					if len(bn[cur_vertex]) > len(g_bn[0][min_vertex]):
						is_tie = -1
						break
					elif len(bn[cur_vertex]) < len(g_bn[0][min_vertex]):
						is_tie = 1
						break

			# Break ties by prioritizing the order that starts with the vertex with a great degree value.
			if is_tie == 0:
				for i in range(vertex_number):
					min_vertex = g_optimal_order[0][i]
					cur_vertex = order[i]

					if len(adj[cur_vertex]) > len(adj[min_vertex]):
						is_tie = -1
						break
					elif len(adj[cur_vertex]) < len(adj[min_vertex]):
						is_tie = 1
						break

			# Break ties by prioritizing the order with the lower bfs order.
			if is_tie == 0:
				for i in range(vertex_number):
					min_vertex = g_optimal_order[0][i]
					cur_vertex = order[i]

					if cur_bfs_level[cur_vertex] < g_bfs_level[0][min_vertex]:
						is_tie = -1
						break
					elif cur_bfs_level[cur_vertex] > g_bfs_level[0][min_vertex]:
						is_tie = 1
						break

			if is_tie < 0:
				g_optimal_order[:] = []
				g_optimal_operation_order[:] = []
				g_bn[:] = []
				g_si_nbr[:] = []
				g_si_c[:] = []
				g_bfs_level[:] = []
				g_enable_local_cache[:] = []
			if is_tie <= 0:
				g_optimal_order.append(copy.deepcopy(order))
				g_optimal_operation_order.append(operation_order)
				g_bn.append(bn)
				g_si_nbr.append(si_nbr)
				g_si_c.append(si_c)
				g_bfs_level.append(cur_bfs_level)
				g_enable_local_cache.append(enable_local_cache)
		return
	for vertex in range(vertex_number):
		# check if the nbr is valid.
		if not visited[vertex]:
			valid = True
			for constraint_vertex in dag[vertex]:
				if not visited[constraint_vertex]:
					valid = False
					break
			if valid:
				valid = False
				for nbr in adj[vertex]:
					if visited[nbr]:
						valid = True
						break
			if valid:
				visited[vertex] = True
				order.append(vertex)
				recursive_permutation(order, visited, vertex_number, adj, dag, cur_bfs_level)
				order.pop()
				visited[vertex] = False


def generate_permutation(adj, dag):
	vertex_number = len(adj)
	order = []
	visited = [False] * vertex_number
	cur_bfs_level = [0] * vertex_number
	for vertex in range(vertex_number):
		valid = True
		for constraint_vertex in dag[vertex]:
			if not visited[constraint_vertex]:
				valid = False
				break
		if valid:
			bfs_level(vertex, adj, cur_bfs_level)
			order.append(vertex)
			visited[vertex] = True
			recursive_permutation(order, visited, vertex_number, adj, dag, cur_bfs_level)
			order.pop()
			visited[vertex] = False


def print_query_plan(enumeration_order, operation_order, backward_neighbors, rules, si_s1, si_s2, enable_local_cache):
	print 'Enumeration Order: {0}'.format(enumeration_order)
	print 'Operation Order: {0}'.format(operation_order)
	print 'Enable Local Cache: {0}'.format(enable_local_cache)
	for i in range(len(enumeration_order)):
		vertex = enumeration_order[i]
		print 'Vertex {0}:'.format(vertex)
		print '{0} Backward Neighbors: {1}'.format(len(backward_neighbors[vertex]), backward_neighbors[vertex])
		print '{0} Set Intersections: Neighbor Sets of {1}, Candidate Sets of {2}'.format(
			max(len(si_s1[vertex]) + len(si_s2[vertex]) - 1, 0), si_s1[vertex], si_s2[vertex])
		print '{0} Symmetry Breaking Rules: {1}'.format(len(rules[vertex]), rules[vertex])
		print '-' * 5


adj, dag = create_pattern(p5_vertex_num, p5_edge_list, p5_rules)
generate_permutation(adj, dag)

print '{0} Optimal Solutions'.format(len(g_optimal_order))
for i in range(len(g_optimal_order)):
	print '-' * 20
	print 'Solution {0}'.format(i)
	print_query_plan(g_optimal_order[i], g_optimal_operation_order[i],
					 g_bn[i], dag, g_si_nbr[i], g_si_c[i], g_enable_local_cache[i])
