#include "graph_utils.h"
#include "RS.h"

void kNN_search(int query_id, const Point& query, const Tree& kdTree,
	const Distance& tr_dist, int k, std::vector<int>& neighbors,
	std::vector<float>& neighbor_distance, bool isContain) {
	if (!isContain)
		k -= 1;
	K_neighbor_search search(kdTree, boost::tuple<Point, int>(query, 0), k + 1);
	for (K_neighbor_search::iterator it = search.begin(); it != search.end(); it++) {
		if (boost::get<1>(it->first) == query_id && isContain)
			continue;
		neighbor_distance.push_back(tr_dist.inverse_of_transformed_distance(it->second));
		neighbors.push_back(boost::get<1>(it->first));
	}
	return;
}

void kNN_search(int query_id, const Point& query, const Tree& kdTree,
	const Distance& tr_dist, float radius, std::vector<int>& neighbors,
	std::vector<float>& neighbor_distance, bool isContain) {

	Sphere s(query, radius, 0.);
	
	std::vector<boost::tuple<Point, int>> out_search;
	kdTree.search(std::back_inserter(out_search), s);

	for (auto& result : out_search) {
		if (result.tail.head == query_id && isContain)
			continue;
		neighbors.push_back(result.tail.head);
		neighbor_distance.push_back(norm(result.head - query));
	}
	return;
}

float find_components(std::vector<Point>& vertices,
	std::vector<std::vector<Point>>& component_vertices,
	std::vector<Point>& smoothed_v,
	std::vector<std::vector<Point>>& component_smoothed_v,
	std::vector<Vector>& normals,
	std::vector<std::vector<Vector>>& component_normals,
	const Tree& kdTree, const Distance& tr_dist,
	int k, bool isGTNormal, bool isEuclidean) {
	float avg_edge_length = 0;
	Graph components(smoothed_v.size());
	int this_idx = 0;
	std::set<int> dup_remove;
	// Construct graph
	for (auto& vertex : smoothed_v) {
		if (dup_remove.find(this_idx) != dup_remove.end()) {
			this_idx++;
			continue;
		}

		std::vector<int> neighbors;
		std::vector<float> neighbor_distance;
		kNN_search(this_idx, vertex, kdTree, tr_dist, k, neighbors, neighbor_distance);

		// Filter out cross connection
		{
			std::vector<int> temp;
			Vector this_normal = normals[this_idx];
			for (int j = 0; j < neighbors.size(); j++) {
				int idx = neighbors[j];
				Vector neighbor_normal = normals[idx];
				float cos_theta = this_normal * neighbor_normal /
					std::sqrt(this_normal.squared_length()) /
					std::sqrt(neighbor_normal.squared_length());
				float cos_thresh = std::cos(60. / 180. * CGAL_PI);
				if (isEuclidean)
					cos_thresh = 0.;
				if (cos_theta >= cos_thresh) {
					temp.push_back(idx);
				}
			}
			if (temp.size() == 0) {
				neighbors.clear();
				//dup_remove.insert(this_idx);
				//this_idx++;
				//continue;
			}
			else {
				neighbors.clear();
				neighbors = temp;
			}
			//neighbors = temp;
		}

		for (int i = 0; i < neighbors.size();i++) {
			int idx = neighbors[i];
			float length = neighbor_distance[i];

			if (this_idx == idx)
				continue;

			// Remove duplicate vertices
			if (length < 1e-8 && this_idx != idx) {
				dup_remove.insert(idx);
				//std::cout << "Duplicate vertex " << idx << " is removed" << std::endl;
			}

			avg_edge_length += length;
			
			for (int i = 0; i < k; i++) {
				if (boost::edge(this_idx, idx, components).second)
					continue;
				boost::add_edge(this_idx, idx, components);
			}
		}
		this_idx++;
	}
	std::cout << std::to_string(dup_remove.size()) << " duplicate vertices will be removed." << std::endl;
	float thresh_r = avg_edge_length / boost::num_edges(components) * 20;

	// Remove Edge Longer than threshold
	auto es = boost::edges(components);
	std::vector<int> edge_rm_v_id1, edge_rm_v_id2;
	for (auto eit = es.first; eit != es.second; eit++) {
		int vertex1 = boost::source(*eit, components);
		int vertex2 = boost::target(*eit, components);
		float edge_length = norm(vertices[vertex1] -vertices[vertex2]);
		if (dup_remove.find(vertex1) != dup_remove.end()) {
			//boost::remove_edge(*eit, components);
			edge_rm_v_id1.push_back(vertex1);
			edge_rm_v_id2.push_back(vertex2);
			continue;
		}	
		if (dup_remove.find(vertex2) != dup_remove.end()) {
			edge_rm_v_id1.push_back(vertex1);
			edge_rm_v_id2.push_back(vertex2);
			continue;
		}	
		if (edge_length > thresh_r){
			edge_rm_v_id1.push_back(vertex1);
			edge_rm_v_id2.push_back(vertex2);
        }
	}

	for (int i = 0; i < edge_rm_v_id1.size(); i++) {
		boost::remove_edge(edge_rm_v_id1[i], edge_rm_v_id2[i], components);
	}

	// Find Components
	std::vector<int> component_id(boost::num_vertices(components));
	int num = boost::connected_components(components, &component_id[0]);
	std::vector<int> num_board(num);
	for (int i = 0; i < component_id.size(); i++) {
		int id = component_id[i];
		num_board[id]++;
	}
	std::cout << "The input contains " << num << " connected components." << std::endl;

	// Valid Components
	int valid_num = 0;
	std::map<int, int> valid_id_map;
	int idx = 0;
	int threshold = std::min<int>(vertices.size(), 100);
	for (auto id_num : num_board) {
		if (id_num >= threshold) {
			valid_id_map[idx] = valid_num;
			valid_num++;
		}
		idx++;
	}
	std::cout << std::to_string(valid_num) << " of them will be reconstructed." << std::endl;
	component_normals = std::vector<std::vector<Vector>>(valid_num, std::vector<Vector>{});
	component_vertices = std::vector<std::vector<Point>>(valid_num, std::vector<Point>{});
	component_smoothed_v = std::vector<std::vector<Point>>(valid_num, std::vector<Point>{});

	// New vector for components
	for (int i = 0; i < component_id.size();i++) {
		int id = component_id[i];
		if (dup_remove.find(i) != dup_remove.end())
			continue;
		if (valid_id_map.find(id)!= valid_id_map.end()) {
			component_vertices[valid_id_map[id]].push_back(vertices[i]);
			component_smoothed_v[valid_id_map[id]].push_back(smoothed_v[i]);
			if(normals.size()==vertices.size())
				component_normals[valid_id_map[id]].push_back(normals[i]);
		}
	}
	components.clear();
	return thresh_r;
}

void build_dist_angle_matrix(const std::vector<Point>& vertices, 
	const std::vector<Vector>& normals, const Tree& kdTree, const Distance& tr_dist,
	int k, std::vector<std::vector<float>>& dist_matrix,
	std::vector<std::vector<int>>& id_matrix,
	std::vector<std::vector<float>>& angle_matrix,
	bool isGTNormal, bool isEuclidean) {
	int i = 0;
	for (auto vertex : vertices) {
		Vector this_normal = normals[i];

		//K_neighbor_search search(kdTree, query, k + 1, 0, true, tr_dist);
		std::vector<int> neighbors;
		std::vector<float> dists;
		kNN_search(i, vertex, kdTree, tr_dist, k, neighbors, dists);

		// Filter out cross connection
		if (isGTNormal) {
			std::vector<int> temp;
			for (int j = 0; j < neighbors.size(); j++) {
				int idx = neighbors[j];
				Vector neighbor_normal = normals[idx];
				float cos_theta = this_normal * neighbor_normal /
					std::sqrt(this_normal.squared_length()) /
					std::sqrt(neighbor_normal.squared_length());
				float cos_thresh = std::cos(120. / 180. * CGAL_PI);
				if (cos_theta >= cos_thresh) {
					temp.push_back(idx);
				}
			}
			if (temp.size() == 0)
				std::cout << "Bad normal input" << std::endl;
			else {
				neighbors.clear();
				neighbors = temp;
			}
		}

		std::vector<float> dist_vector;
		std::vector<int> id_vector;
		std::vector<float> angle_vector;
		for (int j = 0; j < neighbors.size(); j++) {
			int idx = neighbors[j];
			Vector neighbor_normal = normals[idx];
			Point neighbor_pos = vertices[idx];
			Vector edge = neighbor_pos - vertex;
			float Euclidean_dist = norm(edge);
			if (isEuclidean)
				dist_vector.push_back(Euclidean_dist);
			else {
				float neighbor_normal_length = edge * normalize_vector(neighbor_normal);
				float normal_length = edge * normalize_vector(this_normal);
				float projection_dist = sqrtf((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
				projection_dist += sqrtf((Euclidean_dist * Euclidean_dist) -
					(neighbor_normal_length * neighbor_normal_length));
				projection_dist /= 2.;
				if (std::abs(normalize_vector(this_normal) * normalize_vector(neighbor_normal)) < 0.71)
					projection_dist = Euclidean_dist;
				dist_vector.push_back(projection_dist);
			}

			float angle_weight = cal_angle_based_weight(normalize_vector(this_normal),
				normalize_vector(neighbor_normal));
			if (angle_weight < 0)
				std::cout << "error" << std::endl;
			angle_vector.push_back(angle_weight);
			id_vector.push_back(idx);
		}
		dist_matrix.push_back(dist_vector);
		id_matrix.push_back(id_vector);
		angle_matrix.push_back(angle_vector);
		/*if (i == 253059) {
			Vector edge1 = vertices[60179] - vertex;
			float Eucli_dist1 = norm(edge1);
			std::cout << Eucli_dist1 << std::endl;
			float neighbor_normal_length = edge1 * normalize_vector(normals[60179]);
			std::cout << neighbor_normal_length << std::endl;
			float normal_length = edge1 * normalize_vector(this_normal);
			std::cout << normal_length << std::endl;
			float projection_dist = std::sqrtf((Eucli_dist1 * Eucli_dist1) - (normal_length * normal_length));
			projection_dist += std::sqrtf((Eucli_dist1 * Eucli_dist1) -
				(neighbor_normal_length * neighbor_normal_length));
			projection_dist /= 2.;
			std::cout << projection_dist << std::endl;
			if (normalize_vector(this_normal) * normalize_vector(normals[60179]) < 0.71)
				projection_dist = Eucli_dist1;

			Vector edge2 = vertices[253060] - vertex;
			float Eucli_dist2 = norm(edge2);
			std::cout << Eucli_dist2 << std::endl;
			neighbor_normal_length = edge2 * normalize_vector(normals[253060]);
			std::cout << neighbor_normal_length << std::endl;
			normal_length = edge2 * normalize_vector(this_normal);
			std::cout << normal_length << std::endl;
			projection_dist = std::sqrtf((Eucli_dist2 * Eucli_dist2) - (normal_length * normal_length));
			projection_dist += std::sqrtf((Eucli_dist2 * Eucli_dist2) -
				(neighbor_normal_length * neighbor_normal_length));
			projection_dist /= 2.;
			std::cout << projection_dist << std::endl;
			if (normalize_vector(this_normal) * normalize_vector(normals[253060]) < 0.71)
				projection_dist = Eucli_dist2;
		}*/
		i++;
	}
	return;
}

void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
	const std::vector<Vector>& normals, const Tree& kdTree, const Distance& tr_dist,
	int k, s_Graph& dist_graph, s_weightMap& weightmap,
	bool isGTNormal, bool isEuclidean, std::vector<float>& max_length,
	int exp_genus, std::vector<float>& pre_max_length) {
	dist_graph = s_Graph(vertices.size());
	int this_k = k;
	if (exp_genus != 0) {
		this_k = 18;
	}

	int i = 0;
	for (auto& vertex : vertices) {
		Vector this_normal = normals[i];

		std::vector<int> neighbors;
		std::vector<float> dists;
		kNN_search(i, smoothed_v[i], kdTree, tr_dist, k, neighbors, dists);
		pre_max_length[i] = dists[int(k * 2. / 3.)];

		// Filter out cross connection
		{
			std::vector<int> temp;
			for (int j = 0; j < neighbors.size(); j++) {
				int idx = neighbors[j];
				Vector neighbor_normal = normals[idx];
				float cos_theta = this_normal * neighbor_normal /
					std::sqrt(this_normal.squared_length()) /
					std::sqrt(neighbor_normal.squared_length());
				float cos_thresh = std::cos(60. / 180. * CGAL_PI);
				if (isEuclidean)
					cos_thresh = 0.;
				if (cos_theta >= cos_thresh) {
					temp.push_back(idx);
				}
			}
			if (temp.size() == 0)
				std::cout << "Bad normal input" << std::endl;
			else {
				neighbors.clear();
				neighbors = temp;
			}
		}

		for (int j = 0; j < neighbors.size(); j++) {
			int idx = neighbors[j];
			if (boost::edge(i, idx, dist_graph).second)
				continue;
			if (idx == i) {
				std::cout << "Vertex " << idx << " connect back to its own." << std::endl;
				continue;
			}
			Vector neighbor_normal = normals[idx];
			Point neighbor_pos = vertices[idx];
			Vector edge = neighbor_pos - vertex;
			float Euclidean_dist = norm(edge);
			float weight = Euclidean_dist;
			if (!isEuclidean) {
				weight = cal_proj_dist(edge, this_normal, neighbor_normal);
			}
			if (weight > max_length[i])
				max_length[i] = weight;
			if (weight > max_length[idx])
				max_length[idx] = weight;
			if (weight < 1e-8)
				std::cout << "error" << std::endl;

			boost::graph_traits< s_Graph >::edge_descriptor e;
			bool inserted;
			boost::tie(e, inserted) = boost::add_edge(i, idx, dist_graph);
			weightmap[e] = weight;
		}
		i++;
	}
	return;
}

void unordered_set_intersection(std::unordered_set<Vertex>& first,
	std::unordered_set<Vertex>& second,
	std::vector<Vertex>& intersection) {
	for (auto i = first.begin(); i != first.end(); i++) {
		if (second.find(*i) != second.end()) 
			intersection.push_back(*i);
	}
}

void unordered_set_intersection(std::unordered_set<int>& first,
	std::unordered_set<int>& second,
	std::vector<int>& intersection) {
	for (auto i = first.begin(); i != first.end(); i++) {
		if (second.find(*i) != second.end())
			intersection.push_back(*i);
	}
}

void unordered_set_intersection(std::unordered_set<int>& first,
	std::unordered_set<int>& second,
	std::unordered_set<int>& intersection) {
	for (auto i = first.begin(); i != first.end(); i++) {
		if (second.find(*i) != second.end())
			intersection.insert(*i);
	}
}

void print_path(std::vector<Vertex>& pred, Vertex target, std::vector<Vertex>& path) {
	if (pred[target] == target) {
		path.push_back(target);
		return;
	}
	print_path(pred, pred[target], path);
	path.push_back(target);
}


int find_shortest_path(const m_Graph& mst, Vertex this_v, Vertex neighbor, int threshold, std::vector<Vertex>& path) {
	// Init
	int shortest_path_step = -1;
	std::vector<Vertex> pred(boost::num_vertices(mst.graph), mst.graph.null_vertex());
	std::vector<int> dist(boost::num_vertices(mst.graph), -1);
	my_visitor vis{ neighbor, threshold, dist.data()};
	auto predmap = pred.data();
	auto distmap = dist.data();

	try {
		boost::dijkstra_shortest_paths(mst.graph, this_v,
			boost::distance_map(distmap).predecessor_map(predmap).visitor(vis).
			weight_map(boost::get(&EdgeProperty::count_weight, mst.graph))
		);

		std::cout << "No path found" << std::endl;
	}
	catch (my_visitor::found const&) {
		shortest_path_step = distmap[neighbor];
		print_path(pred, neighbor, path);
	}
	catch (my_visitor::enough const&) {
		shortest_path_step = -1;
	}
	return shortest_path_step;
}

void build_dist_map(const m_Graph& mst, std::vector<int>& dist, Vertex this_v, int threshold) {
	// Init
	my_visitor vis{ 0, threshold, dist.data() };
	auto distmap = dist.data();

	try {
		boost::dijkstra_shortest_paths(mst.graph, this_v,
			boost::distance_map(distmap).visitor(vis).
			weight_map(boost::get(&EdgeProperty::count_weight, mst.graph))
		);

		std::cout << "No path found" << std::endl;
	}
	catch (my_visitor::enough const&) {
		return;
	}
}

Vertex travel_ccw(m_Graph& g, Vertex current_vertex, Vertex last_vertex) {
	/*std::vector<int> current_neighbors = g.graph[current_vertex].ordered_neighbors;
	int idx = -1;
	for (int i = 0; i < current_neighbors.size(); i++) {
		if (current_neighbors[i] == int(last_vertex)) {
			idx = i;
			break;
		}
	}
	if (idx == -1)
		throw std::runtime_error("No neighbor found (caused by wrong ordered_neighbors)");
	int this_v_id = g.graph[current_vertex].
		ordered_neighbors[python_mod((idx + 1), current_neighbors.size())];
	return this_v_id;
	 */
    return successor(g,current_vertex,last_vertex).v;
}

Vertex travel_cw(m_Graph& g, Vertex current_vertex, Vertex last_vertex) {
	/*std::vector<int> current_neighbors = g.graph[current_vertex].ordered_neighbors;
	int idx = -1;
	for (int i = 0; i < current_neighbors.size(); i++) {
		if (current_neighbors[i] == int(last_vertex)) {
			idx = i;
			break;
		}
	}
	if (idx == -1)
		throw std::runtime_error("No neighbor found (caused by wrong ordered_neighbors)");
	int this_v_id = g.graph[current_vertex].
		ordered_neighbors[python_mod((idx - 1), current_neighbors.size())];
	return this_v_id;
	 */
    return predecessor(g,current_vertex,last_vertex).v;
}

float cal_proj_dist(const Vector& edge, Vector& this_normal, Vector& neighbor_normal) {
	float Euclidean_dist = norm(edge);
	/*Vector sum_normal;
	if (this_normal * neighbor_normal > 0)
		sum_normal = normalize_vector(normalize_vector(this_normal) +
			normalize_vector(neighbor_normal));
	else
		sum_normal = normalize_vector(normalize_vector(this_normal) -
			normalize_vector(neighbor_normal));
	float normal_length = edge * sum_normal;
	float projection_dist = sqrtf((Euclidean_dist * Euclidean_dist)
		- (normal_length * normal_length));*/
	float neighbor_normal_length = edge * normalize_vector(neighbor_normal);
	float normal_length = edge * normalize_vector(this_normal);
	float projection_dist = sqrtf((Euclidean_dist * Euclidean_dist) - (normal_length * normal_length));
	projection_dist += sqrtf((Euclidean_dist * Euclidean_dist) -
		(neighbor_normal_length * neighbor_normal_length));
	projection_dist /= 2.;
	if (std::abs(normalize_vector(this_normal) * normalize_vector(neighbor_normal)) < std::cos(15./180.* CGAL_PI))
		projection_dist = Euclidean_dist;
	return projection_dist;
}