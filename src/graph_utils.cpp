#include "graph_utils.h"
#include "RS.h"

/**
 * @brief k nearest neighbor search
 *
 * @param query_id: the index of the point to be queried
 * @param query: the coordinate of the point to be queried
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param k: number of nearest neighbors to be queried
 * @param neighbors: [OUT] indices of k nearest neighbors
 * @param neighbor_distance: [OUT] corresponding distance to the query point
 * @param isContain: does the query point itself count as a neighbor
 *
 * @return None
 */
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

/**
 * @brief neighbor search within a specific radius
 *
 * @param query_id: the index of the point to be queried
 * @param query: the coordinate of the point to be queried
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param radius: the radius of the search ball
 * @param neighbors: [OUT] indices of k nearest neighbors
 * @param neighbor_distance: [OUT] corresponding distance to the query point
 * @param isContain: does the query point itself count as a neighbor
 *
 * @return None
 */
void NN_search(int query_id, const Point& query, const Tree& kdTree,
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

/**
 * @brief Find the number of connected components and separate them
 *
 * @param vertices: vertices of the point cloud
 * @param component_vertices: [OUT] point cloud vertices in different connected components
 * @param smoothed_v: smoothed vertices of the point cloud
 * @param component_smoothed_v: [OUT] smoothed point cloud vertices in different connected components
 * @param normals: normal of the point cloud vertices
 * @param component_normals: [OUT] normal of the point cloud vertices in different components
 * @param kdTree: kd-tree for neighbor query
 * @param tr_dist: distance container
 * @param k: number of nearest neighbors to be queried
 * @param isEuclidean: if Euclidean distance is used as edge weight
 * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
 * @param outlier_thresh: threshold to remove outlier
 * 
 *
 * @return None
 */
float find_components(std::vector<Point>& vertices,
	std::vector<std::vector<Point>>& component_vertices,
	std::vector<Point>& smoothed_v,
	std::vector<std::vector<Point>>& component_smoothed_v,
	std::vector<Vector>& normals,
	std::vector<std::vector<Vector>>& component_normals,
	const Tree& kdTree, const Distance& tr_dist,
	int k, bool isEuclidean, float cross_conn_thresh, float outlier_thresh) {
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
				float cos_thresh = std::cos(cross_conn_thresh / 180. * CGAL_PI);
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
	float thresh_r = avg_edge_length / boost::num_edges(components) * outlier_thresh;

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

/**
 * @brief initialize the graph and related information
 *
 * @param vertices: vertices of the componnet
 * @param smoothed_v: smoothed vertices of the component
 * @param normals: normal of the component vertices
 * @param kdTree: kd-tree for neighbor query
 * @param tr_dist: distance container
 * @param k: number of nearest neighbors to be queried
 * @param dist_graph: [OUT] a light-weight graph with essential connection for building MST
 * @param weightmap: [OUT] edge weight of dist_graph
 * @param isEuclidean: if it is using Euclidean distance
 * @param max_length: [OUT] the distance of the longest connection each vertex involved
 * @param exp_genus: user-specified expected genus number
 * @param pre_max_length: [OUT] the maximum length of connection before connecting handles (conservative connection)
 * @param cross_conn_thresh: angle threshold to avoid connecting vertices on different surface
 *
 *
 * @return None
 */
void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
	const std::vector<Vector>& normals, const Tree& kdTree, const Distance& tr_dist,
	int k, s_Graph& dist_graph, s_weightMap& weightmap, bool isEuclidean, std::vector<float>& max_length,
	int exp_genus, std::vector<float>& pre_max_length, float cross_conn_thresh) {
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

/**
 * @brief Retreive the shortest path in the graph
 *
 * @param pred: information map
 * @param target: target vertex
 * @param path: [OUT] vertex indices of the shortest path
 *
 * @return None
 */
void print_path(std::vector<Vertex>& pred, Vertex target, std::vector<Vertex>& path) {
	if (pred[target] == target) {
		path.push_back(target);
		return;
	}
	print_path(pred, pred[target], path);
	path.push_back(target);
}

/**
 * @brief Find the shortest path from one vertex to another in the graph
 *
 * @param mst: the graph
 * @param this_v: the source vertex
 * @param neighbor: the target vertex
 * @param threshold: the step threshold, if longer than this threshold, the algorithm early stop
 * @param path: stores the indices of vertex in the shortest path (for visualization)
 *
 * @return the number of steps of the shortest path
 */
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

/**
 * @brief weighted smoothing method using defined neighborhood with tangential distance weighted
 *
 * @param vertices: vertices of the point cloud
 * @param smoothed_v: [OUT] vertices after smoothing
 * @param normals: normal of the point cloud
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param diagonal_length: the diagonal length of the point cloud
 *
 * @return None
 */
void weighted_smooth(const std::vector<Point>& vertices,
	std::vector<Point>& smoothed_v, const std::vector<Vector>& normals,
	const Tree& kdTree, const Distance& tr_dist, float diagonal_length) {
	int idx = 0;
	for (auto& vertex : vertices) {
		std::vector<int> neighbors;
		std::vector<float> neighbor_dist;
		Vector normal = normals[idx];

		int neighbor_num = 192;
		kNN_search(idx, vertex, kdTree, tr_dist, neighbor_num,
			neighbors, neighbor_dist);

		float weight_sum = 0.;
		float amp_sum = 0.;
		float max_dist = 0.;
		int added = 0;
		std::vector<float> vertical_length;
		std::vector<float> weights;
		for (auto& neighbor : neighbors) {
			Point neighbor_pos = vertices[neighbor];
			Vector n2this = neighbor_pos - vertex;
			if (normals[neighbor] * normal < std::cos(30. / 180. * CGAL_PI)) {
				continue;
			}
			float vertical = n2this * normal;
			float n_dist = norm(neighbor_pos - vertex);

			float tangential_square = n_dist * n_dist -
				vertical * vertical;
			float tangential_dist = 0.;
			if (tangential_square > 0.)
				tangential_dist = std::sqrt(tangential_square);
			if (!isfinite(tangential_dist)) {
				std::cout << n_dist << " " << vertical << std::endl;
				std::cout << "error" << std::endl;
			}
			float weight = -tangential_dist;
			if (tangential_dist > max_dist)
				max_dist = tangential_dist;

			weights.push_back(weight);
			vertical_length.push_back(vertical);
			added++;
		}
		for (int i = 0; i < vertical_length.size(); i++) {
			amp_sum += vertical_length[i] * (weights[i] + max_dist);
			weight_sum += weights[i] + max_dist;
		}

		if (weight_sum == 0.)
			weight_sum = 1.;
		amp_sum /= weight_sum;
		if (!isfinite(amp_sum))
			std::cout << "error" << std::endl;
		Vector move = amp_sum * normal;
		smoothed_v.push_back(vertex + move);
		idx++;
	}
}