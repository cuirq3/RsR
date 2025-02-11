#include "mst.h"

/**
 * @brief Build minimum spanning tree (MST)
 *
 * @param mst: [OUT] constructed MST
 * @param p: connection information of the mst
 * @param isEuclidean: if to use Euclidean distance
 * @param vertices: coordinates of the point cloud
 * @param normals: normal of the point cloud
 *
 * @return None
 */
void build_mst(m_Graph& mst, std::vector<Vertex>& p,
	bool isEuclidean, const std::vector<Point>& vertices,
	std::vector<Vector>& normals ) {
	// Init distance mst
	std::pair<vertex_iter, vertex_iter> vp = boost::vertices(mst.graph);
	mst.vi_start = vp.first;
	mst.vi_end = vp.second;

    for (int i = 0; i < p.size(); i++) {
        mst.graph[i].normal = normals[i];
        mst.graph[i].coords = vertices[i];
        mst.graph[i].id = i;
    }

	Graph mst_temp(boost::num_vertices(mst.graph));
	for (int i = 0; i < p.size(); i++) {
		if (p[i] == i) {
			continue;
		}
		boost::add_edge(i, p[i], mst_temp);
	}
	// Fix strong ambiguous points
	std::vector<std::pair<Vertex, Point>> org_points;
	if (!isEuclidean) {
		auto es = boost::edges(mst_temp);
		for (auto eit = es.first; eit != es.second; eit++) {
			Vertex vertex1 = boost::source(*eit, mst_temp);
			Vertex vertex2 = boost::target(*eit, mst_temp);
			Vector normal1 = normalize_vector(mst.graph[vertex1].normal);
			Vector normal2 = normalize_vector(mst.graph[vertex2].normal);
			Point pos1 = mst.graph[vertex1].coords;
			Point pos2 = mst.graph[vertex2].coords;
			if (boost::degree(vertex1, mst_temp) >= 2 && boost::degree(vertex2, mst_temp) >= 2)
				continue;
			Vector edge = pos2 - pos1;
			float cos_angle = std::abs(edge * normalize_vector(normal1 + normal2) / norm(edge));
			if (cos_angle > std::cos(10. / 180. * CGAL_PI)) {
				Vertex leaf, parent;
				if (boost::degree(vertex1, mst_temp) == 1) {
					mst.graph[vertex1].normal = mst.graph[vertex2].normal;
					parent = vertex2;
					leaf = vertex1;
				}
				else {
					mst.graph[vertex2].normal = mst.graph[vertex1].normal;
					parent = vertex1;
					leaf = vertex2;
				}

				auto neighbors = boost::adjacent_vertices(parent, mst_temp);
				for (auto neighbor : make_iterator_range(neighbors)) {
					if (mst.graph[neighbor].normal_rep == -1) {
						mst.graph[neighbor].normal_rep = parent;
					}
					else {
						// Collision!
						mst.graph[neighbor].normal_rep = -2;
					}
				}
			}
		}
		for (int i = 0; i < boost::num_vertices(mst.graph); i++) {
			if (mst.graph[i].normal_rep >= 0)
				mst.graph[i].normal = mst.graph[mst.graph[i].normal_rep].normal;
		}
	}

	for (int i = 0; i < p.size(); i++) {
		if (p[i] == i) {
			if (p[i] == i)
				std::cout << p[i] << std::endl;
			std::cout << "vertex " + std::to_string(i) + " may be root vertex or seperated" << std::endl;
			continue;
		}
		else {
			Vector edge = mst.graph[i].coords - mst.graph[p[i]].coords;
			float Euclidean_dist = norm(edge);
			float projection_dist = cal_proj_dist(edge, mst.graph[i].normal, mst.graph[p[i]].normal);
			if (std::isnan(projection_dist) || std::isnan(Euclidean_dist))
				std::cout << "debug" << std::endl;

			if (isEuclidean)
				mst.add_edge(i, p[i], Euclidean_dist, true);
			else
				mst.add_edge(i, p[i], projection_dist, true);
		}
	}
	return;
}

/**
 * @brief Check if two segments are intersecting on both planes (defined by their normal) they belong
 *
 * @param mst: graph and vertex information
 * @param v1: 1st vertex of segment 1
 * @param v2: 2nd vertex of segment 1
 * @param v3: 1st vertex of segment 2
 * @param v4: 2nd vertex of segment 2
 *
 * @return if they are intersecting with each other
 */
bool isIntersecting(m_Graph& mst, Vertex v1, Vertex v2, Vertex v3, Vertex v4) {
	Point p1 = mst.graph[v1].coords;
	Point p2 = mst.graph[v2].coords;
	Vector n1 = mst.graph[v1].normal;
	Vector n2 = mst.graph[v2].normal;
	Point midpoint_12 = p1 + (p2 - p1) / 2.;
	Vector normal_12 = (n1 + n2) / 2.;

	Point p3 = mst.graph[v3].coords;
	Point p4 = mst.graph[v4].coords;
	Vector n3 = mst.graph[v3].normal;
	Vector n4 = mst.graph[v4].normal;
	Point midpoint_34 = p3 + (p4 - p3) / 2.;
	Vector normal_34 = (n3 + n4) / 2.;

	// On the plane of edge 12
	{
		bool isIntersecting = true;
		Vector edge1 = p1 - midpoint_12;
		Vector edge2 = p3 - midpoint_12;
		Vector edge3 = p4 - midpoint_12;
		Vector proj_edge1 = projected_vector(edge1, normal_12);
		Vector proj_edge2 = projected_vector(edge2, normal_12);
		Vector proj_edge3 = projected_vector(edge3, normal_12);
		Vector pro1 = CGAL::cross_product(proj_edge2, proj_edge1);
		Vector pro2 = CGAL::cross_product(proj_edge3, proj_edge1);
		if (pro1 * pro2 > 0)
			isIntersecting = false;
		if (isIntersecting) {
			edge1 = p3 - midpoint_34;
			edge2 = p1 - midpoint_34;
			edge3 = p2 - midpoint_34;
			proj_edge1 = projected_vector(edge1, normal_12);
			proj_edge2 = projected_vector(edge2, normal_12);
			proj_edge3 = projected_vector(edge3, normal_12);
			pro1 = CGAL::cross_product(proj_edge2, proj_edge1);
			pro2 = CGAL::cross_product(proj_edge3, proj_edge1);
			if (pro1 * pro2 > 0)
				isIntersecting = false;
		}
		if (isIntersecting)
			return true;
	}

	// On the plane of edge 34
	if(true){
		bool isIntersecting = true;
		Vector edge1 = p1 - midpoint_12;
		Vector edge2 = p3 - midpoint_12;
		Vector edge3 = p4 - midpoint_12;
		Vector proj_edge1 = projected_vector(edge1, normal_34);
		Vector proj_edge2 = projected_vector(edge2, normal_34);
		Vector proj_edge3 = projected_vector(edge3, normal_34);
		Vector pro1 = CGAL::cross_product(proj_edge2, proj_edge1);
		Vector pro2 = CGAL::cross_product(proj_edge3, proj_edge1);
		if (pro1 * pro2 > 0)
			isIntersecting = false;
		if (isIntersecting) {
			edge1 = p3 - midpoint_34;
			edge2 = p1 - midpoint_34;
			edge3 = p2 - midpoint_34;
			proj_edge1 = projected_vector(edge1, normal_34);
			proj_edge2 = projected_vector(edge2, normal_34);
			proj_edge3 = projected_vector(edge3, normal_34);
			pro1 = CGAL::cross_product(proj_edge2, proj_edge1);
			pro2 = CGAL::cross_product(proj_edge3, proj_edge1);
			if (pro1 * pro2 > 0)
				isIntersecting = false;
		}
		if (isIntersecting)
			return true;
	}
	return false;
}

/**
 * @brief Geometry check for connection
 *
 * @param mst: graph and vertex information
 * @param candidate: the edge to be examed
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 *
 * @return if the candidate pass the check
 */
bool geometry_check(m_Graph& mst, m_Edge& candidate, Tree& kdTree,
	Distance& tr_dist) {
	Vertex v1 = candidate.first;
	Vertex v2 = candidate.second;
	Point p1 = mst.graph[v1].coords;
	Point p2 = mst.graph[v2].coords;
	Vector n1 = mst.graph[v1].normal;
	Vector n2 = mst.graph[v2].normal;

	Vector mean_normal = (n1 + n2) / 2.;
	mean_normal = normalize_vector(mean_normal);

	Point search_center = p1 + (p2 - p1) / 2.;
	float radius = std::sqrt((p2-p1).squared_length())/2.;
	std::vector<int> neighbors;
	std::vector<float> distance;
	NN_search(-1, search_center, kdTree, tr_dist, float(radius * 3.), neighbors, distance, false);

	/*if ((v1 == 30045 && v2 == 69461) || v1 == 69461 && v2 == 30045)
		std::cout << "debug here" << std::endl;*/
	// Start neighbors check
	float query_radian1 = cal_radians_3d(p1 - search_center, mean_normal);
	float query_radian2 = cal_radians_3d(p2 - search_center, mean_normal);
	std::set<int> rejection_neighbor_set;
	for (int i = 0; i < neighbors.size(); i++) {
		if (neighbors[i] == v1 || neighbors[i] == v2)
			continue;
		if (mst.graph[neighbors[i]].normal * mean_normal > std::cos(60. / 180. * CGAL_PI))
			rejection_neighbor_set.insert(neighbors[i]);
	}
	for (int i = 0; i < neighbors.size(); i++) {
		if (rejection_neighbor_set.find(neighbors[i]) ==
			rejection_neighbor_set.end())
			continue;
		//if (neighbors[i] == v1 || neighbors[i] == v2)
		//	continue;
		Vertex rejection_neighbor = neighbors[i];
		Point rej_neighbor_pos = mst.graph[rejection_neighbor].coords;
		float min_radian, max_radian;

		for (auto& rej_neighbor_neighbor : mst.graph[rejection_neighbor].ordered_neighbors) {
			if (rejection_neighbor_set.find(rej_neighbor_neighbor.v) ==
				rejection_neighbor_set.end())
				continue;

			if (false) {
				min_radian = cal_radians_3d(rej_neighbor_pos - search_center, mean_normal);
				
				Point rej_neighbor_neighbor_pos = mst.graph[Vertex(rej_neighbor_neighbor.v)].coords;
				max_radian = cal_radians_3d(rej_neighbor_neighbor_pos - search_center,
					mean_normal);

				if (max_radian < min_radian) {
					std::swap(max_radian, min_radian);
				}
				if (max_radian - min_radian > CGAL_PI)
					std::swap(max_radian, min_radian);

				bool is_in_between = false;
				if (max_radian < min_radian &&
					(query_radian1 > min_radian || query_radian1 < max_radian))
					is_in_between = true;
				if (max_radian > min_radian &&
					(query_radian1 < max_radian && query_radian1 > min_radian))
					is_in_between = true;
				if (max_radian < min_radian &&
					(query_radian2 > min_radian || query_radian2 < max_radian))
					is_in_between = true;
				if (max_radian > min_radian &&
					(query_radian2 < max_radian && query_radian2 > min_radian))
					is_in_between = true;

				if (is_in_between) {
					Vector edge1 = p1 - rej_neighbor_pos;
					Vector edge2 = p2 - rej_neighbor_pos;
					Vector edge3 = rej_neighbor_neighbor_pos - rej_neighbor_pos;
					Vector proj_edge1 = projected_vector(edge1, mean_normal);
					Vector proj_edge2 = projected_vector(edge2, mean_normal);
					Vector proj_edge3 = projected_vector(edge3, mean_normal);
					Vector pro1 = CGAL::cross_product(proj_edge1, proj_edge3);
					Vector pro2 = CGAL::cross_product(proj_edge2, proj_edge3);
					if (pro1 * pro2 <= 0)
						return false;
				}
			}
			else {
				bool result = isIntersecting(mst, v1, v2, rejection_neighbor, rej_neighbor_neighbor.v);
				if (result)
					return false;
			}
		}
		rejection_neighbor_set.erase(neighbors[i]);
	}
	return true;
}

/**
 * @brief Topology check for connection and entrance of geometry check
 *
 * @param mst: graph and vertex information
 * @param candidate: the edge to be examed
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 *
 * @return if the candidate pass both checks
 */
bool Vanilla_check(m_Graph& mst, m_Edge& candidate, Tree& kdTree, 
	Distance& tr_dist) {
	Vertex neighbor = candidate.second;
	Vertex this_v = candidate.first;
	Vector this_normal = normalize_vector(mst.graph[this_v].normal);
	Vector neighbor_normal = normalize_vector(mst.graph[neighbor].normal);

	//if ((this_v == 297879 && neighbor == 298084) ||
	//	(this_v == 298084 && neighbor == 297879)) {
	//	/*Vertex temp_v = predecessor(mst, neighbor, this_v).v;
	//	Vertex last_v = neighbor;
	//	for (int i = 0; i < 20; i++) {
	//		Vertex temp = predecessor(mst, temp_v, last_v).v;
	//		last_v = temp_v;
	//		temp_v = temp;
	//		std::cout << temp_v << std::endl;
	//	}*/
	//	std::cout << "debug here" << std::endl;
	//}

	// Topology check
	if(true) {
		auto this_v_tree = predecessor(mst, this_v, neighbor).tree_id;
		auto neighbor_tree = predecessor(mst, neighbor, this_v).tree_id;

		if (!mst.etf.connected(this_v_tree, neighbor_tree)) {
			return false;
		}
	}

	return geometry_check(mst, candidate, kdTree, tr_dist);
	//return true;
}

/**
 * @brief Connect handle to raise the genus number
 * 
 * @param smoothed_v: smoothed vertices of the point cloud
 * @param mst: graph and vertex information
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param connected_handle_root: [OUT] log the connected handles
 * @param betti: [OUT] log betti number changes
 * @param k: number of kNN search
 * @param isEuclidean: if to use Euclidean distance
 * @param step_thresh: step threshold for shortest distance path early stop
 *
 * @return None
 */
void connect_handle(const std::vector<Point>& smoothed_v, Tree& KDTree, Distance& tr_dist,
	m_Graph& mst, std::vector<Vertex>& connected_handle_root, std::vector<int>& betti, int k, bool isEuclidean, int step_thresh) {

	std::vector<Vertex> imp_node;
	int num = 0;
	int edge_num = 0;
	// Collect vertices w/ an open angle larger than pi
	{
		for (int i = 0; i < boost::num_vertices(mst.graph); i++) {
			Vertex this_v = i;
			std::set<Neighbor>& neighbors = mst.graph[this_v].ordered_neighbors;
			float last_angle = (*(--neighbors.end())).angle;
			float this_angle;

			for (auto& neighbor : neighbors) {
				this_angle = neighbor.angle;
				float angle_diff = this_angle - last_angle;
				if (angle_diff < 0)
					angle_diff += 2 * CGAL_PI;
				if (angle_diff > CGAL_PI)
					imp_node.push_back(this_v);
				last_angle = this_angle;
			}
		}
	}

	std::vector<Vertex> connect_p;
	std::vector<Vertex> to_connect_p;
	std::vector<uint> tree_id;
	std::vector<uint> to_tree_id;
	for (auto& this_v : imp_node) {
		Point query = mst.graph[this_v].coords;
		Vector query_normal = mst.graph[this_v].normal;
		std::vector<int> neighbors;
		std::vector<float> dists;
		
		// Potential handle collection
		uint tree, to_tree;
		int validIdx = -1;

		kNN_search(this_v, smoothed_v[int(this_v)], KDTree, tr_dist, k, neighbors, dists);
		for (int i = 0; i < neighbors.size(); i++) {
			int neighbor = neighbors[i];
			m_Edge candidate(this_v, neighbor);
			if (boost::edge(this_v, neighbor, mst.graph).second)
				continue;
			tree = mst.etf.representative((predecessor(mst, this_v, neighbor).tree_id));
			to_tree = mst.etf.representative(predecessor(mst, neighbor, this_v).tree_id);
			if (geometry_check(mst, candidate, KDTree, tr_dist) && tree != to_tree) {
				validIdx = i;
				break;
			}
		}
		// TODO: Check if any tree shares root, and return corresponding edges

		if (validIdx != -1) {
			connect_p.push_back(this_v);
			to_connect_p.push_back(neighbors[validIdx]);
			tree_id.push_back(tree);
			to_tree_id.push_back(to_tree);
		}
	}

	// Select one handle
	std::map<std::string, std::vector<int>> face_connections;
	for (int i = 0; i < connect_p.size(); i++) {
		uint tree = tree_id[i];
		uint to_tree = to_tree_id[i];
		if (to_tree > tree)
			std::swap(tree, to_tree);
		std::string key = std::to_string(tree) + "+" + std::to_string(to_tree);
		if (face_connections.find(key) == face_connections.end())
			face_connections[key] = std::vector<int>{ i };
		else {
			face_connections[key].push_back(i);
		}
	}

	// Sort
	std::vector<m_face_pair> sorted_face;
	for (auto key = face_connections.begin(); key != face_connections.end(); key++) {
		int length = face_connections[key->first].size();
		sorted_face.push_back(m_face_pair(length, key->first));
	}
	std::sort(sorted_face.begin(), sorted_face.end(), face_comparator);
	for (int i = 0; i < sorted_face.size(); i++) {
		std::string key = sorted_face[i].second;
		std::vector<int> idx_vec = face_connections[key];
		if (idx_vec.size() <= 5)
			break;
		if (mst.exp_genus >= 0 && num >= mst.exp_genus)
			break;
		Point query;
		Vertex connected_neighbor, this_v;
		Edge added_edge;
		bool isFind = false;
		for (int idx : idx_vec) {
			this_v = connect_p[idx];
			query = mst.graph[this_v].coords;
			connected_neighbor = to_connect_p[idx];
			std::vector<Vertex> path;
			int steps = find_shortest_path(mst, this_v, connected_neighbor, step_thresh, path);
			if (steps < 0) {
			//if(steps >= 9){
				//std::cout << "This is connected" << std::endl;
				isFind = true;
				m_Edge candidate(this_v, connected_neighbor);
				if (geometry_check(mst, candidate, KDTree, tr_dist)) {
					Vector edge = query - mst.graph[connected_neighbor].coords;
					float Euclidean_dist = norm(edge);
					float projection_dist = cal_proj_dist(edge, mst.graph[this_v].normal,
						mst.graph[connected_neighbor].normal);

					if (isEuclidean) {
						added_edge = mst.add_edge(this_v, connected_neighbor, Euclidean_dist, true);
						connected_handle_root.push_back(this_v);
						connected_handle_root.push_back(connected_neighbor);
					}
					else {
						added_edge = mst.add_edge(this_v, connected_neighbor, projection_dist, true);
						connected_handle_root.push_back(this_v);
						connected_handle_root.push_back(connected_neighbor);
					}

					bettiNum_1++;
					edge_num++;
					betti.push_back(bettiNum_1);
					//fs::path out_edge_path("C:/Projects_output/letters/edge_" + std::to_string(edge_num) + ".obj");
					//export_continuous_edges(mst, path, out_edge_path);
				}
			}
		}
		if (isFind) {
			num++;
		}
	}

	std::cout << "Handle Connection done :)" << std::endl;
	std::cout << std::to_string(num) << " pairs of faces are connected." << std::endl;
	std::cout << std::to_string(edge_num) << " edges are connected." << std::endl;
	return;
}
