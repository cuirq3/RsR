#include "triangulation.h"

void sortVertices(int& a, int& b, int& c) {
	if (a > b) {
		std::swap(a, b);
	}
	if (b > c) {
		std::swap(b, c);
	}
	if (a > b) {
		std::swap(a, b);
	}
}

bool explore(m_Graph& G, int i, m_priority_queue& queue,
	std::unordered_set<std::string>& faces_in_queue, float avg_edge_length
	, std::vector<float>& length_thresh) {
	Vertex v_i = i;
	bool isFound = false;
	//if (v_i == 1439845)
	//	std::cout << "debug here" << std::endl;
	for (auto& neighbor : G.graph[v_i].ordered_neighbors)
	{
		Vertex v_u = neighbor.v;
		Vertex v_w = successor(G, v_i, v_u).v;

		//if ((v_u == 1751375 && v_w == 1888150) ||
		//	(v_w == 1751375 && v_u == 1888150))
		//	std::cout << "debug here" << std::endl;

		//if ((v_u == 33020 && v_w == 30023) || (v_w == 33020 && v_u == 30023)) {
		//	std::cout << "debug here" << std::endl;
		//}

		Point w_pos = G.graph[v_w].coords;
		Point u_pos = G.graph[v_u].coords;
		Point i_pos = G.graph[v_i].coords;
		Vector i_normal = G.graph[v_i].normal;
		Vector u_normal = G.graph[v_u].normal;
		Vector w_normal = G.graph[v_w].normal;
		float angle = cal_radians_3d(w_pos - i_pos, i_normal,
			u_pos - i_pos);
		bool isLargerThanPi = angle < CGAL_PI;
		std::vector<Vertex> face_vector{ v_i, v_u, v_w };
		if (v_u != v_w && isLargerThanPi) {
			if (!boost::edge(v_u, v_w, G.graph).second) {
				float score = (G.graph[v_u].coords - G.graph[v_w].coords).squared_length();
				if (!G.isEuclidean) {
					score = cal_proj_dist(G.graph[v_u].coords - G.graph[v_w].coords,
						u_normal, w_normal);
				}
				else
					score = std::sqrt(score);
				if (score > length_thresh[v_u] || score > length_thresh[v_w])
					continue;
				if (score >= 0){
					std::pair<std::vector<Vertex>, float> queue_item(face_vector, score);
					queue.push(queue_item);
					isFound = true;
				}
			}
		}
	}
	
	return isFound;
}

bool check_face_overlap(m_Graph& G, std::vector<Vertex>& face,
	const Tree& KDTree, const Distance& tr_dist, int k) {
	Vertex v_i = face[0];
	Vertex v_u = face[1];
	Vertex v_w = face[2];
	Point pos_i = G.graph[v_i].coords;
	Point pos_u = G.graph[v_u].coords;
	Point pos_w = G.graph[v_w].coords;
	Vector normal_i = G.graph[v_i].normal;
	Vector normal_u = G.graph[v_u].normal;
	Vector normal_w = G.graph[v_w].normal;

	Point centroid((pos_i.x() + pos_u.x() + pos_w.x()) / 3,
		(pos_i.y() + pos_u.y() + pos_w.y()) / 3,
		(pos_i.z() + pos_u.z() + pos_w.z()) / 3);
	std::vector<Point> triangle_pos{ pos_i, pos_u, pos_w };

	std::vector<int> neighbors;
	std::vector<float> dists;
	kNN_search(v_i, centroid, KDTree, tr_dist, k, neighbors, dists, false);

	//Vector face_normal = CGAL::cross_product(pos_i - pos_u, pos_w - pos_u);
	Vector face_normal = triangle_mean_normal(normal_i, normal_u, normal_w);
	face_normal = normalize_vector(face_normal);
	if (normal_i * face_normal < 0)
		face_normal *= -1;

	for (int idx = 0; idx < neighbors.size(); idx++) {
		int neighbor = neighbors[idx];
		Vertex v_neighbor = neighbor;
		Point neighbor_pos = G.graph[v_neighbor].coords;
		Vector neighbor_normal = G.graph[v_neighbor].normal;
		if (neighbor == v_i || neighbor == v_u ||
			neighbor == v_w || neighbor_normal * face_normal < 0)
			continue;
		//std::vector<int> neighbors_neighbor = G.graph[v_neighbor].ordered_neighbors;

		//// Ignore closed vertex
		//if (neighbors_neighbor.size() > 2) {
		//	bool isClosed = true;
		//	for (int j = 0; j < neighbors_neighbor.size(); j++) {
		//		if (!boost::edge(get_vertex(G, neighbors_neighbor[j]), get_vertex(G, neighbors_neighbor[python_mod(j + 1, neighbors_neighbor.size())]), G.graph).second) {
		//			isClosed = false;
		//			break;
		//		}
		//	}
		//	if (isClosed)
		//		continue;
		//}

		if (point_in_triangle(neighbor_pos, triangle_pos, face_normal))
			return true;
	}
	return false;
}

bool check_branch_validity(m_Graph& G, Vertex root, Vertex branch1, Vertex branch2) {
	Point pos_i = G.graph[root].coords;
	Point pos_u = G.graph[branch1].coords;
	Point pos_w = G.graph[branch2].coords;
	Vector normal_i = G.graph[root].normal;
	Vector normal_u = G.graph[branch1].normal;
	Vector normal_w = G.graph[branch2].normal;

	// Option 1
	std::vector<Point> triangle_pos{ pos_i, pos_u, pos_w };
	//Vector face_normal = CGAL::cross_product(pos_i - pos_u, pos_w - pos_u);
	Vector face_normal = triangle_mean_normal(normal_i, normal_u, normal_w);

	float angle_thresh = 0. / 180. * CGAL_PI;

	/*if(!G.isFinalize){
		for (auto neighbor : G.graph[branch1].ordered_neighbors) {
			if (neighbor == root)
				continue;
			Point pos = G.graph[neighbor].coords;
			if (point_in_triangle(pos, triangle_pos, face_normal))
				return false;
		}

		for (auto neighbor : G.graph[branch2].ordered_neighbors) {
			if (neighbor == root)
				continue;
			Point pos = G.graph[neighbor].coords;
			if (point_in_triangle(pos, triangle_pos, face_normal))
				return false;
		}
	}*/

	// Check u's RS validity
	float this_radian = cal_radians_3d(pos_w - pos_u, normal_u);
    auto former = predecessor(G,branch1,branch2);
    auto next = successor(G,branch1,branch2);
    if(G.isFinalize){
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * CGAL_PI;
            if (diff < CGAL_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * CGAL_PI;
            if (diff < CGAL_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    }
    else {
        float diff = next.angle - this_radian;
        if(diff < 0)
            diff += 2 * CGAL_PI;
        if (next.v != root || diff > CGAL_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = this_radian - former.angle;
        if (diff_angle_thresh < 0)
            diff_angle_thresh += CGAL_PI * 2.;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }

	//Check w
	this_radian = cal_radians_3d(pos_u - pos_w, normal_w);
    former = predecessor(G,branch2,branch1);
    next = successor(G,branch2,branch1);
    if(G.isFinalize){
        bool isValid = false;
        if (next.v == root) {
            float diff = next.angle - this_radian;
            if (diff < 0)
                diff += 2 * CGAL_PI;
            if (diff < CGAL_PI)
                isValid = true;
        }
        if (former.v == root) {
            float diff = -former.angle + this_radian;
            if (diff < 0)
                diff += 2 * CGAL_PI;
            if (diff < CGAL_PI)
                isValid = true;
        }
        if (!isValid)
            return false;
    }
    else {
        float diff = -former.angle + this_radian ;
        if (diff < 0)
            diff += 2 * CGAL_PI;
        if (former.v != root || diff > CGAL_PI) {
            return false;
        }
    }

    // Thresh on angle
    {
        float diff_angle_thresh = -this_radian + next.angle + 2. * CGAL_PI;
        if (diff_angle_thresh < angle_thresh)
            return false;
    }


	//// Option 2
	//// Find former neighbor and latter neighbor
	//float radian = cal_radians_3d(pos_w - pos_u, normal_u);
	//int former, latter;
	//std::vector<float> ordered_radians = G.graph[branch1].ordered_radians;
	//for (int idx = 0; idx < ordered_radians.size(); idx++){
	//	if (radian < ordered_radians[idx]) {
	//		latter = idx;
	//		former = python_mod(idx - 1, ordered_radians.size());
	//		break;
	//	}
	//}
	//if (radian >= ordered_radians[ordered_radians.size() - 1]) {
	//	former = ordered_radians.size() - 1;
	//	latter = 0;
	//}
	//Vertex v_former = get_vertex(G, G.graph[branch1].ordered_neighbors[former]);
	//Vertex v_latter = get_vertex(G, G.graph[branch1].ordered_neighbors[latter]);
	//Point pos_latter = G.graph[v_latter].coords;
	//Point pos_former = G.graph[v_former].coords;
	//// Check
	//if (boost::edge(v_former, v_latter, G.graph).second) {
	//	if (cal_radians_3d(pos_latter - pos_u, normal_u, pos_former - pos_u) < CGAL_PI) {
	//		/*if (v_former != root && v_latter != root)
	//			return false;
	//		std::vector<Point> triangle_pos{ pos_u, pos_former, pos_latter };
	//		Vector face_normal = CGAL::cross_product(pos_u - pos_former, pos_latter - pos_former);
	//		face_normal = normalize_vector(face_normal);
	//		if (!point_in_triangle(pos_w, triangle_pos, face_normal) && ordered_radians.size() > 2)
	//			return false;*/
	//		return false;
	//	}
	//}
	return true;
}

bool check_validity(m_Graph& G, std::pair<std::vector<Vertex>, float>& item,
	const Tree& KDTree, const Distance& tr_dist, bool isFaceloop, bool isFinalize) {
	int i = item.first[0];
	int u = item.first[1];
	int w = item.first[2];
	Vertex v_i = i;
	Vertex v_u = u;
	Vertex v_w = w;
	Point pos_i = G.graph[v_i].coords;
	Point pos_u = G.graph[v_u].coords;
	Point pos_w = G.graph[v_w].coords;
	Vector normal_i = G.graph[v_i].normal;
	Vector normal_u = G.graph[v_u].normal;
	Vector normal_w = G.graph[v_w].normal;

	//if (i == 388237)
	//	std::cout << u << " " << w << std::endl;

	if (boost::edge(v_u, v_w, G.graph).second)
		return false;

    // Non-manifold edge check
	if (G.graph[boost::edge(v_i, v_u, G.graph).first].appear_time == 2 ||
		G.graph[boost::edge(v_i, v_w, G.graph).first].appear_time == 2)
		return false;

	// Check this rotation system
	bool isValid = (successor(G, v_i, v_u).v == v_w);
	float angle = cal_radians_3d(pos_w - pos_i, normal_i, pos_u - pos_i);
	if (angle > CGAL_PI)
		isValid = false;

	if (!isValid)
		return false;

	if (!isFinalize) {
		// Check rotation system's validity of branch nodes
		if (!check_branch_validity(G, v_i, v_u, v_w))
			return false;
	}

	// Check face overlap
	/*if (check_face_overlap(G, item.first, KDTree, tr_dist))
		return false;*/

	return true;
}

void add_face(m_Graph& G, std::vector<Vertex>& item,
	std::vector<Face>& faces) {
	Vertex v_i = item[0];
	Vertex v_u = item[1];
	Vertex v_w = item[2];

	// Maintain face exist for detecting holes
	{
        get_neighbor_info(G,v_i,v_u).faceExist = true;
        get_neighbor_info(G,v_u,v_w).faceExist = true;
        get_neighbor_info(G,v_w,v_i).faceExist = true;
	}

	G.graph[boost::edge(v_u, v_w, G.graph).first].appear_time += 1;
	G.graph[boost::edge(v_i, v_w, G.graph).first].appear_time += 1;
	G.graph[boost::edge(v_u, v_i, G.graph).first].appear_time += 1;
	faces.push_back(Face(v_i, v_u, v_w));
	bettiNum_1--;
	return;
}

void checkAndForce(Vertex v_u, Vertex v_w, m_Graph& G, m_priority_queue& queue, 
	std::vector<float>& length_thresh) {
	std::vector<Vertex> check_v{ v_u,v_w };
	for (int i = 0; i < check_v.size();i++) {
		Vertex v_i = check_v[i];
		std::vector<int> connection_count(G.graph[v_i].ordered_neighbors.size(), 0);
		std::vector<Vertex> neighbor_ids;
		std::set<Vertex> neighbor_ids_set;

		for (auto& neighbor : G.graph[v_i].ordered_neighbors) {
			neighbor_ids_set.insert(neighbor.v);
			neighbor_ids.push_back(neighbor.v);
		}

		int idx = 0;
		for (auto& neighbor : G.graph[v_i].ordered_neighbors) {
			for (auto& neighbor_neighbor : G.graph[neighbor.v].ordered_neighbors) {
				if (neighbor_ids_set.find(neighbor_neighbor.v) != neighbor_ids_set.end()) {
					connection_count[idx] += 1;
				}
			}
			idx++;
		}

		std::vector<Vertex> not_full;
		for (int i = 0; i < connection_count.size(); i++) {
			if (connection_count[i] < 2) {
				not_full.push_back(i);
			}
		}

		if (not_full.size() == 2) {
			Vertex v_u = neighbor_ids[not_full[0]];
			Vertex v_w = neighbor_ids[not_full[1]];
			Vector u_normal = G.graph[v_u].normal;
			Vector w_normal = G.graph[v_w].normal;
			float score = (G.graph[v_u].coords - G.graph[v_w].coords).squared_length();
			if (!G.isEuclidean) {
				score = cal_proj_dist(G.graph[v_u].coords - G.graph[v_w].coords,
					u_normal, w_normal);
			}
			else
				score = std::sqrt(score);
			if (score > length_thresh[v_u] || score > length_thresh[v_w])
				return;
			std::vector<Vertex> face_vector{ v_i, neighbor_ids[not_full[0]], neighbor_ids[not_full[1]] };
			std::pair<std::vector<Vertex>, float> queue_item(face_vector, -1);
			queue.push(queue_item);
			break;
		}
	}
	return;
}

void triangulate(std::vector<Face>& faces, m_Graph& G,
	const Tree& KDTree, const Distance& tr_dist, bool isFaceLoop, bool isEuclidean
	, std::vector<float>& length_thresh, std::vector<Vertex>& connected_handle_root, std::vector<int>& betti, bool isFinalize) {

	std::unordered_set<std::string> faces_in_queue;
	std::unordered_set<Vertex> to_visit;

	m_priority_queue queue;

	float avg_edge_length = G.total_edge_length / boost::num_edges(G.graph);

	// Init priority queue
	for (int i = 0; i < boost::num_vertices(G.graph); i++) {
		to_visit.insert(i);
	}

	for (int i = 0; i < connected_handle_root.size(); i++) {
		bool result = explore(G, connected_handle_root[i], queue, faces_in_queue, avg_edge_length, length_thresh);
		to_visit.erase(connected_handle_root[i]);
	}

	std::cout << "Global init done :)" << std::endl;

	int loop_time = 0;

	while (to_visit.size() > 0) {
		while (!queue.empty()) {
			loop_time += 1;
			std::pair<std::vector<Vertex>, float> item = queue.top();
			queue.pop();

			if (item.second >= 0) {
				// Validity check
				bool isValid = check_validity(G, item, KDTree, tr_dist, isFaceLoop, isFinalize);
				if (!isValid)
					continue;
			}

			// Add the edge
			Vertex v_i = item.first[0];
			Vertex v_u = item.first[1];
			Vertex v_w = item.first[2];
			Point pos_i = G.graph[v_i].coords;
			Point pos_u = G.graph[v_u].coords;
			Point pos_w = G.graph[v_w].coords;
			Vector normal_i = G.graph[v_i].normal;
			Vector normal_u = G.graph[v_u].normal;
			Vector normal_w = G.graph[v_w].normal;

			float dist = norm(pos_u - pos_w);
			Vector edge = pos_u - pos_w;
			float Euclidean_dist = norm(edge);
			float projection_dist = cal_proj_dist(edge,
				G.graph[v_u].normal, G.graph[v_w].normal);

			if (!boost::edge(v_u, v_w, G.graph).second) {
				Edge added_edge;
				if (isEuclidean)
					added_edge = G.add_edge(v_u, v_w, Euclidean_dist);
				else
					added_edge = G.add_edge(v_u, v_w, projection_dist);

				avg_edge_length = G.total_edge_length / boost::num_edges(G.graph);

				bettiNum_1++;

				add_face(G, item.first, faces);

				betti.push_back(bettiNum_1);
			}
			else
				continue;

			// Deal with incident triangles
			{
				std::vector<Vertex> share_neighbors;
				find_common_neighbor(v_u, v_w, share_neighbors, G);
				for (int idx = 0; idx < share_neighbors.size(); idx++) {
					Vertex incident_root = share_neighbors[idx];
					if (incident_root == v_i)
						continue;
					std::vector<Vertex> face{ incident_root,v_w,v_u };

					// Non-manifold edge check
					int time1 = G.graph[boost::edge(incident_root, v_u, G.graph).first].appear_time;
					int time2 = G.graph[boost::edge(incident_root, v_w, G.graph).first].appear_time;
					int time3 = G.graph[boost::edge(v_u, v_w, G.graph).first].appear_time;
					if (time1 == 2 || time2 == 2 || time3 == 2)
						continue;

					add_face(G, face, faces);
				}
			}

			to_visit.erase(v_u);
			to_visit.erase(v_w);

			// Explore and sanity check
			bool isFound = false;
			bool result = explore(G, v_u, queue, faces_in_queue, avg_edge_length, length_thresh);
			isFound = isFound || result;
			result = explore(G, v_w, queue, faces_in_queue, avg_edge_length, length_thresh);
			isFound = isFound || result;

			if (isFinalize) {
				if ((!isFound)) {
					checkAndForce(v_u, v_w, G, queue, length_thresh);
				}
			}
		}
		
		if (!to_visit.empty()) {
			Vertex pick = *to_visit.begin();
			to_visit.erase(pick);
			bool result = explore(G, pick, queue, faces_in_queue, avg_edge_length, length_thresh);
		}
	}
}


bool routine_check(m_Graph& mst, std::vector<Vertex>& triangle) {

	Vertex v1 = triangle[0];
	Vertex v2 = triangle[1];
	Vertex v3 = triangle[2];
	Point p1 = mst.graph[v1].coords;
	Point p2 = mst.graph[v2].coords;
	Point p3 = mst.graph[v3].coords;
	Vector n1 = mst.graph[v1].normal;
	Vector n2 = mst.graph[v2].normal;
	Vector n3 = mst.graph[v3].normal;

	//Vector face_normal = normalize_vector(CGAL::cross_product(p2 - p1, p3 - p1));
	//bool isValid = (n1 * face_normal * (n2 * face_normal) < 0 ||
	//		n1 * face_normal * (n3 * face_normal) < 0);
	
	{
		float len_ui = norm(p1 - p2);
		float len_wi = norm(p3 - p2);
		float len_uw = norm(p1 - p3);

		float max_value = std::acos(boost::algorithm::clamp(
			(p3 - p2) * (p1 - p2) / len_ui /
			len_wi, -1, 1));
		float radian = std::acos(boost::algorithm::clamp(
			(p2 - p1) * (p3 - p1) / len_ui /
			len_uw, -1, 1));
		if (radian > max_value)
			max_value = radian;
		radian = std::acos(boost::algorithm::clamp(
			(p1 - p3) * (p2 - p3) / len_uw /
			len_wi, -1, 1));
		if (radian > max_value)
			max_value = radian;
		if (max_value > 175. / 180. * CGAL_PI)
			return true;
	}


	if (mst.graph[boost::edge(v1, v3, mst.graph).first].appear_time == 2 ||
		mst.graph[boost::edge(v2, v3, mst.graph).first].appear_time == 2)
		return true;

	return false;
}

bool register_face(m_Graph& mst, Vertex v1, Vertex v2, std::vector<Face>& faces,
	Tree& KDTree, Distance& tr_dist, float edge_length) {

	Vector v1_n = mst.graph[v1].normal;
	Vector v2_n = mst.graph[v2].normal;
	Point p1 = mst.graph[v1].coords;
	Point p2 = mst.graph[v2].coords;

	if (boost::edge(v1, v2, mst.graph).second)
		return false;

	std::vector<Vertex> share_neighbors;
	find_common_neighbor(v1, v2, share_neighbors, mst);
	if (share_neighbors.size() == 0) {
		Edge added_edge = mst.add_edge(v1, v2,
			edge_length);
		maintain_face_loop(mst, added_edge);
		return true;
	}
	//if ((v1 == 30045 && v2 == 69461) || v1 == 69461 && v2 == 30045)
	//	std::cout << "debug here" << std::endl;

	auto possible_root1 = predecessor(mst, v1, v2).v;
	float angle1 = cal_radians_3d(p1 - mst.graph[possible_root1].coords,
		mst.graph[possible_root1].normal, p2 - mst.graph[possible_root1].coords);
	auto possible_root2 = predecessor(mst, v2, v1).v;
	float angle2 = cal_radians_3d(p2 - mst.graph[possible_root2].coords,
		mst.graph[possible_root2].normal, p1 - mst.graph[possible_root2].coords);

	bool isValid = true;
	std::vector<Face> temp;
	for (auto v3 : share_neighbors) {
		//if (v3 == 388212)
		//	std::cout << "debug" << std::endl;
		std::vector<Vertex> triangle{ v1,v2,v3 };
		if (v3 == possible_root1 && angle1 < CGAL_PI) {
			if (routine_check(mst, triangle)) {
				isValid = false;
				break;
			}
			if (successor(mst, v2, v1).v != v3) {
				isValid = false;
				break;
			}
			if (successor(mst, possible_root1, v2).v != v1) {
				isValid = false;
				break;
			}
			/*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
				isValid = false;
				break;
			}*/
			temp.push_back(Face(v1, v3, v2));
		}

		if (v3 == possible_root2 && angle2 < CGAL_PI) {
			if (routine_check(mst, triangle)) {
				isValid = false;
				break;
			}
			if (successor(mst, v1, v2).v != v3) {
				isValid = false;
				break;
			}
			if (successor(mst, possible_root2, v1).v != v2) {
				isValid = false;
				break;
			}
			/*if (check_face_overlap(mst, triangle, KDTree, tr_dist)) {
				isValid = false;
				break;
			}*/
			temp.push_back(Face(v1, v2, v3));
		}
	}
	
	if (temp.size() == 0)
		isValid = false;

	if (isValid) {
		Edge added_edge = mst.add_edge(v1, v2,
			edge_length);
		maintain_face_loop(mst, added_edge);
		for (auto& face : temp) {
			std::vector<Vertex> triangle{ face.ids[0],face.ids[1],face.ids[2] };
			add_face(mst, triangle, faces);
		}
	}
	
	return isValid;
}