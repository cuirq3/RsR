#include"normal.h"

/**
 * @brief Estimate normal using CGAL built-in function
 *
 * @param vertices: vertices of the point cloud
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param isGTNormal: boolean determining if the program is using provided normal
 * @param normals: [IN/OUT] normal of the point cloud
 * @param zero_normal_id: [IN/OUT] list of indices of invalid normal values ([0,0,0])
 * @param diagonal_length: [OUT] the diagonal length of the point cloud
 * 
 * @return None
 */
void estimate_normal(const std::vector<Point>& vertices,
	const Tree& kdTree, const Distance& tr_dist, bool isGTNormal,
	std::vector<Vector>& normals, std::vector<int>& zero_normal_id,
	float& diagonal_length) {
	std::vector<PointWithNormal> pointsWithProperties;
	int idx = 0;
	// Data type transfer & Cal diagonal size
	std::vector<float> min{ INFINITY, INFINITY, INFINITY },
		max{ -INFINITY, -INFINITY, -INFINITY };
	for (auto& point : vertices) {
		if (isGTNormal) {
			Vector normal = normals[idx];
			if (norm(normal) == 0) {
				zero_normal_id.push_back(idx);
			}
			else {
				normals[idx] = normalize_vector(normals[idx]);
			}
		}
		else
			pointsWithProperties.push_back(PointWithNormal(point));
		if (point.x() < min[0])
			min[0] = point.x();
		if (point.y() < min[1])
			min[1] = point.y();
		if (point.z() < min[2])
			min[2] = point.z();
		if (point.x() > max[0])
			max[0] = point.x();
		if (point.y() > max[1])
			max[1] = point.y();
		if (point.z() > max[2])
			max[2] = point.z();
		idx++;
	}
	Point min_p(min[0], min[1], min[2]), max_p(max[0], max[1], max[2]);
	diagonal_length = norm(max_p - min_p);

	if (!isGTNormal) {
		int neighbor_num = std::max<int>(int(vertices.size() / 2000.), 192);
		// Estimate normal
		CGAL::pca_estimate_normals<Concurrency_tag>
			(pointsWithProperties, neighbor_num,
				CGAL::parameters::point_map(CGAL::Nth_of_tuple_property_map<0, PointWithNormal>())
				.normal_map(CGAL::Nth_of_tuple_property_map<1, PointWithNormal>()));

		normals.clear();

		for (const auto& vertex_w_normal : pointsWithProperties) {
			if (norm(vertex_w_normal.get<1>()) == 0.
				|| !isfinite(norm(vertex_w_normal.get<1>())))
				std::cout << "error" << std::endl;
			normals.push_back(normalize_vector(vertex_w_normal.get<1>()));
		}
	}
	return;
}

/**
 * @brief Fix invalid normal value - [0, 0, 0] - by re-estimating it with its neighbors
 *
 * @param zero_normal_id: list of indices of invalid normal values ([0,0,0])
 * @param kdTree: kd-tree for knn query
 * @param tr_dist: distance container
 * @param vertices: vertices of the point cloud
 * @param normals: [OUT] normal of the point cloud
 *
 * @return None
 */
void replace_zero_normal(const std::vector<int>& zero_normal_id,
	const Tree& kdTree, const Distance& tr_dist, const std::vector<Point>& vertices,
	std::vector<Vector>& normals) {
	for (int id : zero_normal_id) {
		Point query = vertices[id];
		std::vector<int> neighbors;
		std::vector<float> neighbor_dist;
		kNN_search(id, query, kdTree, tr_dist, 11, neighbors, neighbor_dist);
		Vector normal;
		est_normal_SVD(neighbors, vertices, normal);
		int vote = 0;
		for (int neighbor : neighbors) {
			Vector neighbor_normal = Vector(normals[neighbor][0], normals[neighbor][1], normals[neighbor][2]);
			if (neighbor_normal.squared_length() != 0) {
				if (neighbor_normal * normal > 0)
					vote += 1;
				else
					vote -= 1;
			}
		}
		if (vote < 0)
			normal *= -1;
		normals[id] = normal;
	}
}

/**
 * @brief Determine the normal orientation
 *
 * @param G_angle graph whose edges have angle-based weight
 * @param normals: [OUT] normal of the point cloud with orientation corrected
 *
 * @return None
 */
void correct_normal_orientation(s_Graph& G_angle, std::vector<Vector>& normals) {
	std::vector<bool> visited_vertex(boost::num_vertices(G_angle), false);
	// Start from vertex 0
	std::queue<int> to_visit;
	to_visit.push(0);
	while (!to_visit.empty()) {
		int node_id = to_visit.front();
		to_visit.pop();
		visited_vertex[node_id] = true;
		Vertex node = node_id;
		Vector this_normal = normals[node_id];
		auto neighbours = boost::adjacent_vertices(node, G_angle);
		for (auto vd : make_iterator_range(neighbours)) {
			if (!visited_vertex[int(vd)]) {
				to_visit.push(int(vd));
				Vector neighbor_normal = normals[int(vd)];
				if (this_normal * neighbor_normal < 0) {
					normals[int(vd)] = -normals[int(vd)];
				}
			}
		}
	}
}