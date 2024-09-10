#include"RS.h"

typedef std::pair<float, Vertex> mypair;
bool comparator(const mypair& l, const mypair& r) {
	return l.first < r.first;
}

/*
void init_rotation_system(m_Graph& out_graph) {
	vertex_iter vi = out_graph.vi_start;
	for (int i = 0; i < boost::num_vertices(out_graph.graph); i++)
	{
		auto neighbours = boost::adjacent_vertices(i, out_graph.graph);
		std::vector<mypair> neighbour_infos;

		Vector normal = out_graph.graph[*vi].normal;
		Point root_coords = out_graph.graph[*vi].coords;
		Vector ref_vec;
		calculate_ref_vec(normal, ref_vec);

		int idx = 0;
		for (auto vd : make_iterator_range(neighbours)) {
			if (out_graph.graph[vd].id == out_graph.graph[*vi].id)
				continue;
			Point neighbor_pos = out_graph.graph[vd].coords;
			Vector branch_vec = neighbor_pos - root_coords;
			neighbour_infos.push_back(mypair(cal_radians_3d(branch_vec,
				normal, ref_vec), vd));
			idx++;
		}

		std::sort(neighbour_infos.begin(), neighbour_infos.end());
		// Modify vertex property
		for (auto info : neighbour_infos) {
			out_graph.graph[*vi].ordered_radians.push_back(info.first);
			out_graph.graph[*vi].ordered_neighbors.push_back(info.second);
			out_graph.graph[*vi].faceExist.push_back(false);
		}
		vi++;
	}
	return;
}
*/
/*
void maintain_rotation_system(m_Graph& mst, Edge& added_edge) {
	Vertex source_v = boost::source(added_edge, mst.graph);
	Vertex target_v = boost::target(added_edge, mst.graph);

	Point pos_s = mst.graph[source_v].coords;
	Vector normal_s = mst.graph[source_v].normal;
	Point pos_t = mst.graph[target_v].coords;
	Vector normal_t = mst.graph[target_v].normal;

	float radian_s_t = cal_radians_3d(pos_s - pos_t, normal_t);
	float radian_t_s = cal_radians_3d(pos_t - pos_s, normal_s);

	// Insert t to s
	auto radian_iter_start = mst.graph[source_v].ordered_radians.begin();
	auto radian_iter_end = mst.graph[source_v].ordered_radians.end();
	int idx = -1;
	for (auto radian_iter = radian_iter_start; radian_iter != radian_iter_end; radian_iter++) {
		if (radian_t_s < *radian_iter) {
			idx = radian_iter - radian_iter_start;
			break;
		}
	}
	if (idx < 0) {
		mst.graph[source_v].ordered_radians.push_back(radian_t_s);
		mst.graph[source_v].ordered_neighbors.push_back(int(target_v));
		mst.graph[source_v].faceExist.push_back(false);
	}
	else {
		mst.graph[source_v].ordered_radians.insert(radian_iter_start + idx, radian_t_s);
		auto neighbor_iter_start = mst.graph[source_v].ordered_neighbors.begin();
		mst.graph[source_v].ordered_neighbors.insert(neighbor_iter_start + idx, int(target_v));
		auto faceExist_iter_start = mst.graph[source_v].faceExist.begin();
		mst.graph[source_v].faceExist.insert(faceExist_iter_start + idx, false);
	}

	// Insert s to t
	radian_iter_start = mst.graph[target_v].ordered_radians.begin();
	radian_iter_end = mst.graph[target_v].ordered_radians.end();
	idx = -1;
	for (auto radian_iter = radian_iter_start; radian_iter != radian_iter_end; radian_iter++) {
		if (radian_s_t < *radian_iter) {
			idx = radian_iter - radian_iter_start;
			break;
		}
	}
	if (idx < 0) {
		mst.graph[target_v].ordered_radians.push_back(radian_s_t);
		mst.graph[target_v].ordered_neighbors.push_back(int(source_v));
		mst.graph[target_v].faceExist.push_back(false);
	}
	else {
		mst.graph[target_v].ordered_radians.insert(radian_iter_start + idx, radian_s_t);
		auto neighbor_iter_start = mst.graph[target_v].ordered_neighbors.begin();
		mst.graph[target_v].ordered_neighbors.insert(neighbor_iter_start + idx, int(source_v));
		auto faceExist_iter_start = mst.graph[target_v].faceExist.begin();
		mst.graph[target_v].faceExist.insert(faceExist_iter_start + idx, false);
	}
	return;
}

void find_radians_id(const m_Graph& G, const Vertex& root,
	const Vertex& branch, int& next_id, int& last_id, float& this_radian) {
	Vector root_normal = G.graph[root].normal;
	float radian_branch2Root = cal_radians_3d(G.graph[branch].coords - G.graph[root].coords, root_normal);
	std::vector<float> ordered_radians = G.graph[root].ordered_radians;
	next_id = 0;
	for (auto radian : ordered_radians) {
		if (radian_branch2Root < radian)
			break;
		next_id++;
	}
	next_id = next_id % ordered_radians.size();
	last_id = python_mod(next_id - 1, ordered_radians.size());
	this_radian = radian_branch2Root;
	return;
}

void find_radians_id_no_clamp(const m_Graph& G, const Vertex& root,
                     const Vertex& branch, int& next_id, int& last_id, float& this_radian) {
    Vector root_normal = G.graph[root].normal;
    float radian_branch2Root = cal_radians_3d(G.graph[branch].coords - G.graph[root].coords, root_normal);
    std::vector<float> ordered_radians = G.graph[root].ordered_radians;
    next_id = 0;
    for (auto radian : ordered_radians) {
        if (radian_branch2Root <= radian)
            break;
        next_id++;
    }
    next_id = next_id;
    last_id = python_mod(next_id - 1, ordered_radians.size());
    this_radian = radian_branch2Root;
    return;
}
*/
const Neighbor& get_neighbor_info(const m_Graph& g, const Vertex& root, const Vertex& branch){
    auto& u = g.graph[root];
    auto& v = g.graph[branch];
    auto iter = u.ordered_neighbors.lower_bound({u,v,static_cast<uint>(branch)});
    return (*iter);
}

const Neighbor& successor(const m_Graph &g, const Vertex &root, const Vertex &branch) {
    auto& u = g.graph[root];
    auto& v = g.graph[branch];
    auto iter = u.ordered_neighbors.upper_bound(Neighbor(u,v,branch));
    if(iter == u.ordered_neighbors.end()) iter = u.ordered_neighbors.begin(); // Wrap around
    return (*iter); // This is honestly not good practice - ONLY modification of tree_id
}

const Neighbor& predecessor(const m_Graph& g, const Vertex& root, const Vertex& branch) {
	auto& u = g.graph[root];
	auto& v = g.graph[branch];
	auto iter = u.ordered_neighbors.lower_bound({ u,v,static_cast<uint>(branch) });
	if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
	return (*(std::prev(iter)));
}


void find_common_neighbor(Vertex neighbor, Vertex root,
	std::vector<Vertex>& share_neighbor, m_Graph& mst) {
	auto adj_pair = boost::adjacent_vertices(neighbor, mst.graph);
	std::set<Vertex> neighbor_neighbor(adj_pair.first, adj_pair.second);
	adj_pair = boost::adjacent_vertices(root, mst.graph);
	std::set<Vertex> vertex_neighbor(adj_pair.first, adj_pair.second);
	//std::sort(vertex_neighbor.begin(), vertex_neighbor.end());
	//std::sort(neighbor_neighbor.begin(), neighbor_neighbor.end());
	std::set_intersection(vertex_neighbor.begin(), vertex_neighbor.end(),
		neighbor_neighbor.begin(), neighbor_neighbor.end(),
		std::back_inserter(share_neighbor));
	return;
}