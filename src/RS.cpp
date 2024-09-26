#include"RS.h"

typedef std::pair<float, Vertex> mypair;
bool comparator(const mypair& l, const mypair& r) {
	return l.first < r.first;
}

/**
 * @brief Get the neighbor information
 *
 * @param g: current graph
 * @param root: root vertex index
 * @param branch: the outgoing branch
 * 
 * @return reference to the neighbor struct
 */
const Neighbor& get_neighbor_info(const m_Graph& g, const Vertex& root, const Vertex& branch){
    auto& u = g.graph[root];
    auto& v = g.graph[branch];
    auto iter = u.ordered_neighbors.lower_bound({u,v,static_cast<uint>(branch)});
    return (*iter);
}

/**
 * @brief Get the next neighbor information
 *
 * @param g: current graph
 * @param root: root vertex index
 * @param branch: current outgoing branch
 *
 * @return reference to the next neighbor struct
 */
const Neighbor& successor(const m_Graph &g, const Vertex &root, const Vertex &branch) {
    auto& u = g.graph[root];
    auto& v = g.graph[branch];
    auto iter = u.ordered_neighbors.upper_bound(Neighbor(u,v,branch));
    if(iter == u.ordered_neighbors.end()) iter = u.ordered_neighbors.begin(); // Wrap around
    return (*iter); // This is honestly not good practice - ONLY modification of tree_id
}

/**
 * @brief Get last neighbor information
 *
 * @param g: current graph
 * @param root: root vertex index
 * @param branch: current outgoing branch
 *
 * @return reference to last neighbor struct
 */
const Neighbor& predecessor(const m_Graph& g, const Vertex& root, const Vertex& branch) {
	auto& u = g.graph[root];
	auto& v = g.graph[branch];
	auto iter = u.ordered_neighbors.lower_bound({ u,v,static_cast<uint>(branch) });
	if (iter == u.ordered_neighbors.begin()) iter = u.ordered_neighbors.end(); // Wrap around
	return (*(std::prev(iter)));
}

/**
 * @brief Find the common neighbors that two vertices are sharing
 *
 * @param neighbor: one vertex
 * @param root: the other vertex
 * @param share_neighbor: [OUT] common neighbors these two vertices share
 * @param mst: tcurrent graph
 *
 * @return reference to last neighbor struct
 */
void find_common_neighbor(Vertex neighbor, Vertex root,
	std::vector<Vertex>& share_neighbor, m_Graph& g) {
	auto adj_pair = boost::adjacent_vertices(neighbor, g.graph);
	std::set<Vertex> neighbor_neighbor(adj_pair.first, adj_pair.second);
	adj_pair = boost::adjacent_vertices(root, g.graph);
	std::set<Vertex> vertex_neighbor(adj_pair.first, adj_pair.second);
	std::set_intersection(vertex_neighbor.begin(), vertex_neighbor.end(),
		neighbor_neighbor.begin(), neighbor_neighbor.end(),
		std::back_inserter(share_neighbor));
	return;
}