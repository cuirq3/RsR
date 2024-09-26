#pragma once
#include "graph_utils.h"

#ifndef RS_H
#define RS_H

const Neighbor& get_neighbor_info(const m_Graph& g, const Vertex& root, const Vertex& branch);

const Neighbor& predecessor(const m_Graph& g, const Vertex& root, const Vertex& branch);

const Neighbor& successor(const m_Graph& g, const Vertex& root, const Vertex& branch);

void find_common_neighbor(Vertex, Vertex,
	std::vector<Vertex>&, m_Graph&);

#endif //RS_H