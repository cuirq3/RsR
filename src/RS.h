#pragma once
#include "graph_utils.h"

#ifndef RS_H
#define RS_H

void init_rotation_system(m_Graph&);

void maintain_rotation_system(m_Graph&, Edge&);

void find_radians_id(const m_Graph&, const Vertex&,
	const Vertex&, int&, int&, float&);

void find_radians_id_no_clamp(const m_Graph& G, const Vertex& root,
                              const Vertex& branch, int& next_id, int& last_id, float& this_radian);

const Neighbor& get_neighbor_info(const m_Graph& g, const Vertex& root, const Vertex& branch);

const Neighbor& predecessor(const m_Graph& g, const Vertex& root, const Vertex& branch);

const Neighbor& successor(const m_Graph& g, const Vertex& root, const Vertex& branch);

void find_common_neighbor(Vertex, Vertex,
	std::vector<Vertex>&, m_Graph&);

#endif //RS_H