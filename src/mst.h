#pragma once
//#include "triangulation.h"
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include "graph_utils.h"
#include "RS.h"
#include "timer.h"
#include "faceloop.h"
#include "normal.h"
#include "io.h"

#ifndef MST_H
#define MST_H

typedef std::pair<Vertex, Vertex> m_Edge;
typedef std::pair<float, int> m_Edge_length;
typedef std::pair<int, std::string> m_face_pair;

inline bool edge_comparator(const m_Edge_length& l, const m_Edge_length& r) {
	return l.first < r.first;
}
inline bool face_comparator(const m_face_pair& l, const m_face_pair& r) {
	return l.first > r.first;
}

bool isIntersecting(m_Graph& mst, Vertex v1, Vertex v2, Vertex v3, Vertex v4);

void build_mst(m_Graph&, std::vector<Vertex>&,
	bool, const std::vector<Point>&,
	std::vector<Vector>&);

bool geometry_check(m_Graph& mst, m_Edge& candidate, Tree& kdTree, Distance& tr_dist);

bool Vanilla_check(m_Graph& g, m_Edge& candidate, Tree& kdTree, Distance& tr_dist);

void connect_handle(const std::vector<Point>&, Tree&, Distance&,
	m_Graph&, std::vector<Vertex>& connected_handle_root, std::vector<int>& betti,
	int k = 8, bool isEuclidean = false, int step_thresh = 10);

#endif //MST_H