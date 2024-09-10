#pragma once
#include "RS.h"
#include "faceloop.h"
#include "timer.h"
#include "graph_utils.h"

#ifndef TRIANGULATION_H
#define TRIANGULATION_H

struct m_cmp {
    bool operator()(const std::pair<std::vector<Vertex>, float> &left,
               const std::pair<std::vector<Vertex>, float> &right) {
        return (left.second) > (right.second);
    }
};

typedef std::priority_queue<std::pair<std::vector<Vertex>, float>,
	std::vector<std::pair<std::vector<Vertex>, float>>,
	m_cmp> m_priority_queue;

bool explore(m_Graph& G, int i, m_priority_queue& queue,
	std::unordered_set<std::string>& faces_in_queue, float avg_edge_length,
	std::vector<float>& length_thresh);

void checkAndForce(Vertex v_u, Vertex v_w, m_Graph& G, m_priority_queue& queue,
	std::vector<float>& length_thresh);

bool check_face_overlap(m_Graph&, std::vector<Vertex>&,
	const Tree&, const Distance&, int k = 10);

bool check_branch_validity(m_Graph&, Vertex, Vertex, Vertex);

bool check_validity(m_Graph&, std::pair<std::vector<Vertex>, float>&,
	const Tree&, const Distance&, bool, bool);

void add_face(m_Graph&, std::vector<Vertex>&,
	std::vector<Face>&);

void triangulate(std::vector<Face>&, m_Graph&,
	const Tree&, const Distance&, bool,
	bool isEuclidean, std::vector<float>& length_thresh,
	std::vector<Vertex>& connected_handle_root, std::vector<int>& betti, bool isFinalize=false);

bool routine_check(m_Graph& mst, std::vector<Vertex>& triangle);

bool register_face(m_Graph& mst, Vertex v1, Vertex v2, std::vector<Face>& faces,
	Tree& KDTree, Distance& tr_dist, float);

#endif //TRIANGULATION_H