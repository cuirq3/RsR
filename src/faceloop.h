#pragma once
#include "graph_utils.h"

#ifndef FACELOOP_H
#define FACELOOP_H

extern long t_loop;
extern int bettiNum_1;

std::string generate_path_key(Vertex, Vertex);

// Not used
int generate_face_loop(m_Graph&, Vertex, Vertex,
	std::vector<Vertex>&,
	std::unordered_set<int>&, bool isClockWise = true);

void init_face_loop_label(m_Graph&);

// Not used
void sync_face_loop(m_Graph&);

int update_face_loop(m_Graph& g, Vertex start_v,
	Vertex end_v, bool isTower = false);

void maintain_vertex_face_loop(m_Graph& g, Vertex start_v, Vertex end_v);

uint get_face_loop_id(m_Graph& g, Vertex v, int id);

int get_face_loop_id(m_Graph& g, int& id, Vertex current_v, Vertex last_v);

void maintain_face_loop(m_Graph& g, const Edge& e);

#endif //FACELOOP_H