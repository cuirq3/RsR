#pragma once
#include <iostream>
#include <random>
#include <algorithm>
#include "boost/filesystem/operations.hpp" // includes boost/filesystem/path.hpp
#include "boost/filesystem/fstream.hpp"    // ditto                       // for std::cout
#include <CGAL/IO/OBJ.h>
#include "graph_utils.h"
#include "eval.h"
#include "happly.h"

#ifndef IO_H
#define IO_H

namespace fs = boost::filesystem;

class io_system {
public:
	io_system() {
	}
	void read_obj(fs::path, std::vector<Point>&, std::vector<Vector>&, std::vector<Face>&);
	void export_obj(m_Graph&, fs::path, std::vector<Face>&);
	void export_obj(std::vector<Point>& in_vertices, m_Graph& g, fs::path out_path,
		std::vector<Face>& faces);
	void export_obj(m_Model&, fs::path);
	void export_obj(std::vector<Point>& in_vertices, m_Model& g, fs::path out_path);
	void export_obj(std::vector<Point>&, fs::path);
	void export_obj(std::vector<Point>&, std::vector<Vector>&, fs::path);
	void export_graph(m_Graph&, fs::path);
	void export_graph(m_Graph& g, fs::path out_path, std::vector<Point>& vertices);
	void export_edges(m_Graph& g, std::vector<Vertex>& roots, fs::path out_path);
	void export_betti(std::vector<int>&, fs::path out_path);
	void read_ply(fs::path, std::vector<Point>&, std::vector<Vector>&,
		std::vector<Face>&, std::vector<Point>&);
	void read_pc_ply(fs::path, std::vector<Point>&,
		std::vector<Vector>&, std::vector<Point>&);
};

void export_continuous_edges(m_Graph& g, std::vector<Vertex>& roots, fs::path out_path);

#endif //IO_H