#pragma once
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/properties.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/graph/connected_components.hpp>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/Fuzzy_sphere.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <CGAL/boost/iterator/counting_iterator.hpp>
#include <CGAL/pca_estimate_normals.h>
#include <unordered_map>
#include <unordered_set>
#include <queue>
#include "math_util.h"
#include "etf.h"

#ifndef GRAPH_FUNC_H
#define GRAPH_FUNC_H

typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point;
typedef Kernel::Vector_3 Vector;
typedef boost::tuple<Point, Vector> PointWithNormal;

struct VertexProperty {
	int id = -1;
	int normal_rep = -1;
	bool colored;
	Point coords;
	Vector normal;
	std::vector<bool> faceExist;
	float distance;
    struct Neighbor{
        float angle;
        uint v;
        mutable uint tree_id = 0;
        mutable bool faceExist = false;

        Neighbor(const VertexProperty& u, const VertexProperty& v, uint id){
            this->v = id;
            this->angle = cal_radians_3d(v.coords - u.coords, u.normal);
        }

        bool operator < (const Neighbor& rhs) const{
            return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
        }

        bool operator<(Neighbor& rhs) const{
            return this->angle < rhs.angle || this->angle == rhs.angle && this->v != rhs.v;
        }
    };

    std::set<Neighbor> ordered_neighbors;
};

struct EdgeProperty {
	float weight;
	int appear_time = 0;
	bool bone_edge;
	int count_weight = 1;
	bool hole_visited = false;
};

typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, VertexProperty, EdgeProperty> Graph;
typedef boost::adjacency_list< boost::vecS, boost::vecS, boost::undirectedS,
	boost::property< boost::vertex_distance_t, float >, boost::property< boost::edge_weight_t, float > >
	s_Graph;
typedef boost::property_map< s_Graph, boost::edge_weight_t >::type s_weightMap;
typedef Graph::vertex_descriptor Vertex;
typedef Graph::vertex_iterator vertex_iter;
typedef Graph::edge_iterator edge_iter;
typedef Graph::edge_descriptor Edge;
typedef VertexProperty::Neighbor Neighbor;

class Face {
public:
	Face() {
		ids[0] = 0;
		ids[1] = 0;
		ids[2] = 0;
	}
	Face(Vertex in_id1, Vertex in_id2, Vertex in_id3) {
		ids[0] = in_id1;
		ids[1] = in_id2;
		ids[2] = in_id3;
	}
	Vertex ids[3];
};

struct m_Graph {
	Graph graph;
	vertex_iter vi_start;
	vertex_iter vi_end;
	float total_edge_length = 0;
	int face_loop_id = 0;
	bool isLeastEdge = false;
	bool isEuclidean = false;
	bool isFinalize = false;
	float meanEdgeScore = 0.0;
	int exp_genus = -1;
    ETF etf;
	//std::vector<Vertex> vertex_map;
	//std::unordered_map<int, std::vector<int>> beacon_tower;
	//int tower_step_size = 500;
	void remove_edge(Vertex source, Vertex target) {
		Edge ed = boost::edge(source, target, this->graph).first;
		this->total_edge_length -= this->graph[ed].weight;
		this->graph.remove_edge(ed);
		return;
	}
	void remove_neighbor(Vertex root, Vertex neighbor) {
		/*
		int idx = 0;
		for (auto point : this->graph[root].ordered_neighbors) {
			if (point == int(neighbor))
				break;
			idx++;
		}
		this->graph[root].ordered_neighbors.erase(this->graph[root].ordered_neighbors.begin() + idx);
		this->graph[root].ordered_radians.erase(this->graph[root].ordered_radians.begin() + idx);
		 */
		auto& u = this->graph[root];
		auto& v = this->graph[neighbor];
		u.ordered_neighbors.erase(Neighbor(u, v, neighbor));
		return;
	}
	void insert_neighbor(Vertex root, Vertex neighbor) {
		const auto& u = this->graph[root];
		const auto& v = this->graph[neighbor];
		this->graph[root].ordered_neighbors.insert(Neighbor(u, v, neighbor));
		return;
	}
	Edge add_edge(int source, int target, float weight = 0.,
		bool bone_edge = false) {

		//if ((source == 388587 && target == 387884) || (target == 388587 && source == 387884))
		//	std::cout << "debug" << std::endl;

		std::pair<Edge, bool> result = boost::add_edge(source, target, this->graph);
		//if (weight * weight > maxEdgeScore)
		//	maxEdgeScore = weight * weight;
		//std::cout << "Edge: " << source << " " << target << std::endl;
		this->graph[result.first].weight = weight;
		this->graph[result.first].bone_edge = bone_edge;
		this->total_edge_length += weight;
		insert_neighbor(source, target);
		insert_neighbor(target, source);
		return result.first;
	}
};

struct my_visitor : boost::default_dijkstra_visitor {
	using base = boost::default_dijkstra_visitor;
	struct found {};
	struct enough {};

	my_visitor(Vertex vd, int in_threshold, int* in_dist) : destination(vd),
		threshold(in_threshold), m_dist_map(in_dist) {}

	void examine_edge(Edge e, Graph const& g)
	{
		if (m_dist_map[boost::target(e, g)] >= threshold&&
			m_dist_map[boost::target(e, g)] <= boost::num_vertices(g)) {
			throw enough{};
		}
		base::examine_edge(e, g);
	}

	void finish_vertex(Vertex v, const Graph& g) {

		if (v == destination)
			throw found{};
		base::finish_vertex(v, g);
	}

	int* m_dist_map;
	Vertex destination;
	int threshold;
};

typedef boost::tuple<Point, int>										Point_and_int;
typedef CGAL::Search_traits_3<Kernel>                                   Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,
	CGAL::Nth_of_tuple_property_map<0, Point_and_int>,
	Traits_base>														Traits;
typedef CGAL::Fuzzy_sphere<Traits>										Sphere;
typedef CGAL::Sliding_midpoint<Traits>									Median;
typedef CGAL::Euclidean_distance<Traits>								Distance;
typedef CGAL::Orthogonal_k_neighbor_search<Traits, Distance, Median>	K_neighbor_search;
typedef K_neighbor_search::Tree											Tree;
typedef CGAL::Parallel_if_available_tag Concurrency_tag;

void kNN_search(int, const Point&, const Tree&, const Distance&, int,
	std::vector<int>&, std::vector<float>&, bool isContain = true);

void NN_search(int, const Point&, const Tree&, const Distance&, float,
	std::vector<int>&, std::vector<float>&, bool isContain = true);

float find_components(std::vector<Point>&,
	std::vector<std::vector<Point>>&, std::vector<Point>&,
	std::vector<std::vector<Point>>&, std::vector<Vector>&,
	std::vector<std::vector<Vector>>&, const Tree&, const Distance&,
	int, bool, float, float);

void init_graph(const std::vector<Point>& vertices, const std::vector<Point>& smoothed_v,
                const std::vector<Vector>& normals, const Tree& kdTree, const Distance& tr_dist,
                int k, s_Graph& dist_graph, s_weightMap& weightmap, bool isEuclidean, std::vector<float>& max_length,
				int exp_genus, std::vector<float>& pre_max_length, float cross_conn_thresh);

void print_path(std::vector<Vertex>&, Vertex, std::vector<Vertex>&);

int find_shortest_path(const m_Graph&, Vertex, Vertex, int, std::vector<Vertex>&);

void weighted_smooth(const std::vector<Point>& vertices,
	std::vector<Point>& smoothed_v, const std::vector<Vector>& normals,
	const Tree& kdTree, const Distance& tr_dist,
	float diagonal_length);

#endif