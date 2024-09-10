#include "faceloop.h"
#include "RS.h"

long t_loop = 0;
int bettiNum_1 = 0;
std::string generate_path_key(Vertex v1, Vertex v2) {
	return std::to_string(int(v1)) + "+" + std::to_string(int(v2));
}

void maintain_face_loop(m_Graph& g, const Edge& e){
    Vertex source = e.m_source;
    Vertex target = e.m_target;
    auto this_v_tree = predecessor(g,source, target).tree_id;
    auto neighbor_tree = predecessor(g, target, source).tree_id;

    auto result = g.etf.insert(this_v_tree,neighbor_tree);
	auto u = result.first;
	auto v = result.second;
    get_neighbor_info(g, source, target).tree_id = u;
    get_neighbor_info(g, target, source).tree_id = v;

	bettiNum_1++;
	return;
}

void init_face_loop_label(m_Graph& g) {
	Vertex start_v = 0;
	Vertex last_vertex = start_v;
	int loop_step = 0;
	Vertex current_vertex = g.graph[start_v].ordered_neighbors.begin()->v;
	//g.beacon_tower[0] = std::vector<int>();
	//std::unordered_set<int> faceloop_member;
	std::vector<int> towers;
	do {
        //current_vertex = travel_cw(g, current_vertex, last_vertex);
		/*std::vector<int> current_neighbors = g.graph[current_vertex].ordered_neighbors;
		int idx = -1;
		for (int i = 0; i < current_neighbors.size(); i++) {
			if (current_neighbors[i] == int(last_vertex)) {
				idx = i;
				break;
			}
		}
		if (idx == -1)
			throw std::runtime_error("No neighbor found (caused by wrong ordered_neighbors)");
		int this_v_id = g.graph[temp].
			ordered_neighbors[python_mod((idx - 1), current_neighbors.size())];
		current_vertex = this_v_id;*/

		//// Initialize set and map
		//g.graph[temp].face_loop.insert(0);
		//if (g.graph[temp].face_loop_time.count(0) == 0)
		//	g.graph[temp].face_loop_time[0] = 1;
		//else
		//	g.graph[temp].face_loop_time[0]++;

		//g.graph[temp].path_last.push_back(int(last_vertex));

        auto& next_neighbor = predecessor(g, current_vertex, last_vertex);

        next_neighbor.tree_id = g.etf.accumulate();

		last_vertex = current_vertex;
        current_vertex = next_neighbor.v;

		loop_step++;

		/*faceloop_member.insert(this_v_id);
		if (loop_step % g.tower_step_size == 0) {
			g.graph[current_vertex].faceloop_member[0] = faceloop_member;
			towers.push_back(this_v_id);
			faceloop_member.clear();
		}*/
		//std::cout << current_vertex << std::endl;
	} while (current_vertex != g.graph[start_v].ordered_neighbors.begin()->v || last_vertex != start_v);

	/*if (!faceloop_member.empty()) {
		g.graph[last_vertex].faceloop_member[0] = faceloop_member;
		towers.push_back(0);
	}
	g.beacon_tower[0] = towers;*/

	std::cout << "Loop step initialization finished after " + std::to_string(loop_step) + " steps." << std::endl;
	return;
}
