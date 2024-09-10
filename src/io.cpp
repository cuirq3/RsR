#include "io.h"

void io_system::read_obj(fs::path file_path, std::vector<Point>& vertices,
	std::vector<Vector>& normals,
	std::vector<Face>& faces) {
	fs::ifstream file(file_path);
	std::string line;
	while (std::getline(file, line))
	{
		std::vector<std::string> info;
		int pos = 0;
		while ((pos = line.find(" ")) != std::string::npos) {
			info.push_back(line.substr(0, pos));
			line.erase(0, pos + 1);
		}
		info.push_back(line);
		if (info.size() == 0) {
			continue;
		}
		if (info.at(0) == "v") {
			Point vertex(std::stof(info.at(1)),
				std::stof(info.at(2)), std::stof(info.at(3)));
			vertices.push_back(vertex);
		}
		if (info.at(0) == "vn") {
			Vector normal(std::stof(info.at(1)),
				std::stof(info.at(2)), std::stof(info.at(3)));
			normals.push_back(normal);
		}
		
		if (info.at(0) == "f") {
			Face face(std::stoi(info.at(1)) - 1,
				std::stoi(info.at(2)) - 1, std::stoi(info.at(3)) - 1);
			faces.push_back(face);
		}
	}
	file.close();
	return;
}

//// CGAL can not read vertex normal from .obj file now
//void io_system::read_obj_CGAL(std::string file_path, std::vector<Point>& vertices,
//	std::vector<std::vector<std::size_t>>& faces_ref) {
//	CGAL::IO::read_OBJ(file_path, vertices, faces_ref
//		, CGAL::parameters::point_map(CGAL::First_of_pair_property_map<PointVectorPair>())
//		.normal_map(CGAL::Second_of_pair_property_map<PointVectorPair>()));
//	return;
//}

void io_system::export_obj(m_Graph& g, fs::path out_path, std::vector<Face>& faces) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < boost::num_vertices(g.graph); i++) {
		Vertex this_v = i;
		Point this_coords = g.graph[this_v].coords;
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}
	// Write vertex normal
	file << std::endl;
	file << "# List of vertex normals" << std::endl;
	for (int i = 0; i < boost::num_vertices(g.graph); i++) {
		Vertex this_v = i;
		Vector this_normal = g.graph[this_v].normal;
		file << "vn " << std::to_string(this_normal.x())
			<< " " << std::to_string(this_normal.y())
			<< " " << std::to_string(this_normal.z()) << std::endl;
	}

	// Write faces
	file << std::endl;
	file << "# Polygonal face element" << std::endl;
	for (auto face : faces) {
		file << "f " << std::to_string(face.ids[0] + 1)
			<< " " << std::to_string(face.ids[1] + 1)
			<< " " << std::to_string(face.ids[2] + 1) << std::endl;
	}
	file.close();
	return;
}

void io_system::export_obj(std::vector<Point>& in_vertices, m_Graph& g, fs::path out_path,
	std::vector<Face>& faces) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < in_vertices.size(); i++) {
		Point this_coords = in_vertices[i];
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}
	// Write vertex normal
	file << std::endl;
	file << "# List of vertex normals" << std::endl;
	for (int i = 0; i < boost::num_vertices(g.graph); i++) {
		Vertex this_v = i;
		Vector this_normal = g.graph[this_v].normal;
		file << "vn " << std::to_string(this_normal.x())
			<< " " << std::to_string(this_normal.y())
			<< " " << std::to_string(this_normal.z()) << std::endl;
	}

	// Write faces
	file << std::endl;
	file << "# Polygonal face element" << std::endl;
	for (auto face : faces) {
		file << "f " << std::to_string(face.ids[0] + 1)
			<< " " << std::to_string(face.ids[1] + 1)
			<< " " << std::to_string(face.ids[2] + 1) << std::endl;
	}
	file.close();
	return;
}

void io_system::export_obj(m_Model& g, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < g.vertices.size(); i++) {
		file << "v " << std::to_string(g.vertices.at(i).x())
			<< " " << std::to_string(g.vertices.at(i).y())
			<< " " << std::to_string(g.vertices.at(i).z()) << std::endl;
	}
	// Write vertex normal
	file << std::endl;
	file << "# List of vertex normals" << std::endl;
	for (int i = 0; i < g.normals.size(); i++) {
		Vector this_normal = g.normals[i];
		file << "vn " << std::to_string(this_normal.x())
			<< " " << std::to_string(this_normal.y())
			<< " " << std::to_string(this_normal.z()) << std::endl;
	}

	// Write faces
	file << std::endl;
	file << "# Polygonal face element" << std::endl;
	for (auto face : g.faces) {
		file << "f " << std::to_string(face.ids[0] + 1)
			<< " " << std::to_string(face.ids[1] + 1)
			<< " " << std::to_string(face.ids[2] + 1) << std::endl;
	}
	file.close();
	return;
}

void io_system::export_obj(std::vector<Point>& in_vertices, m_Model& g, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < in_vertices.size(); i++) {
		file << "v " << std::to_string(in_vertices.at(i).x())
			<< " " << std::to_string(in_vertices.at(i).y())
			<< " " << std::to_string(in_vertices.at(i).z()) << std::endl;
	}
	// Write vertex normal
	file << std::endl;
	file << "# List of vertex normals" << std::endl;
	for (int i = 0; i < g.normals.size(); i++) {
		Vector this_normal = g.normals[i];
		file << "vn " << std::to_string(this_normal.x())
			<< " " << std::to_string(this_normal.y())
			<< " " << std::to_string(this_normal.z()) << std::endl;
	}

	// Write faces
	file << std::endl;
	file << "# Polygonal face element" << std::endl;
	for (auto face : g.faces) {
		file << "f " << std::to_string(face.ids[0] + 1)
			<< " " << std::to_string(face.ids[1] + 1)
			<< " " << std::to_string(face.ids[2] + 1) << std::endl;
	}
	file.close();
	return;
}

void io_system::export_obj(std::vector<Point>& vertices, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < vertices.size(); i++) {
		file << "v " << std::to_string(vertices.at(i).x())
			<< " " << std::to_string(vertices.at(i).y())
			<< " " << std::to_string(vertices.at(i).z()) << std::endl;
	}
	file.close();
	return;
}

void io_system::export_obj(std::vector<Point>& vertices,
	std::vector<Vector>& normals, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < vertices.size(); i++) {
		file << "v " << std::to_string(vertices.at(i).x())
			<< " " << std::to_string(vertices.at(i).y())
			<< " " << std::to_string(vertices.at(i).z()) << std::endl;
	}

	// Write vertex normal
	file << std::endl;
	file << "# List of vertex normals" << std::endl;
	for (int i = 0; i < normals.size(); i++) {
		file << "vn " << std::to_string(normals.at(i).x())
			<< " " << std::to_string(normals.at(i).y())
			<< " " << std::to_string(normals.at(i).z()) << std::endl;
	}

	file.close();
	return;
}

void io_system::export_graph(m_Graph& g, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < boost::num_vertices(g.graph); i++){
		Vertex this_v = i;
		Point this_coords = g.graph[this_v].coords;
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}

	// Write lines
	file << std::endl;
	file << "# Line element" << std::endl;
	//std::cout << boost::num_edges(g.graph) << std::endl;
	Graph::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(g.graph); ei != ei_end; ei++)
		file << "l " << std::to_string(boost::source(*ei, g.graph) + 1)
		<< " " << std::to_string(boost::target(*ei, g.graph) + 1) << std::endl;
	file.close();
	return;
}

void io_system::export_graph(m_Graph& g, fs::path out_path, std::vector<Point>& vertices) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < vertices.size(); i++) {
		Point this_coords = vertices[i];
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}

	// Write lines
	file << std::endl;
	file << "# Line element" << std::endl;
	//std::cout << boost::num_edges(g.graph) << std::endl;
	Graph::edge_iterator ei, ei_end;
	for (boost::tie(ei, ei_end) = boost::edges(g.graph); ei != ei_end; ei++)
		file << "l " << std::to_string(boost::source(*ei, g.graph) + 1)
		<< " " << std::to_string(boost::target(*ei, g.graph) + 1) << std::endl;
	file.close();
	return;
}

void io_system::export_edges(m_Graph& g, std::vector<Vertex>& roots, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < roots.size(); i++) {
		Vertex this_v = roots[i];
		Point this_coords = g.graph[this_v].coords;
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}

	// Write lines
	file << std::endl;
	file << "# Line element" << std::endl;
	for (int i = 0; i < roots.size(); i += 2)
		file << "l " << i + 1
		<< " " << i + 2 << std::endl;
	file.close();
	return;
}

void export_continuous_edges(m_Graph& g, std::vector<Vertex>& roots, fs::path out_path) {
	fs::ofstream file(out_path);
	// Write vertices
	file << "# List of geometric vertices" << std::endl;
	for (int i = 0; i < roots.size(); i++) {
		Vertex this_v = roots[i];
		Point this_coords = g.graph[this_v].coords;
		file << "v " << std::to_string(this_coords.x())
			<< " " << std::to_string(this_coords.y())
			<< " " << std::to_string(this_coords.z()) << std::endl;
	}

	// Write lines
	file << std::endl;
	file << "# Line element" << std::endl;
	for (int i = 0; i < roots.size()-1; i++)
		file << "l " << i + 1
		<< " " << i + 2 << std::endl;
	file.close();
	return;
}


void io_system::export_betti(std::vector<int>& betti, fs::path out_path) {
	fs::ofstream file(out_path);
	for (int i = 0; i < betti.size(); i++) {
		file << betti[i] << std::endl;
	}
	return;
}

/*
void io_system::read_ply(fs::path file_path, std::vector<Point>& vertices,
	std::vector<Vector>& normals, std::vector<Face>& GT_faces,
	std::vector<Point>& GT_vertices) {
	open3d::geometry::TriangleMesh mesh;
	open3d::io::ReadTriangleMeshOptions option;
	if (!open3d::io::ReadTriangleMeshFromPLY(file_path.string(), mesh, option)) {
		std::cout << "Error reading ply file." << std::endl;
	}

	for (auto triangle : mesh.triangles_) {
		GT_faces.push_back(Face(triangle.x(), triangle.y(), triangle.z()));
	}

	for (auto point : mesh.vertices_) {
		GT_vertices.push_back(Point(float(point.x()), float(point.y()), float(point.z())));
	}

	int num_sample = 1000000;
	open3d::geometry::PointCloud sample = *mesh.SamplePointsUniformly(num_sample);
	int idx = 0;
	for (auto point : sample.points_) {
		vertices.push_back(Point(float(point.x()), float(point.y()), float(point.z())));
		if (sample.normals_.size() == num_sample) {
			auto this_normal = sample.normals_[idx];
			normals.push_back(Vector(float(this_normal.x()), float(this_normal.y()), float(this_normal.z())));
		}
	}
	return;
}
*/
void io_system::read_pc_ply(fs::path file_path, std::vector<Point>& vertices,
	std::vector<Vector>& normals, std::vector<Point>& GT_vertices) {
	happly::PLYData pc(file_path.string());

	auto vertex_pos = pc.getVertexPositions();
	for (auto point : vertex_pos) {
		if (std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]))
			GT_vertices.push_back(Point(float(point[0]), float(point[1]), float(point[2])));
	}

	int num_sample = 2000000;
	std::vector<int> idx_array(vertex_pos.size());

	if (vertex_pos.size() > 2000000) {
		// Populate the vector
		for (int i = 0; i < vertex_pos.size(); i++) {
			idx_array[i] = i;
		}

		std::random_device rd;
		std::mt19937 g(rd());
		std::shuffle(idx_array.begin(), idx_array.end(), g);
		idx_array.resize(num_sample);
	}

	//std::vector<double> nx = pc.getElement("vertex").getProperty<double>("nx");
	//std::vector<double> ny = pc.getElement("vertex").getProperty<double>("ny");
	//std::vector<double> nz = pc.getElement("vertex").getProperty<double>("nz");

	for (auto idx : idx_array) {
		auto point = vertex_pos[idx];
		if (std::isfinite(point[0]) && std::isfinite(point[1]) && std::isfinite(point[2]))
			vertices.push_back(Point(float(point[0]), float(point[1]), float(point[2])));
		/*if (nx.size() == vertex_pos.size()) {
			if (std::isfinite(nx[idx]) && std::isfinite(ny[idx]) && std::isfinite(nz[idx]))
				normals.push_back(Vector(float(nx[idx]), float(ny[idx]), float(nz[idx])));
		}*/
	}
	return;
}