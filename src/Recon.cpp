#include "Recon.h"

/**
 * @brief Read the input config file and initialize
 *
 * @param config_path Path to the config file
 * @return None
 */
void Reconstructor::read_config(fs::path config_path) {
    fs::ifstream config(config_path);
    string line;
    std::vector<string> tokens;
    while (getline(config, line)) {
        std::size_t pos = 0, found;
        while (true) {
            found = line.find(' ', pos);
            if (found != std::string::npos) {
                std::string token = line.substr(pos, found - pos);
                tokens.push_back(token);
                pos = found + 1;
            }
            else {
                found = line.find('\n', pos);
                std::string token = line.substr(pos, found - pos);
                tokens.push_back(token);
                break;
            }
        }
        string instruct = tokens[0];
        tokens[1].erase(tokens[1].find_last_not_of("\n") + 1);
        if (instruct == "root") {
            root_path = fs::path(tokens[1]);
        }
        if (instruct == "model_name") {
            if (tokens[1] == "all") {
                model_name = tokens[1];
                model_path = root_path;
            }
            else {
                model_name = tokens[1];
                model_path = root_path / tokens[1];
            }
        }
        if (instruct == "out_root") {
            out_root_path = fs::path(tokens[1]);
        }
        if (instruct == "out_name") {
            out_name = tokens[1];
        }
        if (instruct == "isEuclidean") {
            isEuclidean = (tokens[1] == "true");
        }
        if (instruct == "isGTNormal") {
            isGTNormal = (tokens[1] == "true");
        }
        if (instruct == "isDebug") {
            isDebug = (tokens[1] == "true");
        }
        if (instruct == "isLeastEdgeScore") {
            isLeastEdgeScore = (tokens[1] == "true");
        }
        if (instruct == "isProgressive") {
            isProgressive = (tokens[1] == "true");
        }
        if (instruct == "isNoiseExp") {
            isNoiseExperiment = (tokens[1] == "true");
        }
        if (instruct == "k") {
            k = (std::stoi(tokens[1]));
        }
        if (instruct == "genus") {
            exp_genus = (std::stoi(tokens[1]));
        }
        tokens.clear();
    }
}

/**
 * @brief Randomly sample points within a [0,1]^3 box, for experiments
 *
 * @param n Number of samples
 * @return out_path Path to the output file
 */
fs::path Reconstructor::generate_box_samples(int n) {
    io_system IO;
    std::vector<Point> vertices;
    srand((unsigned)time(NULL));

    for (int i = 0; i < n; i++)
    {
        Point vertex((float)rand() / RAND_MAX, (float)rand() / RAND_MAX, (float)rand() / RAND_MAX);
        vertices.push_back(vertex);
    }
    model_name= "box_org.obj";
    fs::path out_path = out_root_path / model_name;
    IO.export_obj(vertices, out_path);
    return out_path;
}

/**
 * @brief Read the input config file and initialize
 *
 * @param config_path Path to the config file
 * @return None
 */
void Reconstructor::noise_experiment() {
    fs::path file_path = model_path / model_name;
    isGTNormal = false;
    std::string experiment_log;
    std::string path_string = file_path.string();
    int start = path_string.find_last_of('\\');
    int str_len = path_string.find_last_of('.') - start;
    std::string this_model_name = path_string.substr(start + 1, str_len - 1);
    out_name = this_model_name;
    fs::path recon_root = out_root_path / "Recon" / noise_type;
    if (!fs::exists(recon_root))
        fs::create_directories(recon_root);
    /*fs::path gt_root = out_root_path / "GT" / noise_type;
    if (!fs::exists(gt_root))
        fs::create_directories(gt_root);*/

    // Load Original GT model
    io_system IO;
    m_Model GT_model_org;
    {
        vector<Point> vertices;
        vector<Vector> normals;
        vector<Face> GT_faces;
        IO.read_obj(file_path, vertices, normals, GT_faces);
        GT_model_org.vertices = vertices;
        GT_model_org.faces = GT_faces;
    }

    // Without noise
    isEuclidean = false;
    //reconstruct_single(Recon_model, GT_model, noise_type, 0.0, 0.0, false);
    //float CD = calculate_chamfer(GT_model_org, Recon_model, chamfer_sample);
    //experiment_log += std::to_string(CD) + "\n";

    for (auto sigma : sigmas) {
        for (auto amplitude : amplitudes) {
            m_Model GT_model;
            m_Model Recon_model;
            isFaceLoop = true;
            std::cout << "Experimenting sigma " + std::to_string(sigma) + " amplitude " + std::to_string(amplitude) + "..." << std::endl;
            std::string this_out_name = this_model_name + std::to_string(sigma) + "_" + std::to_string(amplitude);
            out_root_path = recon_root;
            out_name = this_out_name;
            reconstruct_single(Recon_model, GT_model, noise_type, sigma, amplitude, false);
            /*float CD = calculate_chamfer(GT_model_org, Recon_model, chamfer_sample);
            experiment_log += std::to_string(CD) + " ";*/

            // Output GT models
            //IO.export_obj(Recon_model, recon_root / (this_out_name + ".obj"));
        }
        experiment_log += "\n";
    }
    //std::cout << experiment_log << std::endl;
}

void Reconstructor::reconstruct_single(std::string noise_type, float sigma, float amplitude, bool isStepSave, float thresh_r) {
    m_Model ReconModel = m_Model();
    m_Model GTModel = m_Model();
    reconstruct_single(ReconModel,GTModel, noise_type, sigma, amplitude, isStepSave, thresh_r);
}

void Reconstructor::reconstruct_single(m_Model& ReconModel, m_Model& GTModel,
    std::string noise_type, float sigma, float amplitude, bool isStepSave, float thresh_r) {
    
    isDebug = false;
    recon_timer.start("Initialization");
    // Init
    io_system IO;
    vector<Point> GT_vertices;
    vector<Face> GT_faces;
    std::vector<Point> in_vertices;
    std::vector<Vector> in_normals;
    //vector<float> score_percentages{ .5,.6,.7,.8,.9, 1. };
    //vector<float> score_percentages{ 1.9, 2.0 };

    // Read input
    {
        recon_timer.start("Import obj");
        fs::path file_path = root_path / model_name;
        
        std::string file_end = file_path.extension().string();
        if (file_end == ".obj") {
            IO.read_obj(file_path, in_vertices, in_normals, GT_faces);
            GT_vertices = in_vertices;
        }
        else if (file_end == ".ply"){
            IO.read_pc_ply(file_path, in_vertices, in_normals, GT_vertices);
        }
        else {
            std::cout << "Error file type" << std::endl;
            return;
        }

        if (k >= in_vertices.size())
            k = in_vertices.size() - 1;

        recon_timer.end("Import obj");
    }

    GTModel.vertices = GT_vertices;
    GTModel.faces = GT_faces;

    // Estimate normals & orientation & weighted smoothing
    recon_timer.start("Estimate normals");
    vector<Point> in_smoothed_v;
    {
        std::vector<int> indices(in_vertices.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(in_vertices.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(in_vertices.end(), indices.end())));
        Distance tr_dist;
        float diagonal_length;

        if (isGTNormal && in_normals.size() == 0) {
            std::cout << "No normal can be used!" << std::endl;
            isGTNormal = false;
        }
        std::vector<int> zero_normal_id;
        estimate_normal(in_vertices, kdTree, tr_dist, isGTNormal, in_normals,
            zero_normal_id, diagonal_length);
        // Fix zero normal
        if (zero_normal_id.size() > 0) {
            replace_zero_normal(zero_normal_id, kdTree, tr_dist, in_vertices, in_normals);
        }

        // Add position noise
        {
            if (amplitude != 0.) {
                float avg_edge_length = 0.0;
                for (auto face : GT_faces) {
                    for (int i = 0; i < 3; i++) {
                        Point vertex1 = in_vertices.at(i);
                        Point vertex2 = in_vertices.at((i + 1) % 3);
                        float length = norm(vertex1 - vertex2);
                        avg_edge_length += length;
                    }
                }
                avg_edge_length /= GT_faces.size() * 3;
                add_noise(noise_type, in_vertices, in_normals, sigma, amplitude, avg_edge_length);
            }
        }
        if (true) {
            std::cout << "Start first round smoothing ..." << std::endl;
            if (!isEuclidean)
                weighted_smooth(in_vertices, in_smoothed_v, in_normals, kdTree, tr_dist, diagonal_length);
            else
                in_smoothed_v = in_vertices;

            Tree kdTree1(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
                boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
            Distance tr_dist1;

            estimate_normal(in_smoothed_v, kdTree1, tr_dist1, isGTNormal, in_normals,
                zero_normal_id, diagonal_length);
            if (isDebug)
                IO.export_obj(in_smoothed_v, in_normals, out_root_path / (out_name + "_smoothed_v.obj"));

            // Another round of smoothing
            if (true) {
                if (!isEuclidean) {
                    std::cout << "Start second round smoothing ..." << std::endl;
                    std::vector<Point> temp(in_smoothed_v.begin(), in_smoothed_v.end());
                    in_smoothed_v.clear();
                    weighted_smooth(temp, in_smoothed_v, in_normals, kdTree1, tr_dist1, diagonal_length);

                    Tree kdTree2(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
                        boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
                    Distance tr_dist2;

                    estimate_normal(in_smoothed_v, kdTree2, tr_dist2, isGTNormal, in_normals,
                        zero_normal_id, diagonal_length);
                    if (isDebug)
                        IO.export_obj(in_smoothed_v, in_normals, out_root_path / (out_name + "_smoothed_v2.obj"));
                }
            }
        }
        else {
            in_smoothed_v = in_vertices;
        }
        
        // Add normal noise
        if(false){
            float angle = 25. / 180. * CGAL_PI;
            add_normal_noise(angle, in_normals);
        }
    }
    recon_timer.end("Estimate normals");

    recon_timer.start("algorithm");
    std::cout << "find components" << std::endl;
    // Find components
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<vector<Vector>> component_normals;
    {
        recon_timer.start("Kd and DM");

        // Build kdTree - CGAL
        vector<int> indices(in_smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(in_smoothed_v.end(), indices.end())));
        Distance tr_dist;

        // Correct orientation
        if (!isGTNormal) {
            s_Graph g_angle(in_vertices.size());
            s_weightMap weightmap_a = boost::get(boost::edge_weight, g_angle);
            // Init angle based graph
            for (int i = 0; i < in_vertices.size(); i++) {
                Point vertex = in_vertices[i];
                Vector this_normal = in_normals[i];

                std::vector<int> neighbors;
                std::vector<float> dists;
                kNN_search(i, vertex, kdTree, tr_dist, k, neighbors, dists);
                for (int j = 0; j < neighbors.size(); j++) {
                    if (boost::edge(i, neighbors[j], g_angle).second)
                        continue;
                    Vector neighbor_normal = in_normals[neighbors[j]];
                    boost::graph_traits< s_Graph >::edge_descriptor e;
                    bool inserted;
                    boost::tie(e, inserted) = boost::add_edge(i, neighbors[j], g_angle);
                    float angle_weight = cal_angle_based_weight(this_normal, neighbor_normal);
                    if (i == neighbors[j]) {
                        std::cout << "error" << std::endl;
                    }
                    if (angle_weight < 0)
                        std::cout << "error" << std::endl;
                    weightmap_a[e] = angle_weight;
                }
            }
            std::vector<int> component_id(in_vertices.size());
            int num = boost::connected_components(g_angle, &component_id[0]);
            //std::cout << num << std::endl;

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance_a
                = boost::get(boost::vertex_distance, g_angle);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap_a
                = boost::get(boost::vertex_index, g_angle);
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor>
                p_angle(boost::num_vertices(g_angle));
            boost::prim_minimum_spanning_tree
            (g_angle, *boost::vertices(g_angle).first, &p_angle[0],
                distance_a, weightmap_a, indexmap_a,
                boost::default_dijkstra_visitor());
            g_angle.clear();

            g_angle = s_Graph(in_vertices.size());
            for (int i = 0; i < p_angle.size(); i++) {
                if (p_angle[i] == i) {
                    continue;
                }
                else {
                    boost::add_edge(i, p_angle[i], g_angle);
                }
            }
            correct_normal_orientation(g_angle, in_normals);
        }

        find_components(in_vertices, component_vertices, in_smoothed_v,
            component_smoothed_v, in_normals, component_normals, kdTree, tr_dist, k, isGTNormal, isEuclidean);

        in_vertices.clear();
        in_normals.clear();
        in_smoothed_v.clear();
    }

    recon_timer.end("Initialization");

    for (int component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ..." << std::endl;

        isFaceLoop = true;
        std::vector<Face> faces;
        std::vector<Point> vertices = component_vertices[component_id];
        std::vector<Vector> normals = component_normals[component_id];
        std::vector<Point> smoothed_v = component_smoothed_v[component_id];
        std::string out_component_name = out_name +
            "_component_" + std::to_string(component_id);
        ReconModel.vertices = vertices;

        vector<int> indices(smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(smoothed_v.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(smoothed_v.end(), indices.end())));
        Distance tr_dist;

        recon_timer.end("Kd and DM");

        recon_timer.start("Build MST");

        std::cout << "Init mst" << std::endl;

        // Initial Structure
        m_Graph mst;
        mst.graph = Graph(vertices.size());
        std::vector<m_Edge> full_edges;
        std::vector<m_Edge_length> edge_length;
        std::vector<float> connection_max_length(vertices.size(), 0.);
        std::vector<float> pre_max_length(vertices.size(), 0.);
        mst.isLeastEdge = isLeastEdgeScore;
        mst.isEuclidean = isEuclidean;
        mst.exp_genus = exp_genus;
        {
            s_Graph g;
            s_weightMap weightmap = boost::get(boost::edge_weight, g);
            init_graph(smoothed_v, smoothed_v, normals,
                kdTree, tr_dist, k,
                g, weightmap, isGTNormal,
                isEuclidean, connection_max_length, exp_genus, pre_max_length);

            // Generate MST
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor> p(boost::num_vertices(g));

            //std::cout << weightmap[boost::edge(1405161, 941240, g).first] << std::endl;
            //std::cout << weightmap[boost::edge(1405161, 414857, g).first] << std::endl;
            /*std::cout << weightmap[boost::edge(28879, 28369, g).first] << std::endl;
            std::cout << weightmap[boost::edge(28879, 248370, g).first] << std::endl;*/

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance
                = boost::get(boost::vertex_distance, g);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap
                = boost::get(boost::vertex_index, g);

            boost::prim_minimum_spanning_tree
            (g, *boost::vertices(g).first, &p[0],
                distance, weightmap, indexmap,
                boost::default_dijkstra_visitor());

            build_mst(mst, p, isEuclidean, isGTNormal, smoothed_v, normals);

            // Edge arrays and sort
            if (true) {
                int idx = 0;
                s_Graph::edge_iterator ei, ei_end;
                for (boost::tie(ei, ei_end) = boost::edges(g); ei != ei_end; ei++) {
                    if (weightmap[*ei]>pre_max_length[(*ei).m_source]||
                        weightmap[*ei] > pre_max_length[(*ei).m_target])
                        continue;
                    edge_length.push_back(m_Edge_length
                    (weightmap[*ei], full_edges.size()));
                    full_edges.push_back(m_Edge((*ei).m_source, (*ei).m_target));
                    idx++;
                }
                std::sort(edge_length.begin(), edge_length.end(), edge_comparator);
            }
        }

        recon_timer.end("Build MST");

        // Export MST
        if (isDebug) {
            recon_timer.start("Export MST");
            fs::path out_path = out_root_path / ("MST_" + out_component_name + "_C.obj");
            IO.export_graph(mst, out_path);
            recon_timer.end("Export MST");
        }

        // Initialize face loop label
        mst.etf.reserve(6 * vertices.size() - 11);
        init_face_loop_label(mst);

        // Betti number changes
        std::vector<int> betti_1;

        // Make sure to include all vertices
        if(false){
            for (int i = 0; i < smoothed_v.size(); i++) {
                Vertex id = i;
                if (boost::degree(id, mst.graph) == 1) {
                    Vertex parent = (*mst.graph[id].ordered_neighbors.begin()).v;
                    auto neighbor1 = predecessor(mst, parent, id);
                    auto neighbor2 = successor(mst, parent, id);

                    Vector edge1 = mst.graph[id].coords - mst.graph[neighbor1.v].coords;
                    Vector this_normal = mst.graph[id].normal;
                    Vector neighbor_normal1 = mst.graph[neighbor1.v].normal;
                    float Euclidean_dist1 = norm(edge1);
                    float weight1 = Euclidean_dist1;
                    if (!isEuclidean) {
                        weight1 = cal_proj_dist(edge1, this_normal, neighbor_normal1);
                    }
                    bool isSuccess1 = register_face(mst, id, neighbor1.v, faces, kdTree, tr_dist, weight1);

                    Vector edge2 = mst.graph[id].coords - mst.graph[neighbor2.v].coords;
                    Vector neighbor_normal2 = mst.graph[neighbor2.v].normal;
                    float Euclidean_dist2 = norm(edge2);
                    float weight2 = Euclidean_dist2;
                    if (!isEuclidean) {
                        weight2 = cal_proj_dist(edge2, this_normal, neighbor_normal2);
                    }
                    bool isSuccess2 = register_face(mst, id, neighbor2.v, faces, kdTree, tr_dist, weight2);
                    /*if(!(isSuccess1||isSuccess2)) {
                        std::cout << "Can't find a good connection for leaf node: " << id << std::endl;
                    }*/
                }
            }
        }

        // Vanilla MST imp
        bettiNum_1 = 0;
        betti_1.push_back(bettiNum_1);
        if (true)
        {
            // Edge connection
            for (int i = 0; i < edge_length.size(); i++) {
                if (i % 100000 == 0) {
                    //std::cout << "Step " << i << " / " << edge_length.size() << std::endl;
                    showProgressBar(i / float(edge_length.size()));
                }
                unsigned int edge_idx = edge_length[i].second;
                m_Edge this_edge = full_edges[edge_idx];

                if (boost::edge(this_edge.first, this_edge.second, mst.graph).second)
                    continue;

                bool isValid = Vanilla_check(mst, this_edge, kdTree, tr_dist, edge_length[i].first);

               /* if ((this_edge.first == 30045 && this_edge.second == 69461) ||
                    (this_edge.first == 69461 && this_edge.second == 30045)) {
                    fs::path out_path = out_root_path / ("debug_" + out_name + ".obj");
                    IO.export_obj(mst, out_path, faces);
                    out_path = out_root_path / ("debuggraph_" + out_name + ".obj");
                    IO.export_graph(mst, out_path);
                    std::cout << isValid << std::endl;
                }*/

                if (isValid) {
                    bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces, kdTree, tr_dist, edge_length[i].first);
                    if(isAdded)
                        betti_1.push_back(bettiNum_1);
                }
            }
            showProgressBar(1.0);
            std::cout << std::endl;

            // Output
            if (exp_genus != 0 && isDebug) {
                fs::path out_path = out_root_path / ("Vanilla_beforehandle_" + out_component_name + ".obj");
                IO.export_obj(vertices, mst, out_path, faces);
                out_path = out_root_path / ("VanillaGraph_beforehandle_" + out_component_name + ".obj");
                IO.export_graph(mst, out_path);
            }
        }

        // Create handles & Triangulation
        if (exp_genus != 0) {
            mst.isFinalize = true;
            std::vector<Vertex> connected_handle_root;
            connect_handle(smoothed_v, recon_timer, kdTree, tr_dist, mst, connected_handle_root, betti_1, k, isEuclidean);
            if (isDebug) {
                fs::path out_path = out_root_path / ("handle_" + out_component_name + ".obj");
                IO.export_edges(mst, connected_handle_root, out_path);
            }
            std::cout << bettiNum_1 << std::endl;
            stop_faceloop_check();
            triangulate(faces, mst, kdTree, tr_dist, isFaceLoop, isEuclidean, connection_max_length, connected_handle_root, betti_1);
            //connected_handle_root.clear();
            //triangulate(faces, mst, kdTree, tr_dist, isFaceLoop, isEuclidean, connection_max_length, connected_handle_root, true);
        }

        // Fix non-manifold vertices & unreferenced vertices
        if(false){
            std::vector<bool> isRefered(vertices.size(), false);
            for (auto& face : faces) {
                isRefered[face.ids[0]] = true;
                isRefered[face.ids[1]] = true;
                isRefered[face.ids[2]] = true;
            }
            int idx = 0;
            for (auto refer : isRefered) {
                if (!refer) {
                    if (false) {
                        bool isConnected = false;
                        for (auto& neighbor : mst.graph[idx].ordered_neighbors) {
                            Neighbor neighbor1 = successor(mst, neighbor.v, idx);
                            Neighbor neighbor2 = predecessor(mst, neighbor.v, idx);
                            {
                                Vertex v_u1 = neighbor1.v, v_w = idx, v_u2 = neighbor2.v;
                                if (v_u1 == idx)
                                    continue;
                                if (mst.graph[boost::edge(neighbor.v, v_u1, mst.graph).first].appear_time < 2) {
                                    mst.add_edge(v_u1, v_w);
                                    isRefered[neighbor.v] = true;
                                    isRefered[v_u1] = true;
                                    isRefered[v_w] = true;
                                    std::vector<Vertex> item{ neighbor.v, v_w, v_u1 };
                                    add_face(mst, item, faces);
                                    isConnected = true;
                                    break;
                                }
                                if (mst.graph[boost::edge(neighbor.v, v_u2, mst.graph).first].appear_time < 2) {
                                    mst.add_edge(v_u2, v_w);
                                    isRefered[neighbor.v] = true;
                                    isRefered[v_u2] = true;
                                    isRefered[v_w] = true;
                                    std::vector<Vertex> item{ neighbor.v, v_u2, v_w };
                                    add_face(mst, item, faces);
                                    isConnected = true;
                                    break;
                                }
                            }
                            if (isConnected)
                                break;
                        }
                        if (!isConnected)
                            std::cout << "Vertex " << idx << " is not refered" << std::endl;
                    }
                    else {
                        std::cout << "Vertex " << idx << " is not refered" << std::endl;
                    }
                }
                idx++;
            }
        }

        betti_1.push_back(bettiNum_1);
        // Export betti_1
        //IO.export_betti(betti_1, out_root_path / ("betti_" + out_component_name + ".txt"));

        // Output
        fs::path out_path = out_root_path / ("Vanilla_" + out_component_name + ".obj");
        IO.export_obj(vertices, mst, out_path, faces);
        if (isDebug) {
            out_path = out_root_path / ("VanillaGraph_" + out_component_name + ".obj");
            IO.export_graph(mst, out_path, vertices);
        }
    }
    
    recon_timer.end("algorithm");
    /*int num_sample = 1000000;
    float CD = calculate_chamfer(GTModel, ReconModel, num_sample);
    std::cout << "CD: " << CD << std::endl;*/
    std::string line(40, '=');
    std::cout << line << std::endl << std::endl;
    return;
}

int Reconstructor::reconstruct() {
    // Timer whole
    recon_timer.create("Whole process");
    recon_timer.create("Initialization");
    recon_timer.create("Import obj");
    recon_timer.create("Kd and DM");
    recon_timer.create("Estimate normals");
    recon_timer.create("Build MST");
    recon_timer.create("Remove ambiguous");
    recon_timer.create("Export MST");
    recon_timer.create("Build Rotation System");
    recon_timer.create("Fencing");
    recon_timer.create("Face loop update in Fencing");
    recon_timer.create("Validity check in Fencing");
    recon_timer.create("Triangulation");
    recon_timer.create("algorithm");

    recon_timer.start("Whole process");

    if (mode == "box_sample") {
        int n = 10000;
        fs::path box_model = generate_box_samples(n);
        isGTNormal = true;
        isProgressive = true;
        reconstruct_single();
    }
    else if (mode == "recon_folder") {
        traverse_and_reconstruct(model_path, k);
    }
    else if (mode == "single_file")
        reconstruct_single();
    else if (mode == "dtu") {
        isGTNormal = false;
        isProgressive = true;
        for (int i = 22; i < 129; i++) {
            std::string num = std::to_string(i);
            int n = 3;
            int precision = n - std::min<int>(n, int(num.size()));
            num = std::string(precision, '0').append(num);
            model_name = "stl" + num + "_total.ply";
            fs::path this_path(root_path / model_name);
            out_name = "stl" + num;
            if (fs::exists(this_path)) {
                std::cout << "Start Processing " + model_name + ":D" << std::endl;
                reconstruct_single();
            }
        }
    }
    else if (mode == "noise")
        noise_experiment();
    else
        return -1;
    recon_timer.end("Whole process");
    recon_timer.show();
    return 1;
}

void Reconstructor::traverse_and_reconstruct(const fs::path& dirPath, int k)
{
    if (!fs::exists(dirPath) || !fs::is_directory(dirPath)) {
        std::cerr << "Invalid directory: " << dirPath << std::endl;
        return;
    }

    fs::path tmp_root_path = out_root_path;
    fs::path time_recorder_path = out_root_path / "time.txt";
    fs::ofstream time_file(time_recorder_path, std::ios_base::app);
    long last_time = 0;
    for (fs::directory_iterator it(dirPath); it != fs::directory_iterator(); ++it) {
        const fs::path& filePath = it->path();
        
        if (fs::is_regular_file(filePath)) {
            if (filePath.extension() == ".obj" || filePath.extension()==".ply") {
                model_name = filePath.filename().string();
                out_root_path = tmp_root_path / filePath.stem().string();
                if (!fs::exists(out_root_path))
                    fs::create_directories(out_root_path);
                out_name = filePath.stem().string();
                reconstruct_single();
                long time = recon_timer.log("algorithm");
                time_file << out_name << ": " << time-last_time << std::endl;
                last_time = time;
            }
        }
        
        else if (fs::is_directory(filePath)) {
            traverse_and_reconstruct(filePath, k);
        }
    }
    time_file.close();
    return;
}

void Reconstructor::showProgressBar(float progress) {
    int barWidth = 70;  // Width of the progress bar

    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";  // \r returns to the beginning of the line
    std::cout.flush();  // Flush the output to show the progress
}

int main(int argc, char* argv[]){

    // Check if a parameter was provided
    std::string path;
    if (argc == 2) {
        path = argv[1]; // argv[1] is the first parameter
    }
    else if (argc > 2) {
        std::cout << "Too many arguments!" << std::endl;
        return 1;
    }
    else {
        std::cout << "Please provide the path to the config file as an argument!" << std::endl;
        return 1;
    }

    //fs::path config_path(fs::path("../../configs/tree_config.txt"));
    fs::path config_path(path);
    //std::cout << config_path << std::endl;
    Reconstructor recon(config_path);
    if (!fs::exists(recon.get_model_path())) {
        std::cout << recon.get_model_path() << " does not exist :(" << std::endl;
        return 1;
    }

    int success = recon.reconstruct();
    if (success == -1) {
        std::cout << "Mode Error!" << std::endl;
        return 1;
    }
    return 0;
}