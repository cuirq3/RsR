#include "Recon.h"

/**
 * @brief Read the input config file and initialize
 *
 * @param config_path: Path to the config file
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
        if (instruct == "isNoiseExp") {
            isNoiseExperiment = (tokens[1] == "true");
        }
        if (instruct == "k") {
            k = (std::stoi(tokens[1]));
        }
        if (instruct == "genus") {
            exp_genus = (std::stoi(tokens[1]));
        }
        if (instruct == "r") {
            r = (std::stof(tokens[1]));
        }
        if (instruct == "theta") {
            theta = (std::stof(tokens[1]));
        }
        if (instruct == "n") {
            n = (std::stoi(tokens[1]));
        }
        tokens.clear();
    }
}

/**
 * @brief Perform noise experiments by adding noise to vertex positions and normal
 *
 * @param None
 * @return None
 */
void Reconstructor::noise_experiment() {
    fs::path file_path = model_path / model_name;
    isGTNormal = false;
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
    {
        vector<Point> vertices;
        vector<Vector> normals;
        vector<Face> GT_faces;
        IO.read_obj(file_path, vertices, normals, GT_faces);
    }

    // Without noise
    isEuclidean = false;

    for (auto sigma : sigmas) {
        for (auto amplitude : amplitudes) {
            isFaceLoop = true;
            std::cout << "Experimenting sigma " + std::to_string(sigma) + " amplitude " + std::to_string(amplitude) + "..." << std::endl;
            std::string this_out_name = this_model_name + std::to_string(sigma) + "_" + std::to_string(amplitude);
            out_root_path = recon_root;
            out_name = this_out_name;
            reconstruct_single(noise_type, sigma, amplitude);
        }
    }
}

/**
 * @brief Reconstruct a single file
 *
 * @param noise_type: type of noise added for the noise experiments
 * @param sigma: the standard deviation of added noise
 * @param amplitude: the amplitude of added noise
 * 
 * @return None
 */
void Reconstructor::reconstruct_single(std::string noise_type, float sigma, float amplitude) {

    recon_timer.start("Initialization");
    // Init
    io_system IO;
    vector<Point> GT_vertices;
    vector<Face> GT_faces;
    std::vector<Point> in_vertices;
    std::vector<Vector> in_normals;

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
    recon_timer.end("Initialization");

    recon_timer.start("algorithm");
    std::cout << "find components" << std::endl;
    // Find components
    std::vector<std::vector<Point>> component_vertices;
    std::vector<std::vector<Point>> component_smoothed_v;
    std::vector<vector<Vector>> component_normals;
    {
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
            component_smoothed_v, in_normals, component_normals, kdTree,
            tr_dist, k, isEuclidean, theta, r);

        in_vertices.clear();
        in_normals.clear();
        in_smoothed_v.clear();
    }

    for (int component_id = 0; component_id < component_vertices.size(); component_id++) {
        std::cout << "Reconstructing component " + std::to_string(component_id) + " ..." << std::endl;

        isFaceLoop = true;
        std::vector<Face> faces;
        std::vector<Point> vertices = component_vertices[component_id];
        std::vector<Vector> normals = component_normals[component_id];
        std::vector<Point> smoothed_v = component_smoothed_v[component_id];
        std::string out_component_name = out_name +
            "_component_" + std::to_string(component_id);

        vector<int> indices(smoothed_v.size());
        std::iota(indices.begin(), indices.end(), 0);
        // Insert number_of_data_points in the tree
        Tree kdTree(boost::make_zip_iterator(boost::make_tuple(smoothed_v.begin(), indices.begin())),
            boost::make_zip_iterator(boost::make_tuple(smoothed_v.end(), indices.end())));
        Distance tr_dist;

        recon_timer.start("Build MST");

        std::cout << "Init mst" << std::endl;

        // Initial Structure
        m_Graph mst;
        mst.graph = Graph(vertices.size());
        std::vector<m_Edge> full_edges;
        std::vector<m_Edge_length> edge_length;
        std::vector<float> connection_max_length(vertices.size(), 0.);
        std::vector<float> pre_max_length(vertices.size(), 0.);
        mst.isEuclidean = isEuclidean;
        mst.exp_genus = exp_genus;
        {
            s_Graph g;
            s_weightMap weightmap = boost::get(boost::edge_weight, g);
            init_graph(smoothed_v, smoothed_v, normals,
                kdTree, tr_dist, k,
                g, weightmap, isEuclidean, connection_max_length,
                exp_genus, pre_max_length, theta);

            // Generate MST
            std::vector<boost::graph_traits< s_Graph >::vertex_descriptor> p(boost::num_vertices(g));

            boost::property_map< s_Graph, boost::vertex_distance_t >::type distance
                = boost::get(boost::vertex_distance, g);
            boost::property_map< s_Graph, boost::vertex_index_t >::type indexmap
                = boost::get(boost::vertex_index, g);

            boost::prim_minimum_spanning_tree
            (g, *boost::vertices(g).first, &p[0],
                distance, weightmap, indexmap,
                boost::default_dijkstra_visitor());

            build_mst(mst, p, isEuclidean, smoothed_v, normals);

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
            fs::path out_path = out_root_path / ("MST_" + out_component_name + "_C.obj");
            IO.export_graph(mst, out_path);
        }

        // Initialize face loop label
        mst.etf.reserve(6 * vertices.size() - 11);
        init_face_loop_label(mst);

        // Betti number changes
        std::vector<int> betti_1;

        // Vanilla MST imp
        bettiNum_1 = 0;
        betti_1.push_back(bettiNum_1);
        //int inserted_edge = 0;
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

                bool isValid = Vanilla_check(mst, this_edge, kdTree, tr_dist);

                if (isValid) {
                    bool isAdded = register_face(mst, this_edge.first, this_edge.second, faces, kdTree, tr_dist, edge_length[i].first);
                    if(isAdded)
                        betti_1.push_back(bettiNum_1);
                    //Script to make FF video
                    //inserted_edge++;
                    //if (inserted_edge % int(146985 / 96) == 0) {
                    //    fs::path out_path = out_root_path / ("Frame_" + std::to_string(int(inserted_edge/int(146985/96))) + ".obj");
                    //    IO.export_obj(vertices, mst, out_path, faces);
                    //    out_path = out_root_path / ("Graph_Frame_" + std::to_string(int(inserted_edge / int(146985 / 96))) + ".obj");
                    //    IO.export_graph(mst, out_path);
                    //}
                }
            }
            showProgressBar(1.0);
            std::cout << std::endl;
            //std::cout << inserted_edge << std::endl;

            // Output
            if (exp_genus != 0 && isDebug) {
                fs::path out_path = out_root_path / ("Beforehandle_" + out_component_name + ".obj");
                IO.export_obj(vertices, mst, out_path, faces);
                out_path = out_root_path / ("Graph_beforehandle_" + out_component_name + ".obj");
                IO.export_graph(mst, out_path);
            }
        }

        // Create handles & Triangulation
        if (exp_genus != 0) {
            mst.isFinalize = true;
            std::vector<Vertex> connected_handle_root;
            connect_handle(smoothed_v, kdTree, tr_dist, mst, connected_handle_root, betti_1, k, isEuclidean, n);
            if (isDebug) {
                fs::path out_path = out_root_path / ("handle_" + out_component_name + ".obj");
                IO.export_edges(mst, connected_handle_root, out_path);
            }
            stop_faceloop_check();
            triangulate(faces, mst, kdTree, tr_dist, isFaceLoop, isEuclidean, connection_max_length, connected_handle_root, betti_1);
        }

        betti_1.push_back(bettiNum_1);

        // Export betti_1 : Interesting experiment with betti number.
        //IO.export_betti(betti_1, out_root_path / ("betti_" + out_component_name + ".txt"));

        // Output
        fs::path out_path = out_root_path / (out_component_name + ".obj");
        IO.export_obj(vertices, mst, out_path, faces);
        if (isDebug) {
            out_path = out_root_path / ("Graph_" + out_component_name + ".obj");
            IO.export_graph(mst, out_path, vertices);
        }
    }
    
    recon_timer.end("algorithm");
    std::string line(40, '=');
    std::cout << line << std::endl << std::endl;
    return;
}

/**
 * @brief Entrance of reconstruction algorithm, specific timers are set
 *
 * @param None
 * @return if error happens, return -1. Otherwise return 1
 */
int Reconstructor::reconstruct() {
    // Timer whole
    recon_timer.create("Whole process");
    recon_timer.create("Initialization");
    recon_timer.create("Import obj");
    recon_timer.create("Estimate normals");
    recon_timer.create("Build MST");
    recon_timer.create("Build Rotation System");
    recon_timer.create("algorithm");

    recon_timer.start("Whole process");

    if (mode == "recon_folder") {
        traverse_and_reconstruct(model_path);
    }
    else if (mode == "single_file")
        reconstruct_single();
    else if (mode == "dtu") {
        isGTNormal = false;
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

/**
 * @brief Reconstruct a whole folder (currently only detect .obj and .ply files)
 *
 * @param dirPath: path to the directory
 * 
 * @return None
 */
void Reconstructor::traverse_and_reconstruct(const fs::path& dirPath)
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
            traverse_and_reconstruct(filePath);
        }
    }
    time_file.close();
    return;
}

/**
 * @brief Progress bar indicating the reconstruction process
 *
 * @param progress: current progress
 *
 * @return None
 */
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