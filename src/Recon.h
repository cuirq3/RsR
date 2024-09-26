#pragma once
#include <iostream>
#include <algorithm>
#include <utility>
#include <set>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <sstream>
#include <CGAL/mst_orient_normals.h>
#include <boost/graph/prim_minimum_spanning_tree.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <chrono>

#include "io.h"
#include "math_util.h"
#include "graph_utils.h"
#include "normal.h"
#include "mst.h"
#include "timer.h"
#include "faceloop.h"
#include "triangulation.h"
#include "RS.h"

using namespace std;

// Noise settings
std::vector<float> sigmas{ 1.0 };
std::vector<float> amplitudes{ .1, .2, .3, .4, .5 };
//std::vector<float> amplitudes{0.1};

// Noise type: random, horizontal, vertical
std::string noise_type = "random";

class Reconstructor {
public:
    // Basic use
	Reconstructor() {
		isNoiseExperiment = false;
		isDTUGeneration = false;
		isEuclidean = true;
		isGTNormal = false;
		isDebug = false;
		isFaceLoop = true;
		k = 30;
		r = 20.;
		n = 50;
		theta = 60.;
		exp_genus = -1;
		set_mode();
	}

	// Initializor for reading config file
	Reconstructor(fs::path in_config_path) {
		isNoiseExperiment = false;
		isDTUGeneration = false;
		isEuclidean = true;
		isGTNormal = false;
		isFaceLoop = true;
		isDebug = false;
		k = 30;
		r = 20.;
		n = 50;
		theta = 60.;
		exp_genus = -1;
		config_path = in_config_path;
		read_config(config_path);
		set_mode();
	}
    
    void read_config(fs::path);
    
	void stop_faceloop_check() {
		isFaceLoop = false;
	}

	fs::path get_model_path() {
		return model_path;
	}

	void reconstruct_single(std::string noise_type = "", float sigma = 0, float amplitude = 0);
	
	int reconstruct();

    void traverse_and_reconstruct(const fs::path& dirPath);

	void noise_experiment();

	void showProgressBar(float progress);

	Timer recon_timer;

private:
	bool isNoiseExperiment;
	bool isDTUGeneration;
	bool isEuclidean;
	bool isGTNormal;
	bool isDebug;
	bool isFaceLoop;
	int k;
	float r;
	float theta;
	int n;
	int exp_genus;
	fs::path model_path;
	fs::path root_path;
	string model_name;
	fs::path out_root_path;
	string out_name;
	fs::path config_path;
	string mode;

	void set_mode() {
		if (model_name == "all")
			mode = "recon_folder";
		else if (isNoiseExperiment)
			mode = "noise";
		else if (isDTUGeneration)
			mode = "dtu";
		else
			mode = "single_file";
	}
};
