#pragma once
#include <Eigen/Dense>
#include "graph_utils.h"


struct m_Model {
	std::vector<Point> vertices;
	std::vector<Face> faces;
	std::vector<Vector> normals;
};

float calculate_chamfer(const m_Model& GT_model,const m_Model& Recon_model, int num_sample);