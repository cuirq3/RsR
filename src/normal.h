#pragma once
#include <boost/iterator/zip_iterator.hpp>
#include "graph_utils.h"

#ifndef NORMAL_H
#define NORMAL_H

void estimate_normal(const std::vector<Point>&,
	const Tree&, const Distance&, bool,
	std::vector<Vector>&,
	std::vector<int>&, float&);

void replace_zero_normal(const std::vector<int>&,
	const Tree&, const Distance&, const std::vector<Point>&,
	std::vector<Vector>&);

void correct_normal_orientation(s_Graph&, std::vector<Vector>&);

void weighted_smooth(const std::vector<Point>& vertices,
	std::vector<Point>& smoothed_v, const std::vector<Vector>& normals,
	const Tree& kdTree, const Distance& tr_dist,
	float diagonal_length);

#endif //NORMAL_H